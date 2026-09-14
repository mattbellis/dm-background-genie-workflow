import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

from nudimu import events as E
from nudimu import workflow as W
from nudimu.medium import CRUSTAL_ROCK_MASS_FRACTIONS

from test_dimuon import frame, make_event


def test_parquet_roundtrip_preserves_lists(tmp_path):
    df = frame(
        make_event(event_id=0, fs=[(211, 1.0, 2.0, 100.0), (-321, 0.0, 0.0, 50.0)]),
        make_event(event_id=1, fs=[]),
    )
    path = E.write_events(df, tmp_path / "a.parquet", meta={"e_nu_gev": 1000.0})
    back, meta = E.read_events(path)
    assert meta["stage"] == "genie_events"
    assert meta["n_events"] == 2
    assert meta["e_nu_gev"] == 1000.0
    assert list(back.loc[0, "fs_pdg"]) == [211, -321]
    assert len(back.loc[1, "fs_pdg"]) == 0


def test_footer_key_matches_earthshine_convention(tmp_path):
    df = frame(make_event())
    path = E.write_events(df, tmp_path / "a.parquet")
    meta = pq.read_table(path).schema.metadata
    assert E.FOOTER_KEY in meta
    assert json.loads(meta[E.FOOTER_KEY].decode())["stage"] == "genie_events"


def test_column_pruning(tmp_path):
    df = frame(make_event())
    path = E.write_events(df, tmp_path / "a.parquet")
    back, _ = E.read_events(path, columns=["event_id", "e_nu"])
    assert list(back.columns) == ["event_id", "e_nu"]


def test_audit_final_states_counts():
    df = frame(
        make_event(event_id=0, fs=[(211, 0, 0, 10.0), (211, 0, 0, 20.0)]),
        make_event(event_id=1, fs=[(421, 0, 0, 30.0)]),
    )
    audit = E.audit_final_states(df)
    row = audit[audit.pdg == 211].iloc[0]
    assert row["count"] == 2
    assert row["per_event"] == pytest.approx(1.0)
    assert 421 in set(audit.pdg)


# ------------------------------------------------------------- workflow

def test_target_string_normalises():
    s = W.target_string(CRUSTAL_ROCK_MASS_FRACTIONS)
    fracs = [float(p.split("[")[1].rstrip("]")) for p in s.split(";")]
    assert sum(fracs) == pytest.approx(1.0, abs=1e-5)
    assert len(fracs) == len(CRUSTAL_ROCK_MASS_FRACTIONS)


def test_spline_energy_covers_the_flux():
    cfg = W.GenieConfig()
    assert cfg.spline_max_energy >= max(cfg.energies_gev)


def test_spline_commands_one_per_element_plus_merge():
    cfg = W.GenieConfig(flavours=(14,), sf_warmup=False)
    cmds = W.spline_commands(cfg)
    gmkspl = [c for _, _, c in cmds if c.startswith("gmkspl")]
    gspladd = [c for _, _, c in cmds if c.startswith("gspladd")]
    assert len(gmkspl) == len(cfg.targets)
    assert len(gspladd) == 1
    assert f"-e {cfg.spline_max_energy:g}" in gmkspl[0]
    assert f"--tune {cfg.tune}" in gmkspl[0]


def test_merge_depends_on_every_element():
    cfg = W.GenieConfig(flavours=(14,))
    cmds = W.spline_commands(cfg)
    merged = [c for c in cmds if str(c[0]).endswith("ROCK.xml")][0]
    assert len(merged[1]) == len(cfg.targets)


def test_generation_uses_the_merged_spline_and_a_seed():
    cfg = W.GenieConfig(flavours=(14,), energies_gev=(1000.0,))
    (target, deps, cmd), = W.generation_commands(cfg)
    assert "ROCK.xml" in cmd
    assert "--seed" in cmd
    assert f"--tune {cfg.tune}" in cmd
    assert deps[0].name.endswith("ROCK.xml")


def test_seeds_are_distinct_across_the_matrix():
    cfg = W.GenieConfig()
    seeds = {cfg.seed(p, e) for p in cfg.flavours for e in cfg.energies_gev}
    assert len(seeds) == len(cfg.flavours) * len(cfg.energies_gev)


def test_makefile_is_wellformed(tmp_path):
    cfg = W.GenieConfig(flavours=(14, -14), energies_gev=(1000.0, 10000.0),
                        root=tmp_path / "genie")
    path = W.write_makefile(cfg, tmp_path / "Makefile")
    text = path.read_text()
    assert "all: derived" in text
    assert "gmkspl" in text and "gspladd" in text and "gevgen" in text and "gntpc" in text
    # every recipe line is tab-indented
    for line in text.splitlines():
        if line.startswith("    ") and not line.startswith("\t"):
            raise AssertionError(f"recipe indented with spaces: {line!r}")


def test_makefile_targets_are_unique(tmp_path):
    cfg = W.GenieConfig(root=tmp_path / "genie")
    all_targets = [t for group in (W.spline_commands(cfg),
                                   W.generation_commands(cfg),
                                   W.conversion_commands(cfg))
                   for t, _, _ in group]
    assert len(all_targets) == len(set(map(str, all_targets)))


# ------------------------------------------- HEDIS structure-function warmup

def test_sf_warmup_is_a_prerequisite_of_every_spline():
    cfg = W.GenieConfig(flavours=(14, -14))
    cmds = W.spline_commands(cfg)
    stamp = cfg.root / "splines" / ".hedis_sf_ready"
    gmkspl = [(t, d, c) for t, d, c in cmds if c.startswith("gmkspl")
              and not str(t).endswith(".hedis_sf_ready")]
    assert gmkspl, "no real spline jobs"
    for _, deps, _ in gmkspl:
        assert stamp in deps, "spline job would race the SF cache build"


def test_sf_warmup_is_cheap_and_touches_the_stamp():
    cfg = W.GenieConfig()
    stamp, deps, cmd = W.sf_warmup_command(cfg)
    assert deps == []
    assert "-t 1000010010" in cmd      # free proton, throwaway
    assert "-n 5" in cmd               # few knots
    assert cmd.rstrip().endswith(str(stamp))
    assert f"--tune {cfg.tune}" in cmd


def test_sf_warmup_can_be_disabled():
    cfg = W.GenieConfig(flavours=(14,), sf_warmup=False)
    cmds = W.spline_commands(cfg)
    assert not any(str(t).endswith(".hedis_sf_ready") for t, _, _ in cmds)
    for t, deps, c in cmds:
        if c.startswith("gmkspl"):
            assert deps == []


def test_makefile_exposes_sf_warmup(tmp_path):
    cfg = W.GenieConfig(flavours=(14,), energies_gev=(1000.0,),
                        root=tmp_path / "genie")
    text = W.write_makefile(cfg, tmp_path / "Makefile").read_text()
    assert "sf-warmup:" in text
    assert ".PHONY" in text and "sf-warmup" in text.split(".PHONY")[1].split("\n")[0]
