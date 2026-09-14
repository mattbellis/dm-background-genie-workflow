"""Deterministic generation of the GENIE command lines.

Nothing here runs GENIE.  It emits commands and a Makefile, which is the
whole point: the job list becomes a reviewable artifact, everything is
re-runnable, and the expensive spline step is a proper build dependency
instead of a comment you edit by hand.

Layout produced
---------------
    splines/xsec_<flavour>_<tune>_<target>.xml   one per element, parallel
    splines/xsec_<flavour>_<tune>_ROCK.xml       merged with gspladd
    events/<flavour>_E<energy>_n<N>_s<seed>.ghep.root     source of truth
    derived/<same>.gst.root                                regenerable
    derived/<same>.parquet                                 regenerable

The ghep files are the irreplaceable ones.  gst and parquet are cheap to
rebuild from them, so treat ``derived/`` as disposable exactly as in the
EarthShine signal-side storage layout.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence

from .medium import CRUSTAL_ROCK_MASS_FRACTIONS

__all__ = ["GenieConfig", "FLAVOURS", "spline_commands", "generation_commands",
           "conversion_commands", "write_makefile", "target_string"]

FLAVOURS: Dict[int, str] = {
    12: "nue", -12: "nuebar",
    14: "numu", -14: "numubar",
    16: "nutau", -16: "nutaubar",
}


@dataclass
class GenieConfig:
    """Everything that must be pinned for the output to be reproducible."""

    #: GENIE tune.  For 1-100 TeV you want the high-energy tune, not the
    #: few-GeV default.  G18_* tunes are validated at accelerator energies
    #: and are being extrapolated well past their comfort zone above a few
    #: hundred GeV.  Check `gevgen --help` on your build for the exact
    #: tune/generator-list names available.
    tune: str = "GHE19_00b_00_000"
    event_generator_list: str = "HEDIS"

    #: Spline energy ceiling in GeV.  Must cover the top of the flux you
    #: intend to fold in, with headroom.
    spline_max_energy: float = 1.0e5
    spline_knots: int = 200

    flavours: Sequence[int] = (14, -14, 12, -12, 16, -16)
    energies_gev: Sequence[float] = (
        1.0e3, 2.0e3, 5.0e3, 1.0e4, 2.0e4, 5.0e4, 1.0e5,
    )
    n_events: int = 200_000
    base_seed: int = 20260911

    targets: Dict[int, float] = field(
        default_factory=lambda: dict(CRUSTAL_ROCK_MASS_FRACTIONS)
    )
    root: Path = Path("genie")

    #: Build the HEDIS structure-function cache with a single serial job
    #: before any parallel spline jobs start.  With an empty SF cache,
    #: `make -j8 splines` puts eight processes into the same table files at
    #: once.  Set False only if the tables are already in place.
    sf_warmup: bool = True

    def flavour_name(self, pdg: int) -> str:
        return FLAVOURS[int(pdg)]

    def seed(self, pdg: int, energy: float) -> int:
        return self.base_seed + abs(int(pdg)) * 1000 + int(round(energy)) % 997 + (
            0 if pdg > 0 else 500
        )


def target_string(targets: Dict[int, float]) -> str:
    """GENIE composite-target spec with mass fractions normalised to 1."""
    total = sum(targets.values())
    parts = [f"{pdg}[{frac / total:.6f}]" for pdg, frac in sorted(targets.items())]
    return ";".join(parts)


def _spline_path(cfg: GenieConfig, pdg: int, target: str) -> Path:
    return cfg.root / "splines" / f"xsec_{cfg.flavour_name(pdg)}_{cfg.tune}_{target}.xml"


def merged_spline_path(cfg: GenieConfig, pdg: int) -> Path:
    return cfg.root / "splines" / f"xsec_{cfg.flavour_name(pdg)}_{cfg.tune}_ROCK.xml"


def sf_warmup_command(cfg: GenieConfig) -> tuple:
    """A single cheap gmkspl whose only job is to populate the HEDIS
    structure-function cache.

    The SF tables depend on the tune, not on the target, so one throwaway
    spline on a free proton builds the cache that every later job reads.
    Every real spline target depends on the stamp this produces, so make
    runs it once and serially even under -j.
    """
    stamp = cfg.root / "splines" / ".hedis_sf_ready"
    cmd = " && ".join(
        [
            " ".join(
                [
                    "gmkspl",
                    "-p 14",
                    "-t 1000010010",
                    "-n 5",
                    "-e 1000",
                    f"--tune {cfg.tune}",
                    f'--event-generator-list "{cfg.event_generator_list}"',
                    f"--output-cross-sections {cfg.root}/splines/.sf_warmup.xml",
                ]
            ),
            f"touch {stamp}",
        ]
    )
    return (stamp, [], cmd)


def spline_commands(cfg: GenieConfig) -> List[tuple]:
    """One gmkspl call per (flavour, element), then one gspladd per flavour.

    Per-element is the right granularity: the jobs are independent so they
    parallelise trivially, a failed element costs one job rather than the
    whole set, and adding an element later does not invalidate the rest.
    """
    out: List[tuple] = []
    prereq: List[Path] = []
    if cfg.sf_warmup:
        warm = sf_warmup_command(cfg)
        out.append(warm)
        prereq = [warm[0]]

    for pdg in cfg.flavours:
        pieces = []
        for tgt in sorted(cfg.targets):
            path = _spline_path(cfg, pdg, str(tgt))
            pieces.append(path)
            out.append(
                (
                    path,
                    list(prereq),
                    " ".join(
                        [
                            "gmkspl",
                            f"-p {pdg}",
                            f"-t {tgt}",
                            f"-n {cfg.spline_knots}",
                            f"-e {cfg.spline_max_energy:g}",
                            f"--tune {cfg.tune}",
                            f'--event-generator-list "{cfg.event_generator_list}"',
                            f"--output-cross-sections {path}",
                        ]
                    ),
                )
            )
        merged = merged_spline_path(cfg, pdg)
        out.append(
            (
                merged,
                pieces,
                "gspladd -f " + ",".join(str(p) for p in pieces) + f" -o {merged}",
            )
        )
    return out


def event_stem(cfg: GenieConfig, pdg: int, energy: float) -> str:
    return (
        f"{cfg.flavour_name(pdg)}_E{energy:g}_n{cfg.n_events}"
        f"_s{cfg.seed(pdg, energy)}"
    )


def generation_commands(cfg: GenieConfig) -> List[tuple]:
    out: List[tuple] = []
    tgt = target_string(cfg.targets)
    for pdg in cfg.flavours:
        spline = merged_spline_path(cfg, pdg)
        for energy in cfg.energies_gev:
            stem = event_stem(cfg, pdg, energy)
            outfile = cfg.root / "events" / f"{stem}.ghep.root"
            out.append(
                (
                    outfile,
                    [spline],
                    " ".join(
                        [
                            "gevgen",
                            f"-p {pdg}",
                            f'-t "{tgt}"',
                            f"-n {cfg.n_events}",
                            f"-e {energy:g}",
                            f"--seed {cfg.seed(pdg, energy)}",
                            f"--tune {cfg.tune}",
                            f'--event-generator-list "{cfg.event_generator_list}"',
                            f"--cross-sections {spline}",
                            "--message-thresholds $(GENIE)/config/Messenger_laconic.xml",
                            f"-o {outfile}",
                        ]
                    ),
                )
            )
    return out


def conversion_commands(cfg: GenieConfig) -> List[tuple]:
    """ghep -> gst.  gst is a flat ntuple that uproot reads without ROOT,
    and it already carries every final-state four-vector we need."""
    out: List[tuple] = []
    for pdg in cfg.flavours:
        for energy in cfg.energies_gev:
            stem = event_stem(cfg, pdg, energy)
            ghep = cfg.root / "events" / f"{stem}.ghep.root"
            gst = cfg.root / "derived" / f"{stem}.gst.root"
            out.append((gst, [ghep], f"gntpc -i {ghep} -f gst -o {gst}"))
    return out


_HEADER = """\
# Generated by nudimu.workflow -- do not edit by hand.
# Regenerate with:  python -m nudimu.workflow --out Makefile
#
#   make sf-warmup      build the HEDIS structure-function cache (serial, slow)
#   make splines        build and merge cross-section splines (slow, do once)
#   make events         generate GENIE events
#   make derived        gst conversion
#   make -j8 all        everything, eight jobs at a time
#
# Run `python scripts/check_genie_env.py --tune <tune>` first.  The usual
# first failure is a missing LHAPDF set for the tune's PDF.
#
SHELL := /bin/bash
.DELETE_ON_ERROR:
.SECONDARY:

"""


def write_makefile(cfg: GenieConfig, path="Makefile") -> Path:
    splines = spline_commands(cfg)
    events = generation_commands(cfg)
    derived = conversion_commands(cfg)

    lines = [_HEADER]
    lines.append(".PHONY: all splines events derived sf-warmup clean-derived\n")
    lines.append("all: derived\n")
    if cfg.sf_warmup:
        lines.append(f"sf-warmup: {cfg.root}/splines/.hedis_sf_ready\n")
    lines.append("splines: " + " ".join(str(t) for t, _, _ in splines) + "\n")
    lines.append("events: " + " ".join(str(t) for t, _, _ in events) + "\n")
    lines.append("derived: " + " ".join(str(t) for t, _, _ in derived) + "\n")

    for group in (splines, events, derived):
        for target, deps, cmd in group:
            lines.append(f"{target}: {' '.join(str(d) for d in deps)}")
            lines.append(f"\t@mkdir -p $(dir $@)")
            lines.append(f"\t{cmd}")
            lines.append("")

    lines.append("clean-derived:")
    lines.append(f"\trm -rf {cfg.root}/derived")
    lines.append("")

    path = Path(path)
    path.write_text("\n".join(lines))
    return path


def _main(argv=None):  # pragma: no cover - thin CLI
    import argparse

    ap = argparse.ArgumentParser(description="emit the GENIE build Makefile")
    ap.add_argument("--out", default="Makefile")
    ap.add_argument("--tune", default=GenieConfig.tune)
    ap.add_argument("--generator-list", default=GenieConfig.event_generator_list)
    ap.add_argument("--n-events", type=int, default=GenieConfig.n_events)
    ap.add_argument("--root", default="genie")
    ap.add_argument("--no-sf-warmup", action="store_true",
                    help="skip the serial HEDIS SF cache build")
    args = ap.parse_args(argv)

    cfg = GenieConfig(
        tune=args.tune,
        event_generator_list=args.generator_list,
        n_events=args.n_events,
        root=Path(args.root),
        sf_warmup=not args.no_sf_warmup,
    )
    print(f"wrote {write_makefile(cfg, args.out)}")


if __name__ == "__main__":  # pragma: no cover
    _main()
