"""Readers that turn GENIE output into one canonical event table.

Canonical schema, one row per neutrino interaction:

    event_id     int32
    nu_pdg       int32
    e_nu         float64   GeV
    target_pdg   int32
    is_cc        bool
    mode         str       qel | res | dis | coh | mec | other
    W, Q2, x, y  float64
    lep_pdg      int32     final-state primary lepton
    lep_e, lep_px, lep_py, lep_pz   float64
    fs_pdg       list<int32>    all final-state particles
    fs_e, fs_px, fs_py, fs_pz      list<float64>

The final-state list is *whatever GENIE left undecayed*.  Which particles
that is depends on the decayer configuration of your tune, so
:func:`audit_final_states` is provided to tell you, from the data itself,
whether charmed hadrons are surviving to the final state or have already
been decayed for you.  Get that answer before trusting any charm numbers.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

__all__ = [
    "EVENT_COLUMNS",
    "read_gst",
    "read_nuhepmc",
    "write_events",
    "read_events",
    "audit_final_states",
]

EVENT_COLUMNS = [
    "event_id", "nu_pdg", "e_nu", "target_pdg", "is_cc", "mode",
    "W", "Q2", "x", "y",
    "lep_pdg", "lep_e", "lep_px", "lep_py", "lep_pz",
    "fs_pdg", "fs_e", "fs_px", "fs_py", "fs_pz",
]

FOOTER_KEY = b"earthshine"


# --------------------------------------------------------------------------
# gst  (gntpc -f gst): flat ntuple, readable with uproot, no ROOT needed
# --------------------------------------------------------------------------

_GST_MODE_FLAGS = ["qel", "res", "dis", "coh", "mec", "dfr", "imd", "nuel"]


def read_gst(path, tree: str = "gst", max_events: Optional[int] = None) -> pd.DataFrame:
    """Read a GENIE ``gst`` ntuple into the canonical schema."""
    import uproot

    with uproot.open(f"{path}:{tree}") as t:
        available = set(t.keys())
        wanted = [
            "neu", "Ev", "tgt", "cc", "nc", "W", "Q2", "x", "y",
            "pdgf", "Ef", "pxf", "pyf", "pzf",
            "El", "pxl", "pyl", "pzl",
        ] + _GST_MODE_FLAGS
        wanted = [b for b in wanted if b in available]
        arrs = t.arrays(wanted, library="np", entry_stop=max_events)

    n = len(arrs["Ev"])

    mode = np.full(n, "other", dtype=object)
    for flag in reversed(_GST_MODE_FLAGS):
        if flag in arrs:
            mode[np.asarray(arrs[flag], dtype=bool)] = flag

    nu_pdg = np.asarray(arrs["neu"], dtype=np.int32)
    is_cc = np.asarray(arrs["cc"], dtype=bool)
    # gst has no explicit final-state primary lepton PDG; reconstruct it.
    lep_pdg = np.where(is_cc, np.sign(nu_pdg) * (np.abs(nu_pdg) - 1), nu_pdg)

    df = pd.DataFrame(
        {
            "event_id": np.arange(n, dtype=np.int32),
            "nu_pdg": nu_pdg,
            "e_nu": np.asarray(arrs["Ev"], dtype=float),
            "target_pdg": np.asarray(arrs.get("tgt", np.zeros(n)), dtype=np.int32),
            "is_cc": is_cc,
            "mode": mode,
            "W": np.asarray(arrs.get("W", np.full(n, np.nan)), dtype=float),
            "Q2": np.asarray(arrs.get("Q2", np.full(n, np.nan)), dtype=float),
            "x": np.asarray(arrs.get("x", np.full(n, np.nan)), dtype=float),
            "y": np.asarray(arrs.get("y", np.full(n, np.nan)), dtype=float),
            "lep_pdg": lep_pdg.astype(np.int32),
            "lep_e": np.asarray(arrs["El"], dtype=float),
            "lep_px": np.asarray(arrs["pxl"], dtype=float),
            "lep_py": np.asarray(arrs["pyl"], dtype=float),
            "lep_pz": np.asarray(arrs["pzl"], dtype=float),
        }
    )
    df["fs_pdg"] = [np.asarray(v, dtype=np.int32) for v in arrs["pdgf"]]
    for src, dst in (("Ef", "fs_e"), ("pxf", "fs_px"), ("pyf", "fs_py"), ("pzf", "fs_pz")):
        df[dst] = [np.asarray(v, dtype=float) for v in arrs[src]]
    return df[EVENT_COLUMNS]


# --------------------------------------------------------------------------
# NuHepMC / HepMC3
# --------------------------------------------------------------------------

def read_nuhepmc(path, max_events: Optional[int] = None) -> pd.DataFrame:
    """Read a NuHepMC (HepMC3 ascii) file into the canonical schema.

    Status-code conventions follow NuHepMC: 1 = undecayed physical particle,
    4 = incoming beam, 11 = target.  Some writers use 4 for the target too,
    so the target is identified by nuclear PDG code as a fallback.
    """
    import pyhepmc

    rows: List[Dict] = []
    with pyhepmc.open(str(path)) as reader:
        for i, evt in enumerate(reader):
            if max_events is not None and i >= max_events:
                break
            beam = None
            target = None
            fs = []
            for p in evt.particles:
                pdg, status = int(p.pid), int(p.status)
                if status == 1:
                    fs.append(p)
                elif status == 4 and abs(pdg) in (12, 14, 16):
                    beam = p
                elif status in (11, 4, 20) and abs(pdg) > 1000000000:
                    target = p

            if beam is None:
                raise ValueError(f"no beam neutrino in event {i} of {path}")

            lept = _primary_lepton(fs, int(beam.pid))
            attrs = _hepmc_attributes(evt)
            rows.append(
                {
                    "event_id": np.int32(i),
                    "nu_pdg": np.int32(beam.pid),
                    "e_nu": float(beam.momentum.e),
                    "target_pdg": np.int32(target.pid if target is not None else 0),
                    "is_cc": bool(lept is not None and abs(int(lept.pid)) in (11, 13, 15)),
                    "mode": attrs.get("mode", "other"),
                    "W": attrs.get("W", np.nan),
                    "Q2": attrs.get("Q2", np.nan),
                    "x": attrs.get("x", np.nan),
                    "y": attrs.get("y", np.nan),
                    "lep_pdg": np.int32(lept.pid if lept is not None else 0),
                    "lep_e": float(lept.momentum.e) if lept is not None else np.nan,
                    "lep_px": float(lept.momentum.px) if lept is not None else np.nan,
                    "lep_py": float(lept.momentum.py) if lept is not None else np.nan,
                    "lep_pz": float(lept.momentum.pz) if lept is not None else np.nan,
                    "fs_pdg": np.array([int(p.pid) for p in fs], dtype=np.int32),
                    "fs_e": np.array([float(p.momentum.e) for p in fs]),
                    "fs_px": np.array([float(p.momentum.px) for p in fs]),
                    "fs_py": np.array([float(p.momentum.py) for p in fs]),
                    "fs_pz": np.array([float(p.momentum.pz) for p in fs]),
                }
            )
    return pd.DataFrame(rows, columns=EVENT_COLUMNS)


def _primary_lepton(fs, nu_pdg: int):
    """Highest-energy charged lepton of the matching flavour, else the
    highest-energy outgoing neutrino of the beam flavour (NC)."""
    charged = abs(nu_pdg) - 1
    cands = [p for p in fs if abs(int(p.pid)) == charged]
    if not cands:
        cands = [p for p in fs if int(p.pid) == nu_pdg]
    if not cands:
        return None
    return max(cands, key=lambda p: p.momentum.e)


def _hepmc_attributes(evt) -> Dict:
    out: Dict = {}
    try:
        names = evt.attribute_names()
    except Exception:  # pragma: no cover - depends on pyhepmc version
        return out
    for key in names:
        val = evt.attribute(key)
        if val is None:
            continue
        low = key.lower()
        for target in ("w", "q2", "x", "y"):
            if low.endswith(target):
                try:
                    out[target.upper() if target in ("w", "q2") else target] = float(val)
                except Exception:
                    pass
    return out


# --------------------------------------------------------------------------
# parquet I/O with self-describing footer metadata
# --------------------------------------------------------------------------

def write_events(df: pd.DataFrame, path, meta: Optional[Dict] = None,
                 stage: str = "genie_events") -> Path:
    """Write the event table with an ``earthshine`` footer block."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(df, preserve_index=False)
    payload = dict(meta or {})
    payload.setdefault("stage", stage)
    payload.setdefault("n_events", int(len(df)))
    existing = table.schema.metadata or {}
    table = table.replace_schema_metadata(
        {**existing, FOOTER_KEY: json.dumps(payload).encode()}
    )
    pq.write_table(table, path)
    return path


def read_events(path, columns: Optional[Iterable[str]] = None):
    """Read back an event table.  Returns ``(dataframe, metadata_dict)``."""
    table = pq.read_table(path, columns=list(columns) if columns else None)
    meta = {}
    raw = (table.schema.metadata or {}).get(FOOTER_KEY)
    if raw:
        meta = json.loads(raw.decode())
    return table.to_pandas(), meta


# --------------------------------------------------------------------------
# sanity check on what GENIE actually left in the record
# --------------------------------------------------------------------------

def audit_final_states(df: pd.DataFrame, top: int = 25) -> pd.DataFrame:
    """Count final-state PDG codes.

    Run this first.  If charmed hadrons (411, 421, 431, 4122) appear, GENIE
    left them for you to decay and the charm branch is live.  If they are
    absent but the muon multiplicity per CC event already exceeds one,
    GENIE decayed them itself and you must *not* apply the charm branch on
    top, or you will double count.
    """
    counts: Dict[int, int] = {}
    for arr in df["fs_pdg"]:
        for pdg in np.asarray(arr):
            counts[int(pdg)] = counts.get(int(pdg), 0) + 1
    out = pd.DataFrame(
        sorted(counts.items(), key=lambda kv: -kv[1])[:top],
        columns=["pdg", "count"],
    )
    out["per_event"] = out["count"] / max(len(df), 1)
    return out
