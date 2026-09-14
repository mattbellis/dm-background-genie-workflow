"""End-to-end analysis driver.

    # convert GENIE output to the canonical event parquet
    python scripts/analyse.py convert genie/derived/*.gst.root --out genie/parquet

    # what did GENIE actually leave in the final state?
    python scripts/analyse.py audit genie/parquet/*.parquet

    # opposite-sign dimuon probability per interaction, vs E_nu
    python scripts/analyse.py dimuon genie/parquet/*.parquet \
        --standoff 50 --out dimuon_vs_energy.parquet
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from nudimu import dimuon as D
from nudimu import events as E
from nudimu import transport as T
from nudimu.medium import STANDARD_ROCK

_STEM = re.compile(r"(?P<flavour>nu\w*?)_E(?P<energy>[0-9.e+]+)_n(?P<n>\d+)_s(?P<seed>\d+)")


def _meta_from_name(path: Path) -> dict:
    m = _STEM.search(path.name)
    if not m:
        return {"source_file": path.name}
    return {
        "source_file": path.name,
        "flavour": m.group("flavour"),
        "e_nu_gev": float(m.group("energy")),
        "n_requested": int(m.group("n")),
        "seed": int(m.group("seed")),
    }


def cmd_convert(args) -> None:
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    for p in map(Path, args.inputs):
        if p.suffix == ".root":
            df = E.read_gst(p, max_events=args.max_events)
        else:
            df = E.read_nuhepmc(p, max_events=args.max_events)
        meta = _meta_from_name(p)
        dest = outdir / (p.name.split(".")[0] + ".parquet")
        E.write_events(df, dest, meta=meta)
        print(f"{p.name}: {len(df)} events -> {dest}")


def _load(paths) -> pd.DataFrame:
    frames = []
    for i, p in enumerate(paths):
        df, meta = E.read_events(p)
        # keep event_id unique across files
        df["event_id"] = df["event_id"].astype("int64") + i * 10_000_000
        if "e_nu_gev" in meta:
            df["e_nu_nominal"] = meta["e_nu_gev"]
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def cmd_audit(args) -> None:
    df = _load(args.inputs)
    audit = E.audit_final_states(df, top=args.top)
    print(audit.to_string(index=False))
    charm = audit[audit.pdg.abs().isin([411, 421, 431, 4122])]
    print()
    if len(charm):
        print("charmed hadrons ARE in the final state: keep include_charm=True")
    else:
        n_mu = audit[audit.pdg.abs() == 13]["per_event"].sum()
        print(f"no charmed hadrons in the final state; muons per event = {n_mu:.4f}")
        print("if that is above 1 for a CC sample, GENIE decayed the charm for")
        print("you and you must run with --no-charm to avoid double counting")


def cmd_dimuon(args) -> None:
    df = _load(args.inputs)
    if args.e_mu_min is not None:
        e_min = args.e_mu_min
    else:
        e_min = float(T.muon_threshold_energy(args.standoff, STANDARD_ROCK))
    print(f"standoff {args.standoff} m  ->  muon threshold {e_min:.1f} GeV")

    cand = D.build_candidates(df, e_mu_min=e_min, include_charm=not args.no_charm)
    per = D.event_probabilities(
        cand, df, max_opening_angle=args.max_angle
    )
    summary = D.probability_vs_energy(per, by=("e_nu", "nu_pdg", "is_cc"))
    print()
    print(summary.to_string(index=False))
    print()
    print("contribution by source species:")
    print(D.source_breakdown(cand[cand.kind == "decay"]).head(12).to_string(index=False))

    if args.out:
        E.write_events(summary, args.out, stage="dimuon_probability",
                       meta={"e_mu_min_gev": e_min, "standoff_m": args.standoff,
                             "include_charm": not args.no_charm,
                             "max_opening_angle_rad": args.max_angle})
        print(f"\nwrote {args.out}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("convert")
    c.add_argument("inputs", nargs="+")
    c.add_argument("--out", default="genie/parquet")
    c.add_argument("--max-events", type=int, default=None)
    c.set_defaults(func=cmd_convert)

    a = sub.add_parser("audit")
    a.add_argument("inputs", nargs="+")
    a.add_argument("--top", type=int, default=25)
    a.set_defaults(func=cmd_audit)

    d = sub.add_parser("dimuon")
    d.add_argument("inputs", nargs="+")
    d.add_argument("--standoff", type=float, default=50.0,
                   help="metres of rock between vertex and detector")
    d.add_argument("--e-mu-min", type=float, default=None,
                   help="override the standoff-derived muon threshold")
    d.add_argument("--max-angle", type=float, default=None,
                   help="max opening angle to the primary lepton, radians")
    d.add_argument("--no-charm", action="store_true")
    d.add_argument("--out", default=None)
    d.set_defaults(func=cmd_dimuon)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
