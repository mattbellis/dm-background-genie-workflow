"""Print the physics scales that govern the answer, before any GENIE runs.

    python scripts/scales.py

Useful as a first sanity check and as a reference when choosing the muon
energy threshold and the standoff distance.
"""

from __future__ import annotations

import numpy as np

from nudimu.medium import STANDARD_ROCK as ROCK
from nudimu.species import lookup
from nudimu import transport as T


def main() -> None:
    print(f"medium: {ROCK.name}, rho = {ROCK.density} g/cm^3")
    for key in ("nucleon", "pion", "kaon", "charm"):
        print(f"  lambda_{key:<8s} = {getattr(ROCK, 'lambda_' + key):6.1f} g/cm^2"
              f"  = {ROCK.interaction_length_m(key):.3f} m")
    print()

    print("muon range in rock and the energy needed to cross a given standoff")
    print(f"  {'E [GeV]':>9s} {'range [m]':>11s}   | {'standoff [m]':>13s} {'E_min [GeV]':>12s}")
    for e, d in zip([1, 3, 10, 30, 100, 300, 1000, 3000],
                    [1, 3, 10, 30, 100, 300, 1000, 3000]):
        print(f"  {e:9d} {T.muon_range_m(float(e)):11.1f}   |"
              f" {d:13d} {float(T.muon_threshold_energy(float(d))):12.1f}")
    print()

    print("probability a hadron decays before it is absorbed, and the")
    print("probability it delivers a muon above 10 GeV")
    hdr = f"  {'E [GeV]':>9s}"
    for name in ("pi+", "K+", "K0L", "D+", "D0"):
        hdr += f" {name:>20s}"
    print(hdr)
    for e in (10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0):
        line = f"  {e:9.0f}"
        for pdg in (211, 321, 130, 411, 421):
            sp = lookup(pdg)
            pd_ = float(T.decay_probability(sp, e))
            py = float(T.muon_yield_probability(sp, e, 10.0))
            line += f"  {pd_:8.2e}/{py:8.2e}"
        print(line)
    print()

    print("reading: charm decays essentially always up to ~1 TeV, then starts")
    print("being absorbed in the rock first -- the 'prompt' approximation that")
    print("holds in air does not hold in rock above a TeV.")
    print()

    print("per-particle muon yield above 10 GeV, relative to a pion of the")
    print("same energy:")
    for e in (30.0, 300.0, 3000.0):
        ref = float(T.muon_yield_probability(lookup(211), e, 10.0))
        row = [f"E = {e:6.0f} GeV:"]
        for pdg, name in ((321, "K+"), (130, "K0L"), (411, "D+"), (421, "D0")):
            r = float(T.muon_yield_probability(lookup(pdg), e, 10.0)) / ref
            row.append(f"{name} x{r:9.1f}")
        print("  " + "  ".join(row))


if __name__ == "__main__":
    main()
