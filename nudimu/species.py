"""Hadron species that can feed a secondary muon.

For each species we need four things:

1. mass and c*tau, which set the decay length gamma*c*tau = (E/m) c tau;
2. an interaction-length class, which sets the competing absorption length;
3. the inclusive branching fraction to a muon;
4. the lab-frame muon energy fraction z = E_mu / E_parent.

For the two-body decays pi -> mu nu and K -> mu nu, an isotropically
decaying ultra-relativistic parent gives a spectrum that is *exactly flat*
in z between z_min = (m_mu/m_h)^2 and z_max = 1.  For the semileptonic
three-body modes (K_L -> pi mu nu, D -> mu X) a flat spectrum over the
kinematic range is only an approximation; those entries are marked
``approximate_spectrum=True`` and are a known ~tens-of-percent modelling
uncertainty on the charm branch.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

__all__ = ["Species", "SPECIES", "lookup", "M_MU"]

M_MU = 0.1056584  # GeV


@dataclass(frozen=True)
class Species:
    pdg: int
    name: str
    mass: float          # GeV
    ctau_m: float        # metres
    charge: int
    br_mu: float         # inclusive BR to at least one muon
    z_min: float         # min lab muon energy fraction
    z_max: float         # max lab muon energy fraction
    lambda_key: str      # which Medium.lambda_* applies
    approximate_spectrum: bool = False

    @property
    def is_prompt_scale(self) -> bool:
        """True for charmed hadrons (ctau of order 100 um)."""
        return self.ctau_m < 1e-2


def _two_body_zmin(m_parent: float) -> float:
    return (M_MU / m_parent) ** 2


_PI = 0.1395704
_KCH = 0.493677
_KL = 0.497611

SPECIES: Dict[int, Species] = {}


def _add(s: Species) -> None:
    SPECIES[s.pdg] = s


# --- charged pions: pi -> mu nu, BR ~ 100%, flat z in [(m_mu/m_pi)^2, 1] ---
for _pdg, _q in ((211, +1), (-211, -1)):
    _add(
        Species(
            pdg=_pdg,
            name=("pi+" if _q > 0 else "pi-"),
            mass=_PI,
            ctau_m=7.8045,
            charge=_q,
            br_mu=0.9998770,
            z_min=_two_body_zmin(_PI),
            z_max=1.0,
            lambda_key="pion",
        )
    )

# --- charged kaons: Kmu2 (63.56%) dominates; Kmu3 (3.35%) folded in with a
#     softer spectrum is a small correction, so we use Kmu2 kinematics for
#     the whole 66.9% and flag the spectrum as approximate. ---
for _pdg, _q in ((321, +1), (-321, -1)):
    _add(
        Species(
            pdg=_pdg,
            name=("K+" if _q > 0 else "K-"),
            mass=_KCH,
            ctau_m=3.711,
            charge=_q,
            br_mu=0.6356 + 0.0335,
            z_min=_two_body_zmin(_KCH),
            z_max=1.0,
            lambda_key="kaon",
            approximate_spectrum=True,
        )
    )

# --- K0_L -> pi mu nu, 27.04% summed over charges. Three-body; flat
#     approximation over the kinematic z range. ---
_add(
    Species(
        pdg=130,
        name="K0L",
        mass=_KL,
        ctau_m=15.34,
        charge=0,
        br_mu=0.2704,
        z_min=0.0492,
        z_max=0.9172,
        lambda_key="kaon",
        approximate_spectrum=True,
    )
)

# --- charmed hadrons. Semileptonic; the muon carries <z> ~ 0.3.  Modelled
#     as flat in [0, 0.6].  This is the crudest piece of the chain and is
#     the first thing to replace with a proper decay table if the charm
#     branch turns out to dominate the final number. ---
_CHARM = (
    (411, "D+", 1.86966, 311.8e-6, +1, 0.176),
    (-411, "D-", 1.86966, 311.8e-6, -1, 0.176),
    (421, "D0", 1.86484, 122.9e-6, 0, 0.0680),
    (-421, "D0bar", 1.86484, 122.9e-6, 0, 0.0680),
    (431, "Ds+", 1.96835, 151.2e-6, +1, 0.065),
    (-431, "Ds-", 1.96835, 151.2e-6, -1, 0.065),
    (4122, "Lambda_c+", 2.28646, 60.7e-6, +1, 0.036),
    (-4122, "Lambda_c-", 2.28646, 60.7e-6, -1, 0.036),
)
for _pdg, _name, _m, _ct, _q, _br in _CHARM:
    _add(
        Species(
            pdg=_pdg,
            name=_name,
            mass=_m,
            ctau_m=_ct,
            charge=_q,
            br_mu=_br,
            z_min=0.0,
            z_max=0.6,
            lambda_key="charm",
            approximate_spectrum=True,
        )
    )


def lookup(pdg: int) -> Optional[Species]:
    """Return the Species for ``pdg`` or None if it cannot make a muon."""
    return SPECIES.get(int(pdg))


#: Neutral hadrons still produce muons (K0L, D0) but produce no *charge*
#: information for the muon.  We treat the resulting muon sign as 50/50.
NEUTRAL_PARENTS = frozenset(s.pdg for s in SPECIES.values() if s.charge == 0)
