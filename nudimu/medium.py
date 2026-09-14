"""Properties of the rock surrounding the detector.

Everything that depends on "what the muons and hadrons are moving through"
lives here so that it can be varied in one place for systematics.

Interaction lengths are inelastic *nuclear* interaction lengths in g/cm^2.
The values below are scaled from the usual air/water numbers to A ~ 22
standard rock; they carry a real uncertainty of order 10-20% and are the
single biggest knob on the pi/K branch of the answer.  Vary them.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict

__all__ = ["Medium", "STANDARD_ROCK", "CRUSTAL_ROCK_MASS_FRACTIONS"]


# Continental-crust mass fractions, as used in the existing gevgen target
# string.  These are renormalised to 1 on construction of the target spec;
# the raw numbers sum to 0.985.
CRUSTAL_ROCK_MASS_FRACTIONS: Dict[int, float] = {
    1000080160: 0.466,   # O-16
    1000140280: 0.277,   # Si-28
    1000130270: 0.081,   # Al-27
    1000260560: 0.050,   # Fe-56
    1000200400: 0.036,   # Ca-40
    1000110230: 0.028,   # Na-23
    1000190390: 0.026,   # K-39
    1000120240: 0.021,   # Mg-24
}


@dataclass(frozen=True)
class Medium:
    """A homogeneous medium.

    Parameters
    ----------
    density
        g/cm^3.
    lambda_nucleon, lambda_pion, lambda_kaon, lambda_charm
        Inelastic interaction lengths in g/cm^2.  Pions and kaons have
        smaller inelastic cross sections on nuclei than nucleons do, hence
        the longer lengths.  ``lambda_charm`` is used for D/Ds/Lambda_c and
        only matters above ~1 TeV, where charmed hadrons start to interact
        before they decay.
    a_mu, b_mu
        Muon energy loss  -dE/dX = a + b E  with X in g/cm^2, a in
        GeV cm^2/g and b in cm^2/g.
    """

    name: str = "standard_rock"
    density: float = 2.65
    Z: float = 11.0
    A: float = 22.0

    lambda_nucleon: float = 100.0
    lambda_pion: float = 130.0
    lambda_kaon: float = 150.0
    lambda_charm: float = 130.0

    a_mu: float = 2.0e-3
    b_mu: float = 4.4e-6

    def interaction_length_m(self, key: str) -> float:
        """Interaction length in metres for a species class key."""
        try:
            grammage = getattr(self, f"lambda_{key}")
        except AttributeError as exc:  # pragma: no cover - programmer error
            raise KeyError(f"unknown interaction-length class {key!r}") from exc
        return grammage / self.density * 1e-2  # g/cm^2 -> cm -> m

    def grammage_for_distance_m(self, distance_m: float) -> float:
        """Column density in g/cm^2 traversed over ``distance_m`` metres."""
        return distance_m * 1e2 * self.density

    def scaled(self, factor: float) -> "Medium":
        """Return a copy with all hadronic interaction lengths scaled.

        Convenience for the dominant systematic: ``rock.scaled(1.2)``.
        """
        return replace(
            self,
            name=f"{self.name}_lambda_x{factor:g}",
            lambda_nucleon=self.lambda_nucleon * factor,
            lambda_pion=self.lambda_pion * factor,
            lambda_kaon=self.lambda_kaon * factor,
            lambda_charm=self.lambda_charm * factor,
        )


STANDARD_ROCK = Medium()
