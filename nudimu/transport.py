"""Transport of muons and of muon-producing hadrons through rock.

Two independent pieces:

``muon_range_gcm2`` / ``muon_threshold_energy``
    Continuous-slowing-down range for  -dE/dX = a + b E.  Same functional
    form as Eq. (32) of Feng, Smolinsky & Tanedo, just solved in both
    directions.  If you already trust ``nubkg.muon_transport`` use that
    instead and delete this half of the module; the formula is identical
    and having two copies of it is how they drift apart.

``decay_probability`` / ``muon_yield_probability``
    A hadron produced in the shower either decays or is absorbed.  With
    constant decay length L_d = (E/m) c tau and constant interaction
    length L_i, the probability that decay wins is

        P_decay = L_i / (L_i + L_d) .

    For pions in rock L_i ~ 0.5 m and L_d ~ 56 m x (E/GeV), so P_decay is
    of order 1e-2 / (E/GeV): tiny, which is exactly why this has to be
    computed as a weight rather than sampled.  Sampling it event by event
    would need ~1e8 GENIE events to get percent-level statistics on a 1e-5
    probability; weighting gets the same answer from 1e4 events.
"""

from __future__ import annotations

from typing import Union

import numpy as np

from .medium import STANDARD_ROCK, Medium
from .species import Species

__all__ = [
    "muon_range_gcm2",
    "muon_range_m",
    "muon_threshold_energy",
    "decay_length_m",
    "decay_probability",
    "muon_fraction_above",
    "muon_yield_probability",
]

ArrayLike = Union[float, np.ndarray]


# --------------------------------------------------------------------------
# muons
# --------------------------------------------------------------------------

def muon_range_gcm2(energy: ArrayLike, medium: Medium = STANDARD_ROCK,
                    e_threshold: float = 0.0) -> ArrayLike:
    """Column density a muon of ``energy`` GeV crosses before dropping to
    ``e_threshold``."""
    a, b = medium.a_mu, medium.b_mu
    e = np.asarray(energy, dtype=float)
    return np.log((a + b * e) / (a + b * e_threshold)) / b


def muon_range_m(energy: ArrayLike, medium: Medium = STANDARD_ROCK,
                 e_threshold: float = 0.0) -> ArrayLike:
    """Range in metres."""
    return muon_range_gcm2(energy, medium, e_threshold) / medium.density * 1e-2


def muon_threshold_energy(distance_m: ArrayLike, medium: Medium = STANDARD_ROCK,
                          e_arrival: float = 0.0) -> ArrayLike:
    """Energy a muon needs at creation to still have ``e_arrival`` GeV after
    ``distance_m`` metres of rock.  Inverse of :func:`muon_range_m`."""
    a, b = medium.a_mu, medium.b_mu
    x = medium.grammage_for_distance_m(np.asarray(distance_m, dtype=float))
    return (a / b) * (np.exp(b * x) - 1.0) + e_arrival * np.exp(b * x)


# --------------------------------------------------------------------------
# hadrons
# --------------------------------------------------------------------------

def decay_length_m(species: Species, energy: ArrayLike) -> ArrayLike:
    """Lab-frame decay length gamma*beta*c*tau in metres."""
    e = np.asarray(energy, dtype=float)
    gamma_beta = np.sqrt(np.maximum(e**2 - species.mass**2, 0.0)) / species.mass
    return gamma_beta * species.ctau_m


def decay_probability(species: Species, energy: ArrayLike,
                      medium: Medium = STANDARD_ROCK) -> ArrayLike:
    """Probability the hadron decays before it interacts inelastically."""
    l_dec = decay_length_m(species, energy)
    l_int = medium.interaction_length_m(species.lambda_key)
    return l_int / (l_int + l_dec)


def muon_fraction_above(species: Species, e_parent: ArrayLike,
                        e_mu_min: ArrayLike) -> ArrayLike:
    """Fraction of decays whose muon exceeds ``e_mu_min``.

    Assumes the lab muon energy is flat in z over [z_min, z_max].
    """
    e_parent = np.asarray(e_parent, dtype=float)
    e_mu_min = np.asarray(e_mu_min, dtype=float)
    lo = species.z_min * e_parent
    hi = species.z_max * e_parent
    span = np.where(hi > lo, hi - lo, np.inf)
    frac = (hi - np.maximum(e_mu_min, lo)) / span
    return np.clip(frac, 0.0, 1.0)


def muon_yield_probability(species: Species, e_parent: ArrayLike,
                           e_mu_min: ArrayLike,
                           medium: Medium = STANDARD_ROCK) -> ArrayLike:
    """P(this hadron yields a muon above ``e_mu_min``).

    Product of three factors: decay wins over absorption, the decay is
    muonic, and the muon lands above threshold.
    """
    return (
        decay_probability(species, e_parent, medium)
        * species.br_mu
        * muon_fraction_above(species, e_parent, e_mu_min)
    )


def mean_muon_energy(species: Species, e_parent: ArrayLike,
                     e_mu_min: ArrayLike) -> ArrayLike:
    """Mean muon energy given that it exceeded ``e_mu_min`` (flat spectrum)."""
    e_parent = np.asarray(e_parent, dtype=float)
    lo = np.maximum(np.asarray(e_mu_min, dtype=float), species.z_min * e_parent)
    hi = species.z_max * e_parent
    return np.where(hi > lo, 0.5 * (lo + hi), np.nan)
