import numpy as np
import pytest

from nudimu.medium import STANDARD_ROCK, Medium
from nudimu.species import SPECIES, lookup
from nudimu import transport as T


# --------------------------------------------------------------- muons

def test_muon_range_roundtrip():
    for e in (10.0, 100.0, 1000.0):
        d = T.muon_range_m(e)
        assert T.muon_threshold_energy(d) == pytest.approx(e, rel=1e-9)


def test_muon_range_scale():
    # 10 GeV muon should range roughly 15-25 m in standard rock
    assert 15.0 < T.muon_range_m(10.0) < 25.0
    # 100 GeV roughly 150-200 m
    assert 150.0 < T.muon_range_m(100.0) < 200.0


def test_muon_range_monotonic():
    e = np.logspace(0, 4, 50)
    r = T.muon_range_m(e)
    assert np.all(np.diff(r) > 0)


def test_low_energy_limit_is_ionisation():
    # at E << a/b the range is just E/(a*rho)
    e = 1.0
    expect = e / (STANDARD_ROCK.a_mu * STANDARD_ROCK.density * 1e2)
    assert T.muon_range_m(e) == pytest.approx(expect, rel=1e-2)


# ------------------------------------------------------------- hadrons

def test_pion_decay_length():
    pi = lookup(211)
    # gamma c tau = 55.9 m per GeV
    assert T.decay_length_m(pi, 100.0) == pytest.approx(55.9 * 100.0, rel=2e-3)


def test_kaon_decay_length():
    k = lookup(321)
    assert T.decay_length_m(k, 100.0) == pytest.approx(7.517 * 100.0, rel=2e-3)


def test_kaons_beat_pions():
    """Per particle at fixed energy, kaons are far more likely to decay.

    The ratio is set by the decay lengths (55.9 vs 7.52 m per GeV) times
    the ratio of interaction lengths, since kaons are also absorbed a
    little less readily than pions.
    """
    e = 50.0
    p_pi = T.decay_probability(lookup(211), e)
    p_k = T.decay_probability(lookup(321), e)
    expect = (55.9 / 7.517) * (
        STANDARD_ROCK.lambda_kaon / STANDARD_ROCK.lambda_pion
    )
    assert p_k / p_pi == pytest.approx(expect, rel=0.05)
    assert p_k / p_pi > 5.0


def test_pion_decay_probability_is_small():
    """The whole reason this must be weighted, not sampled."""
    p = T.decay_probability(lookup(211), 20.0)
    assert 1e-5 < p < 1e-3


def test_decay_probability_falls_as_inverse_energy():
    pi = lookup(211)
    assert T.decay_probability(pi, 200.0) == pytest.approx(
        0.1 * T.decay_probability(pi, 20.0), rel=1e-2
    )


def test_charm_is_not_fully_prompt_at_high_energy():
    """D mesons start to be absorbed before decaying above ~1 TeV in rock."""
    d0 = lookup(421)
    assert T.decay_probability(d0, 100.0) > 0.98
    assert 0.7 < T.decay_probability(d0, 1.0e3) < 0.95
    assert 0.25 < T.decay_probability(d0, 1.0e4) < 0.6


def test_denser_rock_suppresses_decay():
    dense = Medium(density=5.0)
    pi = lookup(211)
    assert T.decay_probability(pi, 100.0, dense) < T.decay_probability(pi, 100.0)


def test_longer_interaction_length_enhances_decay():
    soft = STANDARD_ROCK.scaled(2.0)
    pi = lookup(211)
    assert T.decay_probability(pi, 100.0, soft) == pytest.approx(
        2.0 * T.decay_probability(pi, 100.0), rel=1e-2
    )


# ---------------------------------------------------------- muon spectra

def test_pion_muon_fraction_floor():
    """A pion always gives its muon at least (m_mu/m_pi)^2 = 0.573 of E."""
    pi = lookup(211)
    assert T.muon_fraction_above(pi, 100.0, 50.0) == pytest.approx(1.0)
    assert T.muon_fraction_above(pi, 100.0, 57.0) == pytest.approx(1.0)
    assert T.muon_fraction_above(pi, 100.0, 101.0) == pytest.approx(0.0)


def test_pion_muon_fraction_linear_in_between():
    pi = lookup(211)
    zmin = pi.z_min
    # threshold halfway up the allowed band -> half the decays pass
    thr = 0.5 * (zmin + 1.0) * 100.0
    assert T.muon_fraction_above(pi, 100.0, thr) == pytest.approx(0.5, rel=1e-6)


def test_kaon_reaches_lower_fractions():
    k = lookup(321)
    # kaon muon can be as soft as 4.6% of E_K, so at a 10% threshold some
    # decays already fail
    assert T.muon_fraction_above(k, 100.0, 10.0) < 1.0
    assert T.muon_fraction_above(k, 100.0, 4.0) == pytest.approx(1.0)


def test_yield_probability_is_product():
    k = lookup(321)
    e, thr = 200.0, 20.0
    expect = (
        T.decay_probability(k, e) * k.br_mu * T.muon_fraction_above(k, e, thr)
    )
    assert T.muon_yield_probability(k, e, thr) == pytest.approx(expect)


def test_yield_vanishes_below_threshold():
    pi = lookup(211)
    assert T.muon_yield_probability(pi, 10.0, 50.0) == 0.0


def test_all_species_have_sane_table_entries():
    for pdg, sp in SPECIES.items():
        assert sp.mass > 0
        assert sp.ctau_m > 0
        assert 0.0 < sp.br_mu <= 1.0
        assert 0.0 <= sp.z_min < sp.z_max <= 1.0
        assert sp.lambda_key in ("pion", "kaon", "charm", "nucleon")


def test_vectorised():
    pi = lookup(211)
    e = np.array([10.0, 100.0, 1000.0])
    p = T.muon_yield_probability(pi, e, 5.0)
    assert p.shape == (3,)
    assert np.all(np.diff(p) < 0)  # falls with energy
