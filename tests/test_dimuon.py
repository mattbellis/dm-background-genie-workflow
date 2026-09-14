import numpy as np
import pandas as pd
import pytest

from nudimu import dimuon as D
from nudimu.events import EVENT_COLUMNS


def make_event(event_id=0, nu_pdg=14, e_nu=1000.0, is_cc=True,
               lep_pdg=13, lep_p=(0.0, 0.0, 500.0), fs=()):
    """fs: iterable of (pdg, px, py, pz)."""
    lep_e = float(np.linalg.norm(lep_p))
    fs = list(fs)
    row = {
        "event_id": np.int32(event_id), "nu_pdg": np.int32(nu_pdg),
        "e_nu": e_nu, "target_pdg": np.int32(1000140280),
        "is_cc": is_cc, "mode": "dis",
        "W": 10.0, "Q2": 5.0, "x": 0.1, "y": 0.5,
        "lep_pdg": np.int32(lep_pdg), "lep_e": lep_e,
        "lep_px": lep_p[0], "lep_py": lep_p[1], "lep_pz": lep_p[2],
        "fs_pdg": np.array([p[0] for p in fs], dtype=np.int32),
        "fs_e": np.array([float(np.linalg.norm(p[1:])) for p in fs]),
        "fs_px": np.array([p[1] for p in fs]),
        "fs_py": np.array([p[2] for p in fs]),
        "fs_pz": np.array([p[3] for p in fs]),
    }
    return row


def frame(*rows):
    return pd.DataFrame(list(rows), columns=EVENT_COLUMNS)


# -------------------------------------------------------- candidate table

def test_numu_cc_with_no_hadrons_has_no_dimuon():
    ev = frame(make_event())
    c = D.build_candidates(ev, e_mu_min=10.0)
    p = D.event_probabilities(c, ev)
    assert p.loc[0, "p_minus"] == pytest.approx(1.0)   # the primary mu-
    assert p.loc[0, "p_plus"] == 0.0
    assert p.loc[0, "p_os"] == 0.0


def test_pi_plus_supplies_the_opposite_sign_muon():
    ev = frame(make_event(fs=[(211, 0.0, 0.0, 100.0)]))
    c = D.build_candidates(ev, e_mu_min=10.0)
    p = D.event_probabilities(c, ev)
    assert 0.0 < p.loc[0, "p_plus"] < 1e-3
    assert p.loc[0, "p_os"] == pytest.approx(p.loc[0, "p_plus"])


def test_pi_minus_does_not_help_a_numu_event():
    """Same sign as the primary mu-, so no opposite-sign pair."""
    ev = frame(make_event(fs=[(-211, 0.0, 0.0, 100.0)]))
    c = D.build_candidates(ev, e_mu_min=10.0)
    p = D.event_probabilities(c, ev)
    assert p.loc[0, "p_plus"] == 0.0
    assert p.loc[0, "p_os"] == 0.0


def test_antineutrino_flips_which_sign_is_needed():
    ev = frame(make_event(nu_pdg=-14, lep_pdg=-13,
                          fs=[(-211, 0.0, 0.0, 100.0)]))
    c = D.build_candidates(ev, e_mu_min=10.0)
    p = D.event_probabilities(c, ev)
    assert p.loc[0, "p_plus"] == pytest.approx(1.0)   # primary mu+
    assert p.loc[0, "p_minus"] > 0.0
    assert p.loc[0, "p_os"] > 0.0


def test_nc_event_needs_two_decays_and_is_much_rarer():
    cc = frame(make_event(fs=[(211, 0, 0, 100.0)]))
    nc = frame(make_event(nu_pdg=14, is_cc=False, lep_pdg=14,
                          fs=[(211, 0, 0, 100.0), (-211, 0, 0, 100.0)]))
    p_cc = D.event_probabilities(D.build_candidates(cc, 10.0), cc).loc[0, "p_os"]
    p_nc = D.event_probabilities(D.build_candidates(nc, 10.0), nc).loc[0, "p_os"]
    assert p_nc < p_cc * 1e-2


def test_charm_dominates_over_a_pion_of_equal_energy():
    pion = frame(make_event(fs=[(211, 0, 0, 100.0)]))
    charm = frame(make_event(fs=[(411, 0, 0, 100.0)]))
    p_pi = D.event_probabilities(D.build_candidates(pion, 10.0), pion).loc[0, "p_os"]
    p_c = D.event_probabilities(D.build_candidates(charm, 10.0), charm).loc[0, "p_os"]
    assert p_c > 100 * p_pi


def test_include_charm_switch():
    ev = frame(make_event(fs=[(411, 0, 0, 100.0)]))
    on = D.build_candidates(ev, 10.0, include_charm=True)
    off = D.build_candidates(ev, 10.0, include_charm=False)
    assert (on["source_pdg"] == 411).any()
    assert not (off["source_pdg"] == 411).any()


def test_neutral_parent_splits_charge_evenly():
    ev = frame(make_event(fs=[(421, 0, 0, 500.0)]))
    c = D.build_candidates(ev, 10.0)
    d0 = c[c["source_pdg"] == 421]
    assert set(d0["charge"]) == {-1, 1}
    assert d0[d0.charge == 1]["p_mu"].iloc[0] == pytest.approx(
        d0[d0.charge == -1]["p_mu"].iloc[0]
    )


def test_existing_final_state_muon_is_prompt():
    """If GENIE already decayed the charm, the muon shows up directly."""
    ev = frame(make_event(fs=[(-13, 0, 0, 60.0)]))
    c = D.build_candidates(ev, 10.0)
    mu = c[c["source_pdg"] == -13].iloc[0]
    assert mu["kind"] == "prompt"
    assert mu["charge"] == 1
    assert mu["p_mu"] == 1.0
    p = D.event_probabilities(c, ev)
    assert p.loc[0, "p_os"] == pytest.approx(1.0)


def test_primary_lepton_not_double_counted():
    """Primary mu- present in the fs list as well as the lep_ columns."""
    ev = frame(make_event(lep_p=(0, 0, 500.0), fs=[(13, 0, 0, 500.0)]))
    c = D.build_candidates(ev, 10.0)
    assert (c["source_pdg"] == 13).sum() == 1


def test_independent_sources_combine_correctly():
    ev = frame(make_event(fs=[(211, 0, 0, 100.0), (211, 0, 0, 100.0)]))
    c = D.build_candidates(ev, 10.0)
    p1 = c["p_mu"].iloc[0]
    p = D.event_probabilities(c, ev)
    assert p.loc[0, "p_plus"] == pytest.approx(1 - (1 - p1) ** 2)


def test_threshold_kills_soft_hadrons():
    ev = frame(make_event(fs=[(211, 0, 0, 12.0)]))
    lo = D.event_probabilities(D.build_candidates(ev, 5.0), ev).loc[0, "p_os"]
    hi = D.event_probabilities(D.build_candidates(ev, 100.0), ev).loc[0, "p_os"]
    assert lo > 0.0
    assert hi == 0.0


# ------------------------------------------------------------- angular cut

def test_opening_angle_cut_removes_wide_pairs():
    # pi+ at ~0.1 rad from the primary muon
    ev = frame(make_event(lep_p=(0, 0, 500.0), fs=[(211, 10.0, 0.0, 100.0)]))
    c = D.build_candidates(ev, 10.0)
    wide = D.event_probabilities(c, ev, max_opening_angle=0.5).loc[0, "p_os"]
    tight = D.event_probabilities(c, ev, max_opening_angle=0.005).loc[0, "p_os"]
    assert wide > 0.0
    assert tight == 0.0


def test_angle_columns_are_sane():
    ev = frame(make_event(lep_p=(0, 0, 500.0), fs=[(211, 0.0, 0.0, 100.0)]))
    c = D.build_candidates(ev, 10.0)
    assert c["cos_to_lepton"].max() == pytest.approx(1.0)
    assert c["cos_to_beam"].max() == pytest.approx(1.0)


# --------------------------------------------------------------- summaries

def test_probability_vs_energy():
    evs = frame(
        make_event(event_id=0, e_nu=1000.0, fs=[(211, 0, 0, 100.0)]),
        make_event(event_id=1, e_nu=1000.0, fs=[]),
    )
    c = D.build_candidates(evs, 10.0)
    per = D.event_probabilities(c, evs)
    summary = D.probability_vs_energy(per)
    assert len(summary) == 1
    assert summary.loc[0, "n_events"] == 2
    assert summary.loc[0, "p_os"] == pytest.approx(per["p_os"].mean())
    assert summary.loc[0, "p_os_err"] > 0


def test_events_with_zero_candidates_still_appear():
    evs = frame(make_event(event_id=0, lep_pdg=11, nu_pdg=12, fs=[]))
    c = D.build_candidates(evs, 10.0)
    per = D.event_probabilities(c, evs)
    assert len(per) == 1
    assert per.loc[0, "p_os"] == 0.0


def test_source_breakdown_sums_to_one():
    evs = frame(make_event(fs=[(211, 0, 0, 100.0), (321, 0, 0, 100.0),
                               (411, 0, 0, 100.0)]))
    c = D.build_candidates(evs, 10.0)
    bd = D.source_breakdown(c)
    assert bd["fraction"].sum() == pytest.approx(1.0)
    # prompt primary muon dominates the *expected muon count*, so restrict
    decays = bd[bd["kind"] == "decay"]
    assert decays["expected_muons"].sum() > 0
