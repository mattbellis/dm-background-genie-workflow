"""Opposite-sign dimuon probability, event by event.

The calculation is deliberately *not* a Monte Carlo over decays.  Each
final-state hadron contributes an independent Bernoulli probability p_i of
delivering a detectable muon of a definite sign.  Then

    P(at least one mu+) = 1 - prod_{i in +} (1 - p_i)
    P(at least one mu-) = 1 - prod_{i in -} (1 - p_i)
    P(opposite-sign pair) = P(mu+) * P(mu-)

which factorises because the p_i are independent.  The prompt lepton from
a CC vertex enters with p = 1 if it is a muon above threshold and 0
otherwise, so nu_mu CC events reduce to "probability of finding one mu+",
while NC and nu_e CC events need one of each and are correspondingly rarer.

Neutral parents (K0L, D0) yield a muon of either sign with equal
probability, so they are split into two candidate rows of weight p/2.

Threshold model
---------------
A muon is "detectable" if its energy at creation exceeds
``e_mu_min``.  Pass a scalar if you want a fixed cut, or use
``transport.muon_threshold_energy(distance_m, medium)`` to convert a
required standoff distance into an energy.  The muon is assumed to be
created at the interaction vertex; a hadron that travels one interaction
length first shifts the effective standoff by ~0.5 m, which is negligible
compared with the tens to hundreds of metres that set the threshold.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from .medium import STANDARD_ROCK, Medium
from .species import lookup
from .transport import mean_muon_energy, muon_yield_probability

__all__ = [
    "build_candidates",
    "event_probabilities",
    "probability_vs_energy",
    "MUON_PDG",
]

MUON_PDG = 13
_P_TOL = 1e-6


def _unit(px, py, pz):
    p = np.sqrt(px * px + py * py + pz * pz)
    if p <= 0:
        return 0.0, 0.0, 1.0
    return px / p, py / p, pz / p


def build_candidates(events: pd.DataFrame, e_mu_min: float,
                     medium: Medium = STANDARD_ROCK,
                     include_charm: bool = True) -> pd.DataFrame:
    """Long-format table: one row per (event, potential muon source).

    Parameters
    ----------
    events
        Canonical event table from :mod:`nudimu.events`.
    e_mu_min
        Muon energy at creation required for the muon to count.
    include_charm
        Set False if :func:`nudimu.events.audit_final_states` shows that
        GENIE already decayed the charmed hadrons for you.

    Returns
    -------
    DataFrame with columns
        event_id, e_nu, is_cc, nu_pdg, kind, source_pdg, charge,
        e_parent, p_mu, e_mu_mean, cos_to_lepton, cos_to_beam
    """
    rows = []
    for evt in events.itertuples(index=False):
        lep_dir = _unit(evt.lep_px, evt.lep_py, evt.lep_pz)
        lep_is_muon = abs(int(evt.lep_pdg)) == MUON_PDG

        fs_pdg = np.asarray(evt.fs_pdg)
        fs_e = np.asarray(evt.fs_e, dtype=float)
        fs_px = np.asarray(evt.fs_px, dtype=float)
        fs_py = np.asarray(evt.fs_py, dtype=float)
        fs_pz = np.asarray(evt.fs_pz, dtype=float)

        primary_in_fs = False

        for j in range(len(fs_pdg)):
            pdg = int(fs_pdg[j])
            e = float(fs_e[j])
            d = _unit(fs_px[j], fs_py[j], fs_pz[j])
            cos_lep = float(np.dot(d, lep_dir))
            cos_beam = d[2]

            if abs(pdg) == MUON_PDG:
                # already a muon: prompt, no decay weight
                if lep_is_muon and abs(e - float(evt.lep_e)) < _P_TOL:
                    primary_in_fs = True
                rows.append(
                    dict(
                        event_id=int(evt.event_id), e_nu=float(evt.e_nu),
                        is_cc=bool(evt.is_cc), nu_pdg=int(evt.nu_pdg),
                        kind="prompt", source_pdg=pdg,
                        charge=-1 if pdg > 0 else +1,   # pdg 13 is mu-
                        e_parent=e, p_mu=float(e >= e_mu_min),
                        e_mu_mean=e, cos_to_lepton=cos_lep, cos_to_beam=cos_beam,
                    )
                )
                continue

            sp = lookup(pdg)
            if sp is None:
                continue
            if sp.lambda_key == "charm" and not include_charm:
                continue

            p = float(muon_yield_probability(sp, e, e_mu_min, medium))
            if p <= 0.0:
                continue
            e_mu = float(mean_muon_energy(sp, e, e_mu_min))

            if sp.charge == 0:
                charges = ((+1, 0.5), (-1, 0.5))
            else:
                # pi+ -> mu+, K+ -> mu+, D+ -> mu+ ; sign follows the parent
                charges = ((sp.charge, 1.0),)

            for q, frac in charges:
                rows.append(
                    dict(
                        event_id=int(evt.event_id), e_nu=float(evt.e_nu),
                        is_cc=bool(evt.is_cc), nu_pdg=int(evt.nu_pdg),
                        kind="decay", source_pdg=pdg, charge=q,
                        e_parent=e, p_mu=p * frac, e_mu_mean=e_mu,
                        cos_to_lepton=cos_lep, cos_to_beam=cos_beam,
                    )
                )

        # primary lepton not present in the final-state list (gst sometimes
        # stores it separately): add it explicitly.
        if lep_is_muon and not primary_in_fs:
            pdg = int(evt.lep_pdg)
            rows.append(
                dict(
                    event_id=int(evt.event_id), e_nu=float(evt.e_nu),
                    is_cc=bool(evt.is_cc), nu_pdg=int(evt.nu_pdg),
                    kind="prompt", source_pdg=pdg,
                    charge=-1 if pdg > 0 else +1,
                    e_parent=float(evt.lep_e),
                    p_mu=float(float(evt.lep_e) >= e_mu_min),
                    e_mu_mean=float(evt.lep_e),
                    cos_to_lepton=1.0, cos_to_beam=lep_dir[2],
                )
            )

    cols = ["event_id", "e_nu", "is_cc", "nu_pdg", "kind", "source_pdg",
            "charge", "e_parent", "p_mu", "e_mu_mean",
            "cos_to_lepton", "cos_to_beam"]
    return pd.DataFrame(rows, columns=cols)


def event_probabilities(candidates: pd.DataFrame, events: pd.DataFrame,
                        max_opening_angle: Optional[float] = None) -> pd.DataFrame:
    """Per-event P(mu+), P(mu-) and P(opposite-sign pair).

    Parameters
    ----------
    max_opening_angle
        If given (radians), drop candidates whose direction is further than
        this from the primary lepton.  This is how you ask the question
        that actually matters for EarthShine: the A' signal gives two
        tracks a few mrad apart, so a background pair with a 100 mrad
        opening angle is not a background at all.
    """
    cand = candidates
    if max_opening_angle is not None:
        cand = cand[cand["cos_to_lepton"] >= np.cos(max_opening_angle)]

    if len(cand) == 0:
        base = events[["event_id", "e_nu", "is_cc", "nu_pdg"]].copy()
        base["p_plus"] = 0.0
        base["p_minus"] = 0.0
        base["p_os"] = 0.0
        return base

    cand = cand.assign(_q=np.log1p(-np.clip(cand["p_mu"], 0.0, 1.0 - 1e-15)))
    grouped = cand.groupby(["event_id", "charge"])["_q"].sum().unstack(fill_value=0.0)
    for q in (-1, 1):
        if q not in grouped.columns:
            grouped[q] = 0.0

    p_plus = 1.0 - np.exp(grouped[1].to_numpy())
    p_minus = 1.0 - np.exp(grouped[-1].to_numpy())

    out = pd.DataFrame(
        {"event_id": grouped.index.to_numpy(),
         "p_plus": p_plus, "p_minus": p_minus, "p_os": p_plus * p_minus}
    )
    base = events[["event_id", "e_nu", "is_cc", "nu_pdg"]].merge(
        out, on="event_id", how="left"
    )
    return base.fillna({"p_plus": 0.0, "p_minus": 0.0, "p_os": 0.0})


def probability_vs_energy(per_event: pd.DataFrame,
                          by: tuple = ("e_nu", "nu_pdg")) -> pd.DataFrame:
    """Mean opposite-sign dimuon probability per interaction, grouped.

    The statistical error quoted is the standard error on the mean of the
    per-event weights, which is the right uncertainty for a weighted
    estimator and will be far smaller than sqrt(N_pairs)/N would suggest
    for a sampled one.
    """
    g = per_event.groupby(list(by))["p_os"]
    out = g.agg(["mean", "std", "count"]).reset_index()
    out = out.rename(columns={"mean": "p_os", "count": "n_events"})
    out["p_os_err"] = out["std"] / np.sqrt(out["n_events"])
    return out.drop(columns=["std"])


def source_breakdown(candidates: pd.DataFrame) -> pd.DataFrame:
    """Which species is carrying the answer, summed over events.

    Sum of p_mu per species, which to leading order in small p is the
    expected number of detectable muons contributed by that species.
    """
    out = (
        candidates.groupby(["kind", "source_pdg"])["p_mu"]
        .agg(["sum", "count"])
        .reset_index()
        .rename(columns={"sum": "expected_muons", "count": "n_candidates"})
        .sort_values("expected_muons", ascending=False)
    )
    total = out["expected_muons"].sum()
    out["fraction"] = out["expected_muons"] / total if total > 0 else np.nan
    return out
