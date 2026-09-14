# nudimu

Opposite-sign dimuon background from neutrino interactions in the rock
around the EarthShine detector. Companion to `nubkg`, which handles the
single-muon case.

The question: how often does one neutrino interaction in the rock put
**two muons of opposite charge** into the detector at the same time, faking
the A' → μ⁺μ⁻ signature?

---

## 1. Review of the existing scripts

### `gmkspl`

| Issue | Fix |
|---|---|
| `-e 1000` caps the splines at 1 TeV, but the stated range is 1–100 TeV. Above the last knot GENIE extrapolates, silently. | `-e 100000` |
| `-n 100` knots over five decades is thin. | `-n 200` |
| No `--tune`. GENIE v3 requires the same tune at spline time and generation time; leaving it implicit means the answer depends on `GENIE_TUNE` in the environment. | pin it explicitly |
| `--event-generator-list Default` is the few-GeV physics list. Above a few hundred GeV you want the high-energy DIS package (`HEDIS` / a `GHE19_*` tune). The `G18_*` tunes are validated at accelerator energies and are being extrapolated a long way past their comfort zone at 10 TeV. | check `gevgen --help` on your build for the exact names |
| `echo` in front of every line means nothing ran. | remove, or keep and pipe to a job scheduler |

**"Individual elements then add?"** Yes. Per-element is the right
granularity. The jobs are independent so they parallelise trivially, a
failed element costs one job instead of the whole set, and adding an element
later does not invalidate the others. Merge with `gspladd -f a.xml,b.xml -o
merged.xml`. The runtime is identical; the failure modes are much better.

### `gevgen`

- Mass fractions in `-t` sum to 0.985. GENIE probably normalises, but relying
  on that is a bad idea; `workflow.target_string()` normalises explicitly.
- No `--seed`. Without it the output is not reproducible.
- The `-o file,ghep,file2,hepmc` syntax is not standard `gevgen`. Generate
  `ghep` only and convert afterwards with `gntpc`. `ghep` is the lossless
  full event record; everything else is derived from it and cheap to rebuild.
  This is the same `combined/` vs `derived/` split you already use on the
  signal side.
- Only `numu` is uncommented. You need all six: `nue`/`nuebar` and NC events
  can still make two muons from the shower alone, and at these energies
  `nutau` CC gives τ → μ ν ν 17% of the time.
- Fixed-energy generation is the right choice, but for a different reason
  than a flux-weighted run: it gives you **P(dimuon | E_ν) as a response
  function**, which you then fold against the flux in quadrature offline.
  That matches the `nubkg` architecture and lets you change the flux without
  regenerating anything.

### `gntpc -f gst`

You converted these because `gst` is a flat ntuple that `uproot` reads
without a ROOT install, and it already carries `pdgf/Ef/pxf/pyf/pzf` for
every final-state particle plus the primary lepton. That is everything this
analysis needs. Keep doing it; just treat the `gst` files as regenerable.

### The thing the workflow is missing

**GENIE stops at the nuclear boundary.** It applies intranuclear FSI and then
hands back undecayed π±, K±, K⁰_L (and, depending on the decayer config,
undecayed charm) with their four-momenta. It does not decay them in flight
and it does not transport anything through the rock.

So "count the pions and kaons above 10 GeV" is not quite the right question.
A 20 GeV π⁺ in rock has a decay length of 1100 m but an interaction length of
0.5 m, so it is absorbed ~3000 times more often than it decays. The quantity
you need per hadron is

```
P(decay wins) = λ_int / (λ_int + γcτ)
```

times the muonic branching ratio times the chance the muon lands above
threshold. That is what `nudimu.transport` computes.

Two consequences worth flagging before you generate anything:

1. **This cannot be sampled event by event.** The per-pion probabilities are
   ~10⁻⁵. Getting percent-level statistics on the final number by sampling
   decays would take ~10⁸ GENIE events. Computing the probability
   analytically and carrying it as a weight gets the same answer from 10⁴.

2. **Charm almost certainly dominates.** Run `scripts/scales.py`: per
   particle at the same energy, a D⁺ is ~5×10³ times more likely than a π⁺ to
   deliver a muon above 10 GeV. Charm multiplicity is a few percent per CC
   event against tens of pions, so charm still wins by roughly two orders of
   magnitude. This is the well-measured ν_μ → μ⁻μ⁺ opposite-sign dimuon
   process (CCFR / NuTeV / CHORUS / NOMAD), σ(μμ)/σ(CC) ≈ 0.4–0.6% above
   ~100 GeV with a few-GeV cut on the second muon. **Validate against those
   measurements before believing anything.** If GENIE's charm rate at 10 TeV
   is off by a factor of two, so is your background.

3. A related surprise: charmed hadrons are *not* fully prompt in rock. At
   10 TeV a D⁰ has γcτ = 0.66 m against a 0.49 m interaction length, so only
   ~43% of them decay. The formalism handles this uniformly, but it means the
   naive "charm always decays" shortcut breaks exactly where your flux lives.

---

## 2. The rewritten workflow

```
nudimu/
  medium.py      rock: density, interaction lengths, muon dE/dx constants
  species.py     π±, K±, K⁰_L, D, Ds, Λc — mass, cτ, BR(→μ), muon z range
  transport.py   muon range; decay-vs-absorption competition; muon yields
  events.py      gst / NuHepMC → canonical parquet event table
  dimuon.py      candidate table → per-event P(μ⁺), P(μ⁻), P(opposite sign)
  workflow.py    emits the GENIE commands and a Makefile
scripts/
  scales.py      print the governing physics scales (run this first)
  analyse.py     convert / audit / dimuon CLI
tests/           48 tests
```

### Generating

```bash
python -m nudimu.workflow --out Makefile --root /data/earthshine/genie
make -j16 splines      # slow, do once
make -j16 all
```

The Makefile makes the spline step a real build dependency instead of a
comment you edit by hand, and every command is a reviewable artifact. Adjust
the run matrix in `GenieConfig` (flavours, energies, tune, statistics).

### Analysing

```bash
python scripts/analyse.py convert  genie/derived/*.gst.root --out genie/parquet
python scripts/analyse.py audit    genie/parquet/*.parquet
python scripts/analyse.py dimuon   genie/parquet/*.parquet --standoff 50 \
       --max-angle 0.005 --out dimuon_vs_energy.parquet
```

**Run `audit` first.** It tells you from the data whether GENIE left the
charmed hadrons for you (keep `include_charm=True`) or decayed them itself
(pass `--no-charm`, or you will double count).

### The calculation

Each final-state hadron contributes an independent Bernoulli probability
`p_i` of delivering a detectable muon of definite sign, so

```
P(≥1 μ⁺) = 1 − Π₊ (1 − p_i)
P(≥1 μ⁻) = 1 − Π₋ (1 − p_i)
P(opposite-sign pair) = P(≥1 μ⁺) · P(≥1 μ⁻)
```

which factorises because the `p_i` are independent. The prompt CC lepton
enters with `p = 1` if it is a muon above threshold, so a ν_μ CC event
reduces to "probability of finding one μ⁺" while NC and ν_e CC need one of
each and are ~10³ rarer. Neutral parents (K⁰_L, D⁰) are split into two
candidate rows of weight p/2.

### The angular cut is where this gets interesting for EarthShine

`--max-angle` drops candidates further than that from the primary lepton.
This is the discriminator that actually matters: the A' signal gives two
tracks a few mrad apart (Eq. A.7 of Feng–Smolinsky–Tanedo, Δθ ≈ 2 m_A'/m_X ≈
1 mrad for m_X = 1 TeV, m_A' = 500 MeV, ≈ 20 m separation over 2.5 km). A
background pair from a shower decay carries the parent hadron's p_T, which is
O(few hundred MeV), so at 100 GeV the opening angle is tens of mrad. Quote
the background both with and without the cut; the ratio is the headline
number.

### Folding into a rate

`probability_vs_energy()` gives P_OS(E_ν, flavour, CC/NC) per interaction.
Multiply that into the `nubkg` integrand you already have — it is the same
flux × cross-section × target mass quadrature, with one extra factor. Do not
rebuild the flux machinery here.

---

## 3. Known limitations, in priority order

1. **GENIE charm at 10 TeV.** Dominant uncertainty. Validate against NuTeV /
   CHORUS dimuon rates at 100–300 GeV before extrapolating.
2. **First generation only.** A hadron that interacts makes more hadrons,
   which can also decay. Ignoring the cascade underestimates the π/K branch.
   The daughters are much softer so they mostly fail the energy threshold,
   but this is worth bounding with a Z-moment style estimate.
3. **Interaction lengths** (`Medium.lambda_*`) carry ~10–20% uncertainty and
   the π/K branch scales linearly with them. `rock.scaled(1.2)` gives you the
   systematic in one line.
4. **Charm decay spectrum** is flat in z ∈ [0, 0.6], chosen to give ⟨z⟩ ≈ 0.3.
   Since charm dominates, this is the first thing to replace with a real
   decay table if the number matters.
5. **Sharp muon threshold.** A muon is in or out at `e_mu_min`, with no
   fluctuation in the radiative losses. Fine at 10s of metres, softens the
   edge at hundreds.
6. **Hadron travels before decaying**, shifting the standoff by ~0.5 m. Ignored.
7. `transport.muon_range_m` duplicates `nubkg.muon_transport`. Pick one and
   delete the other before they drift apart.
