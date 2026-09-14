# Building ROOT 6.36 + GENIE 3.6.2 (Pythia8) on Debian

A from-scratch recipe: bare Debian 12 to working GENIE event generation.

**Verified working** on Debian 12 (bookworm), September 2026, with:

| Component | Version |
|---|---|
| ROOT | 6.36.14 |
| GENIE | `R-3_06_02` + PR #430 config patch |
| Pythia8 | 8.317 |
| LHAPDF | 6.5.6 |
| GSL | system (Debian) |
| Compiler | Debian g++, C++17 |

Budget 3–5 hours, most of it unattended ROOT compilation.

---

## Step 0 — Environment decision: system libraries, not conda

Build everything against **Debian system libraries**, with conda/micromamba
deactivated. Reasons:

- GENIE's build is hand-rolled Makefiles with no real dependency resolution. It
  finds libraries by path and environment variable. Mixing two library prefixes
  is how you get link-time success and runtime `undefined symbol` failures.
- conda-forge ships its own `libstdc++`. Compile with Debian's `g++` but let the
  loader find conda's `libstdc++.so.6` first, and you get ABI mismatches that
  surface as segfaults inside ROOT dictionaries.
- GSL, libxml2, and log4cpp are all packaged by Debian and are exactly what GENIE
  expects.

**Cost:** PyROOT in this build binds to `/usr/bin/python3` only, and will not be
importable from conda envs. Don't fight it — install ROOT separately in your
analysis env (`micromamba install -c conda-forge root`). Two ROOTs with different
jobs is the right arrangement, not a workaround.

**Check before every step.** This is the single most common cause of mysterious
failures, and it is easy to forget after opening a new terminal:

```bash
micromamba deactivate      # repeat until no env prefix in your prompt
echo $CONDA_PREFIX         # must be empty
which g++ python3 cmake    # must all be under /usr/bin
```

---

## Step 1 — Debian prerequisites

```bash
sudo apt update
sudo apt install -y \
  build-essential cmake git dpkg-dev binutils gfortran \
  libx11-dev libxpm-dev libxft-dev libxext-dev libssl-dev \
  libgsl-dev libxml2-dev python3-dev python3-numpy \
  libpcre3-dev libglu1-mesa-dev libglew-dev libftgl-dev \
  libfftw3-dev libsqlite3-dev libxxhash-dev \
  liblzma-dev libzstd-dev liblz4-dev libtbb-dev \
  nlohmann-json3-dev libjpeg-dev libpng-dev libgif-dev libtiff-dev \
  libgl2ps-dev libcfitsio-dev libvdt-dev \
  libblas-dev liblapack-dev libcurl4-openssl-dev \
  liblog4cpp5-dev
```

Verify the two GENIE cares most about:

```bash
gsl-config --version                      # need >= 2.6
ls /usr/include/log4cpp/Category.hh
```

Note your library directory — on amd64 it is `/usr/lib/x86_64-linux-gnu`.

> **Lesson worth carrying:** when cmake can't find a library that `apt list
> --installed` says is present, you're missing the `-dev` companion package.
> Debian splits runtime libraries from headers and `.pc` files.

---

## Step 2 — Build Pythia8

```bash
cd ~
wget https://pythia.org/download/pythia83/pythia8317.tgz
tar xzf pythia8317.tgz
cd pythia8317
./configure --prefix=$HOME/pythia8317 --cxx-common='-O2 -std=c++17 -fPIC -pthread'
make -j8
make install
```

`-std=c++17` matters: ROOT 6.36 and GENIE both build as C++17, and keeping all
three consistent removes a class of link errors that are painful to diagnose.

> `WARNING: CXX not set, using g++` during configure is **benign** — it's
> reporting the fallback, and the fallback is what you want. Silence it with
> `CXX=g++ ./configure ...` if it bothers you.

Verify, including that `CXX_COMMON` actually took (Pythia accepts malformed
values silently):

```bash
grep -E 'CXX|CXX_COMMON' Makefile.inc
ls ~/pythia8317/lib/libpythia8.so
ls ~/pythia8317/include/Pythia8/Pythia.h
ls ~/pythia8317/share/Pythia8/xmldoc/ | head
```

GENIE 3.6.2 requires Pythia8 >= 8.2.45.

---

## Step 3 — Build LHAPDF 6

**This is mandatory, not optional.** GENIE 3 has a native GRV98LO implementation,
so the standard DIS chain doesn't need LHAPDF — but the HEDIS package (high-energy
DIS) has a hard `#include "LHAPDF/LHAPDF.h"` and is built unconditionally.
`--disable-lhapdf6` does not skip that code; it just leaves the header unfound and
the build dies in `Physics/HEDIS/XSection`.

Debian does not package LHAPDF, so build from source:

```bash
cd ~
wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.6.tar.gz -O LHAPDF-6.5.6.tar.gz
tar xzf LHAPDF-6.5.6.tar.gz
cd LHAPDF-6.5.6
./configure --prefix=$HOME/lhapdf --disable-python
make -j8
make install
```

> `--disable-python` is deliberate. LHAPDF's autoconf macro searches for a
> program literally named `python` on `PATH`; Debian 12 ships only `python3`, so
> configure dies with "Cannot find python in your system path" regardless of what
> you pass to `--with-python` (which sets the *module install prefix*, not the
> interpreter). GENIE links the C++ library only. If you want LHAPDF in Python,
> use `micromamba install -c conda-forge lhapdf` in your analysis env. To build
> the bindings against system Python instead:
> `sudo apt install python-is-python3 cython3` then configure without the flag.

```bash
ls ~/lhapdf/lib/libLHAPDF.so
ls ~/lhapdf/include/LHAPDF/LHAPDF.h
```

GENIE 3.6.2 requires LHAPDF6 >= 6.3.0.

PDF *data sets* are separate from the library and are only needed at runtime for
HEDIS. Fetch them later with `lhapdf install <setname>` once you know which
configuration you want.

---

## Step 4 — Check out ROOT 6.36

GENIE 3.6.2 states ROOT >= 6.22.08, but that's a floor with no ceiling — it means
"tested down to here," not "anything above works." GENIE 3.6.2 shipped July 2025;
6.36 is roughly contemporaneous. Don't use master.

```bash
cd ~
git clone https://github.com/root-project/root.git root_src
cd root_src
git fetch --tags
git tag -l 'v6-36-*'
git checkout v6-36-14        # or the highest 6-36 patch tag listed
git describe --tags          # should print the tag cleanly
```

> If `git describe` says "No names found," you cloned without tags —
> `git fetch --tags` fixes it.

---

## Step 5 — Configure ROOT

```bash
cd ~
cmake -S ~/root_src -B ~/root_build \
  -DCMAKE_INSTALL_PREFIX=$HOME/root_install \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_STANDARD=17 \
  -Dfail-on-missing=ON \
  -Dmathmore=ON \
  -Dgdml=ON \
  -Droofit=ON \
  -Dpythia8=ON \
  -DPYTHIA8_DIR=$HOME/pythia8317 \
  -DPYTHIA8_INCLUDE_DIR=$HOME/pythia8317/include \
  -DPYTHIA8_LIBRARY=$HOME/pythia8317/lib/libpythia8.so \
  -Dpyroot=ON \
  -DPython3_EXECUTABLE=/usr/bin/python3 \
  -Dxrootd=OFF \
  -Ddavix=OFF \
  -Dbuiltin_gsl=OFF
```

| Flag | Why |
|---|---|
| `mathmore=ON` | GSL-backed `ROOT::Math`. GENIE's cross-section integration needs it. |
| `gdml=ON` | GENIE's geometry drivers read GDML detector descriptions. |
| `pythia8=ON` + paths | Builds `libEGPythia8`. |
| `fail-on-missing=ON` | **The important one.** Without it ROOT silently disables options whose dependencies it can't find, and you discover it 98% into the build. Every configure failure you hit because of this flag is a failure you'd otherwise have hit much later and much more confusingly. |
| `xrootd=OFF` | Not needed for GENIE; frequent source of unrelated failures. |
| `davix=OFF` | **Required on Debian.** ROOT wants `davix >= 0.6.4` via pkg-config, but Debian ships only the `davix` runtime package with no `libdavix-dev`, so there's no `.pc` file and no headers. Since davix is just HTTPS/WebDAV remote file access, and you already disabled XRootD, turning it off is the right call rather than building davix from source. |
| `builtin_gsl=OFF` | Use Debian's GSL so ROOT and GENIE link the same one. |

**Do not pass `-Dminuit2=ON`.** It was retired as a build option — Minuit2 is now
compiled in unconditionally. Passing it (either value) gives:
`'minuit2' is no longer part of ROOT 6.36.14 build options`. It won't appear in
`root-config --features` either; confirm with `ls ~/root_install/lib/libMinuit2.so`.

**Do not pass `GSL_CONFIG_EXECUTABLE`.** CMake finds GSL on its own. That variable
wants the path to the `gsl-config` *binary*, and pointing it at an include
directory breaks detection in confusing ways. If you ever need a hint, use
`-DGSL_ROOT_DIR=/usr`.

To flip a single option without re-running the whole line:

```bash
cmake -B ~/root_build -Ddavix=OFF
```

That reuses the cache. Fine for toggling options; do a fresh configure if you
change the compiler or C++ standard.

---

## Step 6 — Build and install ROOT

```bash
cmake --build ~/root_build -j8
cmake --install ~/root_build
```

Expect 1–3 hours.

> **If you hit `Building module 'X' implicitly`:** this is a known CMake
> dependency-ordering bug (ROOT PR #22648). The dependency gets silently dropped,
> so with `-j` on, whether you hit it depends on whether the `.pcm` happens to
> exist yet. Build that dictionary target first and resume:
> ```bash
> cmake --build ~/root_build --target G__Net -j8
> cmake --build ~/root_build -j8
> ```
> Same trick for any module named in the error. `make -j1` is immune but slow.

> **If you hit `Couldn't find module with name 'X' in modulemap!`:** different
> problem. The modulemap is generated at *configure* time, so this means your
> cache and modulemap are out of sync. Reconfigure; if that fails, `rm lib/*.pcm`
> and rebuild.

### Verify before going further

```bash
source ~/root_install/bin/thisroot.sh
root-config --version
root-config --has-mathmore     # yes
root-config --has-gdml         # yes
root-config --has-pythia8      # yes
root-config --features
ls ~/root_install/lib/libMinuit2.so
```

A known-good feature list from this recipe:

```
cxx17 asimage asimage_tiff builtin_clang builtin_cling builtin_llvm builtin_openui5
clad dataframe fitsio gdml geom http imt mathmore opengl pyroot pythia8 roofit root7
rpath runtime_cxxmodules shared spectrum sqlite ssl tmva tmva-cpu tmva-cudnn
tmva-pymva tpython use_gsl_cblas vdt webgui x11 xml
```

Note `use_gsl_cblas` confirms system GSL, and the absence of `davix`/`xrootd` is
intentional. **Fix anything wrong here now** — redoing ROOT later means redoing
GENIE too.

---

## Step 7 — GENIE environment script

Write `~/genie-env.sh`:

```bash
#!/bin/bash
# Deactivate conda/micromamba before sourcing this.

source $HOME/root_install/bin/thisroot.sh

export GENIE=$HOME/Generator
export PYTHIA8=$HOME/pythia8317
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc

export LHAPDF6_INC=$HOME/lhapdf/include
export LHAPDF6_LIB=$HOME/lhapdf/lib
export LHAPDF_DATA_PATH=$HOME/lhapdf/share/LHAPDF

export LOG4CPP_INC=/usr/include
export LOG4CPP_LIB=/usr/lib/x86_64-linux-gnu
export LIBXML2_INC=/usr/include/libxml2
export LIBXML2_LIB=/usr/lib/x86_64-linux-gnu

export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE/lib:$PYTHIA8/lib:$LHAPDF6_LIB:$LD_LIBRARY_PATH
```

```bash
chmod +x ~/genie-env.sh
```

GENIE 3.6.2 locates ROOT via `root-config` rather than `ROOTSYS`, so sourcing
`thisroot.sh` is sufficient. Don't set `ROOTSYS` manually.

---

## Step 8 — Check out GENIE and apply the Pythia8 configuration

```bash
cd ~
git clone https://github.com/GENIE-MC/Generator.git
cd Generator
git checkout R-3_06_02
```

### Why this step exists

In 3.6.2, **Pythia6 is still the default hadronizer.** The Pythia8 *code* is in
the release (PRs #408 and #425), but PR #410 deliberately kept the XML
configuration pointed at Pythia6. The XMLs that switch to Pythia8 live in
**PR #430, which is still open** as of September 2026.

Skip this and GENIE builds cleanly, then fails at *runtime* asking for a Pythia6
algorithm that isn't there.

Pythia6 is not an escape route: the ROOT/Pythia6 interface was removed in ROOT
6.32, so on any modern ROOT, Pythia8 is the only option.

### Applying it

Staying on the release tag and patching is better than moving to master. Master is
14 months past this tag and unreleased; a tag plus a saved patch file is a
reproducible recipe you can describe in a methods section.

```bash
cd ~/Generator
git remote add nusense https://github.com/nusense/Generator.git
git fetch nusense rhatcher/pythia8_xml_config_update

BASE=$(git merge-base HEAD nusense/rhatcher/pythia8_xml_config_update)
git diff --stat $BASE nusense/rhatcher/pythia8_xml_config_update
```

Expected — 11 files, 24 insertions, 16 deletions:

```
 config/AGKY2019.xml                   |  2 +-
 config/AGKYLowW2019.xml               |  4 ++--
 config/DISHadronicSystemGenerator.xml |  2 +-
 config/Messenger.xml                  |  1 +
 config/Messenger_inuke_verbose.xml    |  3 +++
 config/Messenger_laconic.xml          |  2 ++
 config/Messenger_rambling.xml         |  1 +
 config/Messenger_whisper.xml          |  1 +
 config/UnstableParticleDecayer.xml    |  6 +++---
 configure                             | 10 +++++-----
 src/make/Make.include                 |  8 ++++----
```

If you see C++ source files in that list, stop — the branch has moved and this
recipe needs revisiting.

**Take the whole patch, not just `config/`.** The `configure` and `Make.include`
hunks carry the enable-Pythia8-by-default logic and the library link ordering.
XMLs asking for Pythia8 with a build that didn't link it is a bad place to be.

```bash
git diff $BASE nusense/rhatcher/pythia8_xml_config_update > /tmp/pythia8.patch
git apply --check /tmp/pythia8.patch
git apply /tmp/pythia8.patch
git diff --stat                        # must match the 11-file list above
```

> **If `--check` fails with "patch does not apply" or "does not match index":**
> the most likely cause is that part of it is *already applied* from an earlier
> attempt. `git apply` is atomic — if any file fails, nothing is written, even for
> files it reports as clean. Check with
> `git diff nusense/rhatcher/pythia8_xml_config_update -- config/` (empty output
> means the XMLs already match the branch), then apply only what's missing:
> `git diff $BASE nusense/... -- configure src/make/Make.include | git apply`.
> Only if that fails is it a genuine conflict; then `git apply -3` for a real
> three-way merge, or hand-edit — it's only 7 algorithm-name substitutions and 8
> inserted log-stream lines.

> The substitutions are Pythia**6**→Pythia**8** algorithm names
> (`Pythia6Hadro2019` → `Pythia8Hadro2019`), *not* a 2019→2023 year bump. The 2019
> suffix stays. (`AGCharmPythia8Hadro2023` exists separately, which makes a blanket
> year rename tempting and wrong.)

Save your patch — this file plus `R-3_06_02` is the reproducible recipe:

```bash
git diff > ~/genie-pythia8-config.patch
```

---

## Step 9 — Configure GENIE

```bash
micromamba deactivate
echo $CONDA_PREFIX                 # must be empty
source ~/genie-env.sh
which root-config                  # must be ~/root_install/bin/root-config

cd ~/Generator
./configure \
  --prefix=$HOME/genie_install \
  --enable-pythia8 \
  --disable-pythia6 \
  --disable-cernlib \
  --disable-lhapdf5 \
  --enable-lhapdf6 \
  --with-lhapdf6-inc=$LHAPDF6_INC \
  --with-lhapdf6-lib=$LHAPDF6_LIB \
  --enable-flux-drivers \
  --enable-geom-drivers \
  --enable-atmo \
  --disable-masterclass \
  --disable-debug \
  --with-optimiz-level=O2 \
  --with-pythia8-inc=$PYTHIA8/include \
  --with-pythia8-lib=$PYTHIA8/lib \
  --with-libxml2-inc=$LIBXML2_INC \
  --with-libxml2-lib=$LIBXML2_LIB \
  --with-log4cpp-inc=$LOG4CPP_INC \
  --with-log4cpp-lib=$LOG4CPP_LIB
```

**Do not add `--enable-nucleon-decay`.** That package fails to build in this
configuration — `libGPhNDcy` never gets produced, and every app then fails to link
with `cannot find -lGPhNDcy`. It's also not what you want if your interest is
baryon-number violation in heavy-flavour decays; GENIE's module simulates
proton/bound-neutron decay for detectors like Super-K and DUNE.

Optional packages default to off, so simply omitting a flag disables it. Verify:

```bash
grep GOPT_ENABLE src/make/Make.config
```

Confirm `PYTHIA8=YES`, `PYTHIA6=NO`, `LHAPDF6=YES`, `NUCLEON_DECAY=NO`. Also skim
the `GOPT_WITH_*` paths and check ROOT resolves to `$HOME/root_install`.

---

## Step 10 — Build GENIE

```bash
make -j8
make install
```

> **`-j8` hides errors.** A package can fail while seven other directories keep
> printing success, and the failure scrolls away. If `make install` complains or
> `bin/` is empty, don't guess — rebuild the suspect directory serially:
> ```bash
> cd ~/Generator/src/Apps        # or whichever package
> make 2>&1 | tee /tmp/genie-build.log
> ```
> The last thing on screen is the real error.

> **`make install` says "Previous installation exists... Try 'make distclean'
> first":** do **not** run `make distclean` — it wipes your source tree's build
> products and costs you the entire rebuild. The check is a crude
> does-this-path-exist guard. Just `rm -rf ~/genie_install` and re-run
> `make install`.

> **`Building Professor2 was not enabled. Skipping...`** at the end of install is
> informational, not an error. Professor2 is an external model-tuning package you
> didn't enable.

> **Changing configure options:** package enable/disable doesn't need
> `make clean` — just reconfigure and `make -j8`, and it'll pick up the
> difference in a couple of minutes. But changing `--with-optimiz-level`, the C++
> standard, or a library path *does* require `make clean`: GENIE's `.d` files
> track headers, not compiler flags, so you'd otherwise get mixed object files.

### Verify

```bash
ls ~/Generator/bin/ | wc -l           # expect 22
ls ~/genie_install/bin/ | wc -l       # must match
```

You should see `gevgen`, `gmkspl`, `gntpc`, `gevdump`, `genie-config` among them.

---

## Step 11 — Confirm you got a Pythia8 build

```bash
source ~/genie-env.sh
genie-config                    # bare invocation prints usage; --help prints nothing
genie-config --features
genie-config --has-pythia8      # yes
genie-config --has-pythia6      # no
```

This is the direct answer to the question all of Step 8 was about.

Save `genie-config --libs` output somewhere — you'll need it to compile anything
against GENIE, and it's tedious to reconstruct.

---

## Step 12 — Smoke test

```bash
cd ~/scratch                    # gevgen writes its ntuple to $PWD
gevgen -n 10 -p 14 -t 1000060120 -e 1 \
       --event-generator-list Default \
       --seed 12345
```

1 GeV muon neutrinos on carbon-12.

> **Do not pass `--cross-sections none`.** That is not a keyword meaning "compute
> on the fly" — GENIE treats it as a literal filename and dies with
> `Input cross-section file [none] does not exist!`. Omit the flag entirely to
> compute splines on the fly.

Without splines this takes several minutes of apparently-idle CPU before output
appears. That's normal. Then:

```bash
ls -la gntp.0.ghep.root
gevdump -f gntp.0.ghep.root -n 2
```

Hadronic final states in the event records mean Pythia8 is working and the Step 8
patch did its job.

For real work, precompute splines once per probe/target:

```bash
gmkspl -p 14 -t 1000060120 -n 100 -e 10 -o carbon_splines.xml
gevgen -n 1000 -p 14 -t 1000060120 -e 1 \
       --cross-sections carbon_splines.xml --seed 12345
```

---

## Two installs, one authority

You now have GENIE in two places: `~/Generator` (build tree) and
`~/genie_install`. The env script above points `$GENIE` at the build tree. Pick
one and stay with it, or you'll eventually run binaries from one against config
XMLs from the other. Confirm the install got your patched config:

```bash
diff -r ~/Generator/config ~/genie_install/config | head
```

---

## Failure modes, ranked by likelihood

1. **conda/micromamba active.** Check `echo $CONDA_PREFIX` before every build step.
2. **Step 8 skipped or half-applied.** Builds clean, fails at runtime. Verify with `genie-config --has-pythia8` and a real `gevgen` run.
3. **Missing `-dev` package.** `apt list --installed` showing a library means nothing; cmake needs headers and `.pc` files.
4. **`-j8` hiding the real error.** Rebuild the failing directory serially.
5. **`PYTHIA8DATA` unset.** Pythia8 init fails with an unhelpful message about missing XML data.
6. **A second ROOT on `PATH`.** `which root-config` before configuring GENIE.
7. **C++ standard mismatch** across Pythia8 / ROOT / GENIE. Keep all three at C++17.
