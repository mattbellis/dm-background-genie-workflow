#!/usr/bin/env python3
"""Rebuild the nudimu package layout from a flattened unzip.

The zip download flattens directories, so all the .py files land next to
each other and `import nudimu` fails.  Run this once in the directory
containing the unzipped files:

    python setup_layout.py

It moves each file to where it belongs, writes the nudimu/__init__.py that
was missing from the download, and verifies that the package imports.
Safe to re-run.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

# where each file belongs, relative to the project root
LAYOUT = {
    "medium.py": "nudimu",
    "species.py": "nudimu",
    "transport.py": "nudimu",
    "events.py": "nudimu",
    "dimuon.py": "nudimu",
    "workflow.py": "nudimu",
    "scales.py": "scripts",
    "analyse.py": "scripts",
    "test_transport.py": "tests",
    "test_dimuon.py": "tests",
    "test_io_and_workflow.py": "tests",
    "README.md": ".",
    "pyproject.toml": ".",
}

INIT_PY = '''"""nudimu: opposite-sign dimuon background from neutrino interactions in rock.

Companion to the ``nubkg`` single-muon background package.  Where ``nubkg``
asks "how often does a neutrino in the rock put one muon in the detector",
this asks "how often does it put two, of opposite sign, at the same time".

Pipeline
--------
1. :mod:`nudimu.workflow`   emit reproducible GENIE spline/generation commands
2. :mod:`nudimu.events`     GENIE output -> canonical parquet event table
3. :mod:`nudimu.dimuon`     event table -> per-event dimuon probability
4. fold P(E_nu) against your flux with the existing nubkg quadrature

The step that GENIE does *not* do for you is step 3: GENIE hands back
undecayed pi/K (and, depending on tune, undecayed charm) at the nuclear
boundary and stops.  Whether those become muons is a competition between
decay in flight and absorption in the rock, and that competition is what
:mod:`nudimu.transport` computes.
"""

from .medium import STANDARD_ROCK, Medium  # noqa: F401
from .species import SPECIES, Species, lookup  # noqa: F401

# Submodules are imported lazily so that `python -m nudimu.workflow` does not
# trip the "module found in sys.modules before execution" runpy warning, and
# so that importing nudimu does not pull in uproot/pyhepmc.
_SUBMODULES = ("transport", "events", "dimuon", "workflow")


def __getattr__(name):
    if name in _SUBMODULES:
        import importlib
        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(_SUBMODULES))


__version__ = "0.1.0"

__all__ = [
    "STANDARD_ROCK", "Medium", "SPECIES", "Species", "lookup",
    "transport", "events", "dimuon", "workflow", "__version__",
]
'''

PYPROJECT = '''[build-system]
requires = ["setuptools>=64"]
build-backend = "setuptools.build_meta"

[project]
name = "nudimu"
version = "0.1.0"
description = "Opposite-sign dimuon background from neutrino interactions in rock"
requires-python = ">=3.9"
dependencies = ["numpy", "pandas", "pyarrow"]

[project.optional-dependencies]
gst = ["uproot"]
hepmc = ["pyhepmc"]
dev = ["pytest"]

[tool.setuptools.packages.find]
include = ["nudimu*"]

[tool.pytest.ini_options]
pythonpath = [".", "tests"]
testpaths = ["tests"]
'''


def main() -> int:
    root = Path.cwd()
    print(f"repairing layout under {root}")

    for sub in ("nudimu", "tests", "scripts"):
        (root / sub).mkdir(exist_ok=True)

    moved, already, missing = 0, 0, []
    for name, dest in LAYOUT.items():
        src = root / name
        target = root / dest / name if dest != "." else root / name
        if target.exists() and not src.exists():
            already += 1
            continue
        if not src.exists():
            missing.append(name)
            continue
        if src.resolve() == target.resolve():
            already += 1
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(target))
        print(f"  {name} -> {dest}/")
        moved += 1

    init = root / "nudimu" / "__init__.py"
    if not init.exists():
        init.write_text(INIT_PY)
        print("  wrote nudimu/__init__.py  (this was missing from the download)")

    pp = root / "pyproject.toml"
    if not pp.exists():
        pp.write_text(PYPROJECT)
        print("  wrote pyproject.toml")

    print(f"\nmoved {moved}, already in place {already}")
    if missing:
        print(f"not found (fine if you did not download them): {', '.join(missing)}")

    # verify
    sys.path.insert(0, str(root))
    try:
        import nudimu  # noqa: F401
        from nudimu import transport
        r = transport.muon_range_m(100.0)
    except Exception as exc:
        print(f"\nimport still failing: {exc}")
        return 1
    print(f"\nimport ok: nudimu {nudimu.__version__}, "
          f"100 GeV muon range = {r:.1f} m in standard rock")

    print("\nfinal layout:")
    for p in sorted(root.rglob("*")):
        if "__pycache__" in p.parts or p.is_dir():
            continue
        if p.suffix in (".py", ".toml", ".md"):
            print(f"  {p.relative_to(root)}")

    print("\nnext:")
    print("  python -m pytest -q")
    print("  python scripts/scales.py")
    print("  python -m nudimu.workflow --out Makefile --root /data/earthshine/genie")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
