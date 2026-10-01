"""nudimu: opposite-sign dimuon background from neutrino interactions in rock.

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
_SUBMODULES = ("transport", "events", "dimuon", "workflow", "splines")


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
    "transport", "events", "dimuon", "workflow", "splines", "__version__",
]
