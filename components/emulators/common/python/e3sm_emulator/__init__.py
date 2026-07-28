"""Python emulators for E3SM's emulator components.

The C++ side of the bridge lives in
``components/emulators/common/src/inference``; it imports
:mod:`e3sm_emulator.bridge` and calls ``create_emulator(config)``.

Nothing here needs anything but numpy, so the whole bridge can be
exercised -- and unit tested -- on a machine with no ML stack installed.
"""

from .context import Context

__all__ = ["Context"]
