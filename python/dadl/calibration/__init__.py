"""Calibration engine (phase 2).

Each calibration estimates detector parameters from data the experiment
already produces and returns a `Calibration`: point estimates with 95%
intervals, a posterior sample table, fit checks and the hashes of its inputs.
`dadl.artifact` turns it into an immutable, content-addressed artifact.

Calibrations estimate and record only. Applying them to analyses is phase 3
(one forward model, no chained corrections).
"""

from dadl.calibration.accidentals import calibrate_accidentals
from dadl.calibration.base import Calibration
from dadl.calibration.dark import calibrate_dark
from dadl.calibration.dcr import calibrate_dcr_temperature
from dadl.calibration.delay import calibrate_delay
from dadl.calibration.heralding import calibrate_heralding

__all__ = [
    "Calibration",
    "calibrate_accidentals",
    "calibrate_dark",
    "calibrate_dcr_temperature",
    "calibrate_delay",
    "calibrate_heralding",
]
