from __future__ import print_function

import sys
from .core_parameter import CoreParameter


class ZonalMean2dParameter(CoreParameter):
    def __init__(self):
        super(ZonalMean2dParameter, self).__init__()
        self.plevs = [100., 115.47819847, 133.35214322, 153.99265261, 177.827941,
            205.35250265, 237.13737057, 273.84196343, 316.22776602, 365.17412725,
            421.69650343, 486.96752517, 562.34132519, 649.38163158, 749.89420933,
            865.96432336, 1000]
        # Granulating plevs causes duplicate plots in this case.
        # So we don't include that parameter below.
        self.granulate = ['variables', 'seasons', 'regions']