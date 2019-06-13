from __future__ import print_function

import sys
import numpy
from .core_parameter import CoreParameter


class ZonalMean2dParameter(CoreParameter):
    def __init__(self):
        super(ZonalMean2dParameter, self).__init__()
        self.plevs = numpy.logspace(2.0, 3.0, num=17)
        # Granulating plevs causes duplicate plots in this case.
        # So we don't include that parameter below.
        self.granulate = ['variables', 'seasons', 'regions']
