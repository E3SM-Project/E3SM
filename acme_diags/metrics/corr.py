from __future__ import print_function

import numpy
import genutil
import cdp.cdp_metric


class CORR(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(CORR, self).__init__(metric_path)

    def compute(self, model, obs, axis='xy'):
        corr = -numpy.infty
        try:
            corr = float(genutil.statistics.correlation(
                model, obs, axis=axis, weights='generate'))
        except Exception as err:
            print(err)
        return corr
