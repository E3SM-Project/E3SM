from __future__ import print_function

import numpy
import genutil
import cdp.cdp_metric


class RMSE(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(RMSE, self).__init__(metric_path)

    def compute(self, model, obs, axis='xy'):
        rmse = -numpy.infty
        try:
            rmse = float(genutil.statistics.rms(
                model, obs, axis=axis, weights='generate'))
        except Exception as err:
            print(err)
        return rmse
