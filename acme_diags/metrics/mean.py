import cdutil
import cdp.cdp_metric

class Mean(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(Mean, self).__init__(metric_path)

    def compute(self, variable):
        return cdutil.averager(variable, axis='xy', weights='generate')
