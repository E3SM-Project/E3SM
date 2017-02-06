import cdp.cdp_metric

class Min(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(Min, self).__init__(metric_path)

    def compute(self, variable):
        return float(variable.min())
