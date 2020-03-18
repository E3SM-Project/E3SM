from .core_parameter import CoreParameter


class EnsoDiagsParameter(CoreParameter):
    def __init__(self):
        super(EnsoDiagsParameter, self).__init__()
        self.granulate.remove('seasons')
        self.nino_region = 'NINO34'
        self.plot_type = 'map'
        self.print_statements = False
        self.ref_timeseries_input = True
        self.test_timeseries_input = True

    def check_values(self):
        valid_nino_regions = ['NINO3', 'NINO34', 'NINO4']
        if self.nino_region not in valid_nino_regions:
            msg = 'nino_region={} not in {}'.format(
                self.nino_region, valid_nino_regions) 
            raise RuntimeError(msg)

        valid_plot_types = ['map', 'scatter']
        if self.plot_type not in valid_plot_types:
            msg = 'plot_type={} not in {}'.format(
                self.plot_type, valid_plot_types)
            raise RuntimeError(msg)

        if not (hasattr(self, 'start_yr') and hasattr(self, 'end_yr')):
            msg = "You need to define both the 'start_yr' and 'end_yr' parameter."
            raise RuntimeError(msg)
