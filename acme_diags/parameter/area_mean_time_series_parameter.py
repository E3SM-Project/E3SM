from .core_parameter import CoreParameter


class AreaMeanTimeSeriesParameter(CoreParameter):
    def __init__(self):
        super(AreaMeanTimeSeriesParameter, self).__init__()
        # A list of the reference names to run the diags on.
        self.ref_names = []
        self.ref_timeseries_input = True
        self.test_timeseries_input = True
        # Granulating with regions doesn't make sense,
        # because we have multiple regions for each plot.
        # So keep all of the default values except regions. 
        #self.seasons = ['ANN']
        #print(dir(self))
        self.granulate.remove('regions')
        self.granulate.remove('seasons')
        

    def check_values(self):
        if not self.ref_names:
            msg = 'You have no value for ref_names. Caculate test data only'
            print(msg)


        if not (hasattr(self, 'start_yr') and hasattr(self, 'end_yr')):
            msg = "You need to define both the 'start_yr' and 'end_yr' parameter."
            raise RuntimeError(msg)


            #raise RuntimeError(msg)
