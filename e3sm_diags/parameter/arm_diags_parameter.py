from .core_parameter import CoreParameter


class ARMDiagsParameter(CoreParameter):
    def __init__(self):
        super(ARMDiagsParameter, self).__init__()
        # A list of the reference names to run the diags on.
        self.granulate.remove("seasons")
        self.test_timeseries_input = True
        self.ref_timeseries_input = True


#    def check_values(self):
#        if not self.ref_names:
#            msg = 'You must have a value for ref_names.'
#            raise RuntimeError(msg)
