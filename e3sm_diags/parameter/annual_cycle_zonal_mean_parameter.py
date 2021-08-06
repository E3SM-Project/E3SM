from .core_parameter import CoreParameter


class ACzonalmeanParameter(CoreParameter):
    def __init__(self):
        super(ACzonalmeanParameter, self).__init__()
        # A list of the reference names to run the diags on.
        self.granulate.remove("seasons")
