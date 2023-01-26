from .core_parameter import CoreParameter


class ACzonalmeanParameter(CoreParameter):
    def __init__(self):
        super(ACzonalmeanParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.granulate.remove("seasons")
