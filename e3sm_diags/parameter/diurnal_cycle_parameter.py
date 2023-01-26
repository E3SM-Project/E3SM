from .core_parameter import CoreParameter


class DiurnalCycleParameter(CoreParameter):
    def __init__(self):
        super(DiurnalCycleParameter, self).__init__()
        # Custom attributes
        # =============================
        self.normalize_amp_int = 0
        self.print_statements = False
        self.normalize_test_amp = False
