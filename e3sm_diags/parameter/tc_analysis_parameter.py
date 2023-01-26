from .core_parameter import CoreParameter


class TCAnalysisParameter(CoreParameter):
    def __init__(self):
        super(TCAnalysisParameter, self).__init__()
        # Override existing attributes
        # ============================
        self.granulate.remove("seasons")
        self.test_timeseries_input: bool = True
        self.ref_timeseries_input: bool = True

        # Custom attributes
        # =================
        self.ref_title: str = ""

        # TODO: There should be some validation checks because these attributes
        # are used in `tc_analysis_driver.run_diag` when `self.run_type == "model_vs_model"`
        self.ref_start_yr: str = ""
        self.ref_end_yr: str = ""
        self.test_start_yr: str = ""
        self.test_end_yr: str = ""
