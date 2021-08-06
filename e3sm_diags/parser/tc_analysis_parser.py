from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter

from .core_parser import CoreParser


class TCAnalysisParser(CoreParser):
    def __init__(self, *args, **kwargs):
        if "parameter_cls" in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(parameter_cls=TCAnalysisParameter, *args, **kwargs)

    def load_default_args(self, files=[]):
        # This has '-p' and '--parameter' reserved.
        super().load_default_args(files)
