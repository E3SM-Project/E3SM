from acme_diags.parameter.arm_diags_parameter import ARMDiagsParameter

from .core_parser import CoreParser


class ARMDiagsParser(CoreParser):
    def __init__(self, *args, **kwargs):
        if "parameter_cls" in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(parameter_cls=ARMDiagsParameter, *args, **kwargs)

    def load_default_args(self, files=[]):
        # This has '-p' and '--parameter' reserved.
        super().load_default_args(files)

        self.add_argument(
            "--ref_names",
            type=str,
            nargs="+",
            dest="ref_names",
            help="List of reference names.",
            required=False,
        )
