from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter
from e3sm_diags.parser.core_parser import CoreParser


class ARMDiagsParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=ARMDiagsParameter, *args, **kwargs)  # type: ignore

    def add_arguments(self):
        super().add_arguments()

        self.parser.add_argument(
            "--ref_names",
            type=str,
            nargs="+",
            dest="ref_names",
            help="List of reference names.",
            required=False,
        )
