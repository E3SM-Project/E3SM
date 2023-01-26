from typing import Optional

from .core_parameter import CoreParameter


class ARMDiagsParameter(CoreParameter):
    def __init__(self):
        super(ARMDiagsParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.granulate.remove("seasons")
        self.test_timeseries_input = True
        self.ref_timeseries_input = True

        # Custom attributes
        # =============================
        self.test_start_yr: Optional[str] = None
        self.test_end_yr: Optional[str] = None
        self.ref_start_yr: Optional[str] = None
        self.ref_end_yr: Optional[str] = None

        # "Options include: annual_cycle", "diurnal_cycle", "diurnal_cycle_zt",
        # "pdf_daily", "convection_onset"
        self.diags_set: Optional[str] = None

        self.var_name: Optional[str] = None
        self.var_units: Optional[str] = None
