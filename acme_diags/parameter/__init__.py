from .area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
from .arm_diags_parameter import ARMDiagsParameter
from .core_parameter import CoreParameter
from .diurnal_cycle_parameter import DiurnalCycleParameter
from .enso_diags_parameter import EnsoDiagsParameter
from .meridional_mean_2d_parameter import MeridionalMean2dParameter
from .qbo_parameter import QboParameter
from .streamflow_parameter import StreamflowParameter
from .tc_analysis_parameter import TCAnalysisParameter
from .zonal_mean_2d_parameter import ZonalMean2dParameter

SET_TO_PARAMETERS = {
    "zonal_mean_xy": CoreParameter,
    "zonal_mean_2d": ZonalMean2dParameter,
    "meridional_mean_2d": MeridionalMean2dParameter,
    "lat_lon": CoreParameter,
    "polar": CoreParameter,
    "cosp_histogram": CoreParameter,
    "area_mean_time_series": AreaMeanTimeSeriesParameter,
    "enso_diags": EnsoDiagsParameter,
    "qbo": QboParameter,
    "streamflow": StreamflowParameter,
    "diurnal_cycle": DiurnalCycleParameter,
    "arm_diags": ARMDiagsParameter,
    "tc_analysis": TCAnalysisParameter,
}
