from .area_mean_time_series_parser import AreaMeanTimeSeriesParser
from .arm_diags_parser import ARMDiagsParser
from .core_parser import CoreParser
from .diurnal_cycle_parser import DiurnalCycleParser
from .enso_diags_parser import EnsoDiagsParser
from .meridional_mean_2d_parser import MeridionalMean2dParser
from .qbo_parser import QboParser
from .streamflow_parser import StreamflowParser
from .tc_analysis_parser import TCAnalysisParser
from .zonal_mean_2d_parser import ZonalMean2dParser

SET_TO_PARSER = {
    "zonal_mean_xy": CoreParser,
    "zonal_mean_2d": ZonalMean2dParser,
    "meridional_mean_2d": MeridionalMean2dParser,
    "lat_lon": CoreParser,
    "polar": CoreParser,
    "cosp_histogram": CoreParser,
    "area_mean_time_series": AreaMeanTimeSeriesParser,
    "enso_diags": EnsoDiagsParser,
    "qbo": QboParser,
    "streamflow": StreamflowParser,
    "diurnal_cycle": DiurnalCycleParser,
    "arm_diags": ARMDiagsParser,
    "tc_analysis": TCAnalysisParser,
}
