from e3sm_diags.parser.area_mean_time_series_parser import AreaMeanTimeSeriesParser
from e3sm_diags.parser.arm_diags_parser import ARMDiagsParser
from e3sm_diags.parser.core_parser import CoreParser
from e3sm_diags.parser.diurnal_cycle_parser import DiurnalCycleParser
from e3sm_diags.parser.enso_diags_parser import EnsoDiagsParser
from e3sm_diags.parser.meridional_mean_2d_parser import MeridionalMean2dParser
from e3sm_diags.parser.mp_partition_parser import MPpartitionParser
from e3sm_diags.parser.qbo_parser import QboParser
from e3sm_diags.parser.streamflow_parser import StreamflowParser
from e3sm_diags.parser.tc_analysis_parser import TCAnalysisParser
from e3sm_diags.parser.zonal_mean_2d_parser import ZonalMean2dParser
from e3sm_diags.parser.zonal_mean_2d_stratosphere_parser import (
    ZonalMean2dStratosphereParser,
)

SET_TO_PARSER = {
    "zonal_mean_xy": CoreParser,
    "zonal_mean_2d": ZonalMean2dParser,
    "zonal_mean_2d_stratosphere": ZonalMean2dStratosphereParser,
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
    "annual_cycle_zonal_mean": CoreParser,
    "lat_lon_land": CoreParser,
    "lat_lon_river": CoreParser,
    "aerosol_aeronet": CoreParser,
    "aerosol_budget": CoreParser,
    "mp_partition": MPpartitionParser,
}
