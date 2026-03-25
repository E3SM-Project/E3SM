# shared_data.py

# this module defines the shared data structures and functions for the landgen workflow
# this includes data shared between modules, and data that are needed for multiprocessing
# data shared via multiprocessing need to be instantiated as mp.Manager() objects in main()
# custom data strucures associated managers are also defined and registered here

# the shared custom data structures are registered in main()

import multiprocessing as mp
from multiprocessing.managers import BaseManager
import numpy as np

# set some default dimension sizes here for reference; these can be changed at allocation time
# see classes below for descriptions of these dimensions
## todo: add these to the config file and pass them to allocate()

# LtData dimensions
n_pfts_default = 51
n_harvest_default = 5
n_grazing_default = 2
n_elev_default = 61
n_elev_edges_default = 62
n_density_default = 3
n_month_default = 12
n_levurb_default = 10
n_rad_default = 2
n_solar_default = 2
n_vocveg_default = 6

# GridData dimensions
n_cells_default = 1
n_vertices_default = 4

# TopoData dimensions
n_levslp_default = 11

########## shared data structures not needed for multiprocessing
# are there any?

########## multiprocessing shared data structures


# ---------------------------------------------------------------------------
# GridData: custom data structure for grid geometry data shared across
# multiprocessing workers.  All arrays are float64 except cell_id (int64).
# This can be used for any grid, but we currently use a healpix grid
# n_cells: number of 'land' cells in the landgen out grid
#    lon_xy and lat_xy are the cell 'center' coordinates
#    lon_vtx and lat_vtx are the vertex coordinates for each cell
#    n_vertices: number of vertices per cell (4 for quadrilateral cells)
# landfrac is determined by the topography processing, so initialize to 1 here
# ---------------------------------------------------------------------------
class GridData:
    """Per-cell grid geometry container (cell ids, coordinates, landfrac)."""

    def __init__(self):

        self.num_cells   = None    # int     - number of 'land' cells in out grid
        self.num_vertices   = None    # int  - number of vertices per cell
        self.cell_id  = None    # int64   [n_cells]
        self.lon_xy   = None    # float64 [n_cells]
        self.lat_xy   = None    # float64 [n_cells]
        self.cell_area = None   # float64 [n_cells]
        self.landfrac = None    # float64 [n_cells]
        self.lon_vtx  = None    # float64 [n_cells, n_vertices]
        self.lat_vtx  = None    # float64 [n_cells, n_vertices]
        

    def allocate(self, n_cells=n_cells_default, n_vertices=n_vertices_default):
        """Allocate all arrays given dimension sizes."""
        self.num_cells   = n_cells
        self.num_vertices = n_vertices
        self.cell_id  = np.zeros(n_cells,               dtype=np.int64)
        self.lon_xy   = np.zeros(n_cells,               dtype=np.float64)
        self.lat_xy   = np.zeros(n_cells,               dtype=np.float64)
        self.cell_area = np.zeros(n_cells,              dtype=np.float64)
        self.landfrac = np.ones(n_cells,                dtype=np.float64)
        self.lon_vtx  = np.zeros((n_cells, n_vertices), dtype=np.float64)
        self.lat_vtx  = np.zeros((n_cells, n_vertices), dtype=np.float64)


# ---------------------------------------------------------------------------
# GridManager: custom BaseManager that can vend GridData proxy objects to
# worker processes.  Register, instantiate, and start in landgen.main().
# ---------------------------------------------------------------------------
class GridManager(BaseManager):
    pass



# ---------------------------------------------------------------------------
# TopoData: custom data structure for topography data shared across
# multiprocessing workers.  All arrays are float64 except where noted.
# Dimensions are set at allocation time; fields are initialised to None.
# n_cells: number of 'land' cells in the landgen grid
# n_levslp: number of slope percentage levels (11)
# ---------------------------------------------------------------------------
class TopoData:
    """Per-cell topography data container for the landgen workflow."""

    def __init__(self):
        self.topo            = None    # float64 [n_cells]          - topographic height
        self.std_elev        = None    # float32 [n_cells]          - standard deviation of elevation
        self.slope           = None    # float64 [n_cells]          - mean slope
        self.slp_p10         = None    # float64 [n_cells, n_levslp] - slope percentiles
        self.sinsl_cosas     = None    # float64 [n_cells]          - sin(slope)*cos(aspect)
        self.sinsl_sinas     = None    # float64 [n_cells]          - sin(slope)*sin(aspect)
        self.sky_view        = None    # float64 [n_cells]          - sky view factor
        self.terrain_config  = None    # float64 [n_cells]          - terrain configuration factor
        self.fmax            = None    # float64 [n_cells]          - maximum fractional saturated area

    def allocate(self, n_cells=n_cells_default, n_levslp=n_levslp_default):
        """Allocate all arrays given dimension sizes."""
        g = n_cells
        self.topo            = np.zeros(g,              dtype=np.float64)
        self.std_elev        = np.zeros(g,              dtype=np.float32)
        self.slope           = np.zeros(g,              dtype=np.float64)
        self.slp_p10         = np.zeros((g, n_levslp),  dtype=np.float64)
        self.sinsl_cosas     = np.zeros(g,              dtype=np.float64)
        self.sinsl_sinas     = np.zeros(g,              dtype=np.float64)
        self.sky_view        = np.zeros(g,              dtype=np.float64)
        self.terrain_config  = np.zeros(g,              dtype=np.float64)
        self.fmax            = np.zeros(g,              dtype=np.float64)


# ---------------------------------------------------------------------------
# TopoManager: custom BaseManager that can vend TopoData proxy objects to
# worker processes.  Register, instantiate, and start in landgen.main().
# ---------------------------------------------------------------------------
class TopoManager(BaseManager):
    pass




## todo: add some aggregate land cover arrays to ltData for the simple land cover data

# ---------------------------------------------------------------------------
# LtData: custom data structure for land-type data shared across
# multiprocessing workers.  All arrays are float64 except cell_id (int64).
# Dimensions are set at allocation time; fields are initialised to None.
# These data are on the landgen grid defined by GridData 
# n_cells: number of 'land' cells in the landgen grid
# n_pfts: number of plant functional types (51, includes bare and crop functional types)
# n_harvest: number of harvest types (5: luh categories)
# n_grazing: number of grazing types (2: pasture (grass, intensive) and rangeland)
# n_elev: number of elevation bins for glacier cover (currently 61, may change)
# n_elev_edges: number of elevation bin edges for glacier cover (n_elev + 1)
# n_density: number of urban density bins (3: low, medium, high)
# n_month: number of months in the year (12)
# n_levurb: number of vertical levels in the urban canopy model (10: from the ground up)
# n_rad: number of radiation bands (2: shortwave and longwave)
# n_solar: number of solar radiation categories (2: direct and diffuse)
# n_vocveg: number of vegetation types for VOC emissions (6: 3 tree, shrub, grass, crop)
# ---------------------------------------------------------------------------
class LtData:
    """Per-cell landcover data container for the landgen workflow."""

    def __init__(self):

        # 1-D grid arrays  [n_cells]
        self.pct_ocean          = None         # float64
        self.lake_depth         = None         # float64
        self.lake_depth_mask    = None         # float64
        self.pct_lake           = None         # float64
        self.pct_wetland        = None         # float64
        self.pct_glacier        = None         # float64
        self.canyon_hwr         = None         # float64
        self.em_improad         = None         # float64
        self.em_perroad         = None         # float64
        self.em_roof            = None         # float64
        self.em_wall            = None         # float64
        self.ht_roof            = None         # float64
        self.thick_roof         = None         # float64
        self.thick_wall         = None         # float64
        self.t_building_min     = None         # float64
        self.t_building_max     = None         # float64
        self.wind_hgt_canyon    = None         # float64
        self.wtlunit_roof       = None         # float64
        self.wtroad_perv        = None         # float64
        self.nlev_improad       = None         # float64
        self.firrig             = None         # float64
        self.fsurf              = None         # float64
        self.fgrd               = None         # float64
        self.Nmanure_mixed      = None         # float64
        self.Nmanure_pastures   = None         # float64
        self.fract_nitr         = None         # float64
        self.fract_urea         = None         # float64
        self.soilph             = None         # float64
        self.abm                = None         # float64

        # 2-D arrays  [n_cells, second_index]
        self.pct_pft            = None         # [n_cells, n_pfts]
        self.nfert              = None         # [n_cells, n_pfts]
        self.pfert              = None         # [n_cells, n_pfts]
        self.harvest_frac       = None         # [n_cells, n_harvest]
        self.harvest_mass       = None         # [n_cells, n_harvest]
        self.grazing_frac       = None         # [n_cells, n_grazing]
        self.pct_glc_gic        = None         # [n_cells, n_elev]
        self.pct_glc_icesheet   = None         # [n_cells, n_elev]
        self.bin_centers        = None         # [n_elev]
        self.bin_edges          = None         # [n_elev_edges]
        self.pct_urban          = None         # [n_cells, n_density]
        self.monthly_lai        = None         # [n_cells, n_pfts, n_month]
        self.monthly_sai        = None         # [n_cells, n_pfts, n_month]
        self.monthly_height_top = None         # [n_cells, n_pfts, n_month]
        self.monthly_height_bot = None         # [n_cells, n_pfts, n_month]
        self.veg_voc_emit       = None         # [n_cells, n_vocveg]

        # 3-D arrays  [n_cells, n_rad, n_solar]
        self.alb_improad        = None         # [n_cells, n_rad, n_solar]
        self.alb_perroad        = None         # [n_cells, n_rad, n_solar]
        self.alb_roof           = None         # [n_cells, n_rad, n_solar]
        self.alb_wall           = None         # [n_cells, n_rad, n_solar]

        # 2-D urban layer arrays  [n_cells, n_levurb]
        self.tk_roof            = None         # [n_cells, n_levurb]
        self.tk_wall            = None         # [n_cells, n_levurb]
        self.tk_improad         = None         # [n_cells, n_levurb]
        self.cv_roof            = None         # [n_cells, n_levurb]
        self.cv_wall            = None         # [n_cells, n_levurb]
        self.cv_improad         = None         # [n_cells, n_levurb]

    def allocate(self, n_cells=n_cells_default, n_pfts=n_pfts_default, n_harvest=n_harvest_default,
                n_grazing=n_grazing_default, n_elev=n_elev_default, n_elev_edges=n_elev_edges_default,
                n_density=n_density_default, n_month=n_month_default, n_levurb=n_levurb_default,
                n_rad=n_rad_default, n_solar=n_solar_default, n_vocveg=n_vocveg_default):
        """Allocate all arrays given dimension sizes."""
        self.pct_ocean          = np.zeros(n_cells, dtype=np.float64)
        self.lake_depth         = np.zeros(n_cells, dtype=np.float64)
        self.lake_depth_mask    = np.zeros(n_cells, dtype=np.float64)
        self.pct_lake           = np.zeros(n_cells, dtype=np.float64)
        self.pct_wetland        = np.zeros(n_cells, dtype=np.float64)
        self.pct_glacier        = np.zeros(n_cells, dtype=np.float64)
        self.canyon_hwr         = np.zeros(n_cells, dtype=np.float64)
        self.em_improad         = np.zeros(n_cells, dtype=np.float64)
        self.em_perroad         = np.zeros(n_cells, dtype=np.float64)
        self.em_roof            = np.zeros(n_cells, dtype=np.float64)
        self.em_wall            = np.zeros(n_cells, dtype=np.float64)
        self.ht_roof            = np.zeros(n_cells, dtype=np.float64)
        self.thick_roof         = np.zeros(n_cells, dtype=np.float64)
        self.thick_wall         = np.zeros(n_cells, dtype=np.float64)
        self.t_building_min     = np.zeros(n_cells, dtype=np.float64)
        self.t_building_max     = np.zeros(n_cells, dtype=np.float64)
        self.wind_hgt_canyon    = np.zeros(n_cells, dtype=np.float64)
        self.wtlunit_roof       = np.zeros(n_cells, dtype=np.float64)
        self.wtroad_perv        = np.zeros(n_cells, dtype=np.float64)
        self.nlev_improad       = np.zeros(n_cells, dtype=np.float64)
        self.firrig             = np.zeros(n_cells, dtype=np.float64)
        self.fsurf              = np.zeros(n_cells, dtype=np.float64)
        self.fgrd               = np.zeros(n_cells, dtype=np.float64)
        self.Nmanure_mixed      = np.zeros(n_cells, dtype=np.float64)
        self.Nmanure_pastures   = np.zeros(n_cells, dtype=np.float64)
        self.fract_nitr         = np.zeros(n_cells, dtype=np.float64)
        self.fract_urea         = np.zeros(n_cells, dtype=np.float64)
        self.soilph             = np.zeros(n_cells, dtype=np.float64)
        self.abm                = np.zeros(n_cells, dtype=np.float64)
        self.pct_pft            = np.zeros((n_cells, n_pft),      dtype=np.float64)
        self.nfert              = np.zeros((n_cells, n_pft),      dtype=np.float64)
        self.pfert              = np.zeros((n_cells, n_pft),      dtype=np.float64)
        self.harvest_frac       = np.zeros((n_cells, n_harvest),  dtype=np.float64)
        self.harvest_mass       = np.zeros((n_cells, n_harvest),  dtype=np.float64)
        self.grazing_frac       = np.zeros((n_cells, n_grazing),  dtype=np.float64)
        self.pct_glc_gic        = np.zeros((n_cells, n_elev),     dtype=np.float64)
        self.pct_glc_icesheet   = np.zeros((n_cells, n_elev),     dtype=np.float64)
        self.bin_centers        = np.zeros(n_elev,          dtype=np.float64)
        self.bin_edges          = np.zeros(n_elev_edges,    dtype=np.float64)
        self.pct_urban          = np.zeros((n_cells, n_density),  dtype=np.float64)
        self.monthly_lai        = np.zeros((n_cells, n_pft,  n_month), dtype=np.float64)
        self.monthly_sai        = np.zeros((n_cells, n_pft,  n_month), dtype=np.float64)
        self.monthly_height_top = np.zeros((n_cells, n_pft,  n_month), dtype=np.float64)
        self.monthly_height_bot = np.zeros((n_cells, n_pft,  n_month), dtype=np.float64)
        self.veg_voc_emit       = np.zeros((n_cells, numvocveg),  dtype=np.float64)
        self.alb_improad        = np.zeros((n_cells, numrad, numsolar), dtype=np.float64)
        self.alb_perroad        = np.zeros((n_cells, numrad, numsolar), dtype=np.float64)
        self.alb_roof           = np.zeros((n_cells, numrad, numsolar), dtype=np.float64)
        self.alb_wall           = np.zeros((n_cells, numrad, numsolar), dtype=np.float64)
        self.tk_roof            = np.zeros((n_cells, n_levurb),   dtype=np.float64)
        self.tk_wall            = np.zeros((n_cells, n_levurb),   dtype=np.float64)
        self.tk_improad         = np.zeros((n_cells, n_levurb),   dtype=np.float64)
        self.cv_roof            = np.zeros((n_cells, n_levurb),   dtype=np.float64)
        self.cv_wall            = np.zeros((n_cells, n_levurb),   dtype=np.float64)
        self.cv_improad         = np.zeros((n_cells, n_levurb),   dtype=np.float64)


# ---------------------------------------------------------------------------
# LtManager: custom BaseManager that can vend LtData proxy objects to
# worker processes.  Register, instantiate, and start in landgen.main().
# ---------------------------------------------------------------------------
class LtManager(BaseManager):
    pass







