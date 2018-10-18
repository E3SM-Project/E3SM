from __future__ import print_function

import sys
import cdp.cdp_parameter


class ACMEParameter(cdp.cdp_parameter.CDPParameter):
    def __init__(self):
        self.case_id = ''
        # The user must define these, so don't give any defaults.
        # self.reference_data_path = ''
        # self.test_data_path = ''
        self.viewer_descr = {}

        self.sets = ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']
        self.dataset = ''
        self.run_type = 'model_vs_obs'
        self.variables = []
        self.seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']
        self.regions = ['global']
        self.regrid_tool = 'esmf'
        self.regrid_method = 'conservative'
        self.plevs = []
        self.zonal_mean_2d_plevs = []

        # Plotting related.
        self.main_title = ''
        self.backend = 'mpl'
        self.save_netcdf = False
        self.output_format = ['png']
        self.output_format_subplot = []
        self.canvas_size_w = 1212
        self.canvas_size_h = 1628
        self.figsize = [8.5, 11.0]
        self.dpi = 150
        self.arrows = True
        self.logo = False

        self.contour_levels = []
        self.test_name = ''
        self.short_test_name = ''
        self.test_title = ''
        self.test_colormap = 'cet_rainbow.rgb'
        self.test_units = ''

        self.reference_name = ''
        self.reference_title = ''
        self.reference_colormap = 'cet_rainbow.rgb'
        self.reference_units = ''

        self.diff_name = ''
        self.diff_title = 'Model - Observation'
        self.diff_colormap = 'diverging_bwr.rgb'
        self.diff_levels = []
        self.diff_units = ''

        self.multiprocessing = False
        self.distributed = False
        self.num_workers = 4

        self.no_viewer = False
        self.debug = False

        self.granulate = ['variables', 'seasons', 'regions', 'plevs']
        self.selectors = ['sets', 'seasons']

    def check_values(self):
        if not hasattr(self, 'reference_data_path'):
            msg = 'You need to specify reference_data_path in the parameters file or in the command line using --reference_data_path'
            raise RuntimeError(msg)
        if not hasattr(self, 'test_data_path'):
            msg = 'You need to specify test_data_path in the parameters file or in the command line using --test_data_path'
            raise RuntimeError(msg)

