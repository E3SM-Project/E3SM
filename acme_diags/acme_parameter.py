from __future__ import print_function

import cdp.cdp_parameter


class ACMEParameter(cdp.cdp_parameter.CDPParameter):
    def __init__(self):
        self.results_dir = '.'
        self.case_id = ''
        self.reference_data_path = ''
        self.test_data_path = ''
        self.reference_name = ''
        self.test_name = ''

        self.sets = []
        self.variables = []
        self.seasons = []
        self.regions = ['global']
        self.regrid_tool = 'esmf'
        self.regrid_method = 'linear'
        self.plevs = []

        # Plotting related
        self.main_title = 'Main Title'
        self.backend = 'vcs'
        self.save_netcdf = False
        self.output_file = 'output'
        self.output_format = ['png']
        self.canvas_size_w = 1212
        self.canvas_size_h = 1628
        self.arrows = True
        self.logo = False

        self.contour_levels = []  # used both in test and reference
        self.test_name = ''
        self.test_title = ''
        self.test_colormap = ''
        self.test_units = ''

        self.reference_name = ''
        self.reference_title = ''
        self.reference_colormap = ''
        self.reference_units = ''

        self.diff_name = ''
        self.diff_title = 'Model - Observation'
        self.diff_colormap = ''
        self.diff_levels = []
        self.diff_units = ''

        self.num_workers = 1
        self.multiprocessing = False
        self.distributed = False
        self.scheduler_addr = '127.0.0.1:8786'

    def check_values(self):
        if self.reference_data_path == '':
            print('reference_data_path is needed! Define it in the parameters file or in the command line using --reference_data_path')
            quit()
        if self.test_data_path == '':
            print('test_data_path is needed! Define it in the parametersfile or in the command line using --test_data_path')
            quit()
        if self.multiprocessing and self.distributed:
            print("Why are you trying to run the diags multiprocessed and distributedly? You can't do this, only choose one or none.")
            quit()
