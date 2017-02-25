import cdp.cdp_parameter


class ACMEParameter(cdp.cdp_parameter.CDPParameter):
    def __init__(self):
        self.case_id = ''
        self.reference_data_path = ''
        self.test_data_path = ''
        self.reference_data_set = ''
        self.test_data_set = ''

        self.variables = ''
        self.season = ''
        self.region = ''
        self.regrid_tool = 'esmf'
        self.regrid_method = 'linear'

        self.backend = 'vcs'
        self.output_file = 'output.png'

        # VCS Options
        self.main_title = 'Main Title'

        self.test_name = ''
        self.test_title = 'Reference'
        self.test_colormap = ''
        self.test_levels = []
        self.test_units = ''

        self.reference_name = ''
        self.reference_title = 'Observation'
        self.reference_colormap = ''
        self.reference_levels = []
        self.reference_units = ''

        self.diff_name = ''
        self.diff_title = 'Model - Observation'
        self.diff_colormap = ''
        self.diff_levels = []
        self.diff_units = ''

        self.canvas_size_w = 1212
        self.canvas_size_h = 1628
        self.arrows = True
        self.logo = True

    def check_values(self):
        # just check if reference_data_path + reference_data_set and
        # test_data_path + test_data_set is valid. More checks to come.
        # Also check if the below needed attributes exist
        '''
        if self.case_id == '':
            print 'case_id is needed! Define it in the parameter file or in the command line using --case_id'
            quit()
        '''
        if self.reference_data_path == '':
            print 'reference_data_path is needed! Define it in the parameter file or in the command line using --reference_data_path'
            quit()
        if self.test_data_path == '':
            print 'test_data_path is needed! Define it in the parameter file or in the command line using --test_data_path'
            quit()
        '''
        if self.reference_data_set == '':
            print 'reference_data_set is needed! Define it in the parameter file or in the command line using -r or --reference_data_set'
            quit()
        if self.test_data_set == '':
            print 'test_data_set is needed! Define it in the parameter file or in the command line using -t or --test_data_set'
            quit()
        if self.variables == '':
            print 'variables is needed! Define it in the parameter file or in the command line using -v or --variables'
            quit()
        if self.season == '':
            print 'season is needed! Define it in the parameter file or in the command line using -s or --season'
            quit()
        '''
