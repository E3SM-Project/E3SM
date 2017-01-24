import cdp.cdp_parameter


class ACMEParameter(cdp.cdp_parameter.CDPParameter):
    def __init__(self):
        self.reference_data_path = ''
        self.test_data_path = ''
        self.reference_data_set = ''
        self.test_data_set = ''
        self.variables = ''
        self.season = ''
        self.reference_colormap = ''
        self.test_colormap = ''
        self.diff_colormap = ''


    def check_values(self):
        # just check if reference_data_path + reference_data_set and
        # test_data_path + test_data_set is valid. More checks to come.
        # Also check if the below needed attributes exist
        if self.case_id == '':
            print 'case_id is needed! Define it in the parameter file or in the command line using --case_id'
            quit()
        if self.reference_data_path == '':
            print 'reference_data_path is needed! Define it in the parameter file or in the command line using --reference_data_path'
            quit()
        if self.test_data_path == '':
            print 'test_data_path is needed! Define it in the parameter file or in the command line using --test_data_path'
            quit()
        if self.reference_data_set == '':
            print 'reference_data_set is needed! Define it in the parameter file or in the command line using -r or --reference_data_set'
            quit()
        if self.test_data_set == '':
            print 'test_data_set is needed! Define it in the parameter file or in the command line using -t or --test_data_set'
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
