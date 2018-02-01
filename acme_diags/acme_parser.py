import cdp.cdp_parser
import acme_diags.acme_parameter


class ACMEParser(cdp.cdp_parser.CDPParser):
    def __init__(self, *args, **kwargs):
        super(ACMEParser, self).__init__(
            acme_diags.acme_parameter.ACMEParameter, *args, **kwargs)

    def load_default_args(self, files=[]):
        # this has '-p' and '--parameter' reserved
        super(ACMEParser, self).load_default_args(files)

        self.add_argument(
            '-r', '--reference_data_set',
            type=str,
            dest='reference_data_set',
            help='List of observations or models that are used as a ' +
                 'reference against the test_data_set',
            required=False)

        self.add_argument(
            '--reference_data_path',
            dest='reference_data_path',
            help='Path for the reference climitologies',
            required=False)
        
        self.add_argument(
            '--ref_name',
            dest='ref_name',
            help="The string used to locate the reference file. " + 
                 "The reference file starts with the string.",
            required=False)

        self.add_argument(
            '--ref_file',
            dest='ref_file',
            help="Path to the reference file.",
            required=False)

        self.add_argument(
            '-t', '--test_data_set',
            type=str,
            dest='test_data_set',
            help='List of observations or models to test ' +
                 'against the reference_data_set',
            required=False)

        self.add_argument(
            '--test_data_path',
            dest='test_data_path',
            help='Path for the test climitologies',
            required=False)

        self.add_argument(
            '--test_file',
            dest='test_file',
            help="Path to the test file.",
            required=False)

        self.add_argument(
            '--results_dir',
            dest='results_dir',
            help='Path of where to save the results',
            required=False)

        self.add_argument(
            '--sets',
            nargs='+',
            dest='sets',
            help='Sets to use',
            required=False)

        self.add_argument(
            '-D', '--datasets',
            nargs='+',
            dest='datasets',
            help='Datasets to use. Ex: ACME or AMWG',
            required=False)

        self.add_argument(
            '-v', '--variables',
            nargs='+',
            dest='variables',
            help='Variables to use',
            required=False)

        self.add_argument(
            '--plevs',
            type=float,
            nargs='+',
            dest='plevs',
            help='Selected pressure level',
            required=False)

        self.add_argument(
            '-s', '--seasons',
            nargs='+',
            dest='seasons',
            help='Seasons to use',
            required=False)

        self.add_argument(
            '-r', '--region',
            nargs='+',
            dest='region',
            help='region to use',
            required=False)

        self.add_argument(
            '--case_id',
            dest='case_id',
            help='Defines a subdirectory to the metrics output, so multiple' +
                 'cases can be compared',
            required=False)

        self.add_argument(
            '-o', '--output_file',
            dest='output_file',
            help='Name of the output file',
            required=False)

        self.add_argument(
            '--contour_levels',
            type=float,
            nargs='+',
            dest='contour_levels',
            help='Levels for the test and reference',
            required=False)

        self.add_argument(
            '--diff_levels',
            type=float,
            nargs='+',
            dest='diff_levels',
            help='Levels for the difference',
            required=False)

        self.add_argument(
            '--reference_name',
            dest='reference_name',
            help='Name of the reference variable',
            required=False)

        self.add_argument(
            '--test_name',
            dest='test_name',
            help='Name of the test variable',
            required=False)

        self.add_argument(
            '--short_test_name',
            dest='short_test_name',
            help='User-defined test name',
            required=False)

        self.add_argument(
            '--diff_name',
            dest='diff_name',
            help='Name of the difference variable',
            required=False)

        self.add_argument(
            '--main_title',
            dest='main_title',
            help='The big title that appears on the top of the graph',
            required=False)

        self.add_argument(
            '--reference_title',
            dest='reference_title',
            help='Title for the middle graph.',
            required=False)

        self.add_argument(
            '--test_title',
            dest='test_title',
            help='Title for the top graph',
            required=False)

        self.add_argument(
            '--diff_title',
            dest='diff_title',
            help='Title for the bottom graph',
            required=False)

        self.add_argument(
            '--reference_colormap',
            dest='reference_colormap',
            help='Colormap for the middle graph.',
            required=False)

        self.add_argument(
            '--test_colormap',
            dest='test_colormap',
            help='Colormap for the top graph',
            required=False)

        self.add_argument(
            '--diff_colormap',
            dest='diff_colormap',
            help='Colormap for the bottom graph',
            required=False)

        self.add_argument(
            '--backend',
            dest='backend',
            help='Graphical backend to use',
            required=False)

        self.add_argument(
            '--multiprocessing',
            dest='multiprocessing',
            help='Run the diags using multiprocessing',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--distributed',
            dest='distributed',
            help='Run the diags distributedly',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--save_netcdf',
            dest='save_netcdf',
            help='Save the NetCDF files.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--debug',
            dest='debug',
            help='Turns debugging on, allows code to prematurely break.',
            action='store_const',
            const=True,
            required=False)
