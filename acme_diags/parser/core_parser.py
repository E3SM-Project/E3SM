import cdp.cdp_parser
from acme_diags.parameter.core_parameter import CoreParameter


class CoreParser(cdp.cdp_parser.CDPParser):
    def __init__(self, *args, **kwargs):
        if 'parameter_cls' in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(parameter_cls=CoreParameter, *args, **kwargs)

    def load_default_args(self, files=[]):
        # This has '-p' and '--parameter' reserved.
        super().load_default_args(files)

        self.add_argument(
            'set_name',
            type=str,
            help='Name of the diags set to send ' + 
                 'these arguments to.',
            nargs='?')

        self.add_argument(
            '-r', '--reference_data_set',
            type=str,
            dest='reference_data_set',
            help='List of observations or models that are used as a ' +
                 'reference against the test_data_set.',
            required=False)

        self.add_argument(
            '--reference_data_path',
            dest='reference_data_path',
            help='Path for the reference data.',
            required=False)

        self.add_argument(
            '--ref_timeseries_input',
            dest='ref_timeseries_input',
            help='The input reference data are timeseries files.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--ref_start_yr',
            dest='ref_start_yr',
            help="Start year for the reference timeseries files.",
            required=False)

        self.add_argument(
            '--ref_end_yr',
            dest='ref_end_yr',
            help="End year for the reference timeseries files.",
            required=False)

        self.add_argument(
            '--ref_start_time_slice',
            dest='ref_start_time_slice',
            help="Starting time slice year for the reference timeseries files.",
            required=False)

        self.add_argument(
            '--ref_end_time_slice',
            dest='ref_end_time_slice',
            help="Ending time slice year for the reference timeseries files.",
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
                 'against the reference_data_set.',
            required=False)

        self.add_argument(
            '--test_data_path',
            dest='test_data_path',
            help='Path for the test data.',
            required=False)

        self.add_argument(
            '--test_timeseries_input',
            dest='test_timeseries_input',
            help='The input test data are timeseries files.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--test_start_yr',
            dest='test_start_yr',
            help="Start year for the test timeseries files.",
            required=False)

        self.add_argument(
            '--test_end_yr',
            dest='test_end_yr',
            help="End year for the test timeseries files.",
            required=False)

        self.add_argument(
            '--test_start_time_slice',
            dest='test_start_time_slice',
            help="Starting time slice year for the test timeseries files.",
            required=False)

        self.add_argument(
            '--test_end_time_slice',
            dest='test_end_time_slice',
            help="Ending time slice year for the test timeseries files.",
            required=False)

        self.add_argument(
            '--test_file',
            dest='test_file',
            help="Path to the test file.",
            required=False)

        self.add_argument(
            '--results_dir',
            dest='results_dir',
            help='Path of where to save the results.',
            required=False)

        self.add_argument(
            '--sets',
            nargs='+',
            dest='sets',
            help='Sets to use.',
            required=False)

        self.add_argument(
            '-D', '--dataset',
            dest='dataset',
            help="Dataset to use. Ex: 'ACME' or 'AMWG'.",
            required=False)

        self.add_argument(
            '--run_type',
            dest='run_type',
            help="What comparison to do. One of three options: "
            + "'model_vs_obs'/'obs_vs_model', 'model_vs_model', or 'obs_vs_obs'.",
            required=False)

        self.add_argument(
            '-v', '--variables',
            nargs='+',
            dest='variables',
            help='Variables to use.',
            required=False)

        self.add_argument(
            '--plevs',
            type=float,
            nargs='+',
            dest='plevs',
            help='Selected pressure level.',
            required=False)

        self.add_argument(
            '--plot_plevs',
            dest='plot_plevs',
            help='plot specified plevs',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--plot_log_plevs',
            dest='plot_log_plevs',
            help='plot plevs on log-scale',
            action='store_const',
            const=True,
            required=False)


        self.add_argument(
            '-s', '--seasons',
            nargs='+',
            dest='seasons',
            help='Seasons to use.',
            required=False)

        self.add_argument(
            '-r', '--regions',
            nargs='+',
            dest='regions',
            help='regions to use.',
            required=False)

        self.add_argument(
            '--regrid_tool',
            dest='regrid_tool',
            help="What regrid tool to use.",
            required=False)

        self.add_argument(
            '--regrid_method',
            dest='regrid_method',
            help="What regrid method for the regrid tool to use.",
            required=False)

        self.add_argument(
            '--case_id',
            dest='case_id',
            help='Defines a subdirectory to the metrics output, so multiple' +
                 'cases can be compared.',
            required=False)

        self.add_argument(
            '--output_format',
            nargs='+',
            dest='output_format',
            help="What output format the plots should be saved in. "
                 + "Possible values are: ['png', 'pdf', 'svg'].",
            required=False)

        self.add_argument(
            '--output_format_subplot',
            nargs='+',
            dest='output_format_subplot',
            help="What output format the individual subplots should be saved in (leave empty for no subplots)."
                 + "Possible values are: ['png', 'pdf', 'svg'].",
            required=False)

        self.add_argument(
            '--canvas_size_w',
            type=int,
            dest='canvas_size_w',
            help="Size in pixels of the width for the output figure. "
                 + "VCS only.",
            required=False)

        self.add_argument(
            '--canvas_size_h',
            type=int,
            dest='canvas_size_h',
            help="Size in pixels of the height for the output figure. "
                 + "VCS only.",
            required=False)

        self.add_argument(
            '--figsize',
            type=float,
            nargs='+',
            dest='figsize',
            help="Width and height like so: [width, height]. "
                 + "Matplotlib only.",
            required=False)

        self.add_argument(
            '--dpi',
            type=int,
            dest='dpi',
            help="DPI to use. "
                 + "Matplotlib only.",
            required=False)

        self.add_argument(
            '--arrows',
            dest='arrows',
            help='Display arrows on the plot.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--logo',
            dest='logo',
            help='Display the logo. VCS only.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--contour_levels',
            type=float,
            nargs='+',
            dest='contour_levels',
            help='Levels for the test and reference plots.',
            required=False)

        self.add_argument(
            '--diff_levels',
            type=float,
            nargs='+',
            dest='diff_levels',
            help='Levels for the difference plot.',
            required=False)

        self.add_argument(
            '--reference_name',
            dest='reference_name',
            help='Name of the reference variable.',
            required=False)

        self.add_argument(
            '--test_name',
            dest='test_name',
            help='Name of the test variable.',
            required=False)

        self.add_argument(
            '--short_test_name',
            dest='short_test_name',
            help='User-defined test name.',
            required=False)

        self.add_argument(
            '--diff_name',
            dest='diff_name',
            help='Name of the difference variable.',
            required=False)

        self.add_argument(
            '--main_title',
            dest='main_title',
            help='The big title that appears on the top of the graph.',
            required=False)

        self.add_argument(
            '--reference_title',
            dest='reference_title',
            help='Title for the middle graph.',
            required=False)

        self.add_argument(
            '--test_title',
            dest='test_title',
            help='Title for the top graph.',
            required=False)

        self.add_argument(
            '--diff_title',
            dest='diff_title',
            help='Title for the bottom graph.',
            required=False)

        self.add_argument(
            '--reference_colormap',
            dest='reference_colormap',
            help='Colormap for the middle graph.',
            required=False)

        self.add_argument(
            '--test_colormap',
            dest='test_colormap',
            help='Colormap for the top graph.',
            required=False)

        self.add_argument(
            '--diff_colormap',
            dest='diff_colormap',
            help='Colormap for the bottom graph.',
            required=False)

        self.add_argument(
            '--reference_units',
            dest='reference_units',
            help='Units to use for the middle graph.',
            required=False)

        self.add_argument(
            '--test_units',
            dest='test_units',
            help='Units to use for the top graph.',
            required=False)

        self.add_argument(
            '--diff_units',
            dest='diff_units',
            help='Units to use for the bottom graph.',
            required=False)

        self.add_argument(
            '--backend',
            dest='backend',
            help='Graphical backend to use.',
            required=False)

        self.add_argument(
            '--multiprocessing',
            dest='multiprocessing',
            help='Run the diags using multiprocessing.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--distributed',
            dest='distributed',
            help='Run the diags distributedly.',
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
            '--no_viewer',
            dest='no_viewer',
            help="Don't generate the viewer.",
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
