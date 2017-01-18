import ast
import cdp.cdp_parser
import acme_diags.acme_parameter


class ACMEParser(cdp.cdp_parser.CDPParser):
    def __init__(self, *args, **kwargs):
        super(ACMEParser, self).__init__(acme_diags.acme_parameter.ACMEParameter, *args, **kwargs)

    def load_default_args(self):
        # this has '-p' and '--parameter' reserved
        super(ACMEParser, self).load_default_args()

        self.add_argument(
            '--case_id',
            dest='case_id',
            help='Defines a subdirectory to the metrics output, so multiple' +
                 'cases can be compared',
            required=False)

        self.add_argument(
            '-r', '--reference_data_set',
            type=str,
            nargs='+',
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
            '-t', '--test_data_set',
            type=str,
            nargs='+',
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
            '-v', '--variables',
            type=str,
            dest='variables',
            help='Variables to use',
            required=False)

        self.add_argument(
            '-s', '--season',
            type=str,
            nargs='+',
            dest='season',
            help='Season to use',
            required=False)
