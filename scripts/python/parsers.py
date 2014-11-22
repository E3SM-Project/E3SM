import argparse


class cliParser():

    def doCliParse(self, useMsg):
        """ Do CLI parse with one argument
        """
        parser = argparse.ArgumentParser(description=useMsg)

        parser.add_argument('--test',action='store_true',help='run the ctests')

        parser.add_argument('--backtrace', action='store_true',
                            help='show exception backtraces as extra debugging '
                            'output')

        parser.add_argument('--debug', action='store_true',
                            help='extra debugging output')

        parser.add_argument('--mach', nargs=1, 
                            help='target machine name')

        parser.add_argument('--compiler', nargs=1, required=False,
                            help='target compiler name')

        parser.add_argument('--xmlpath', nargs=1, required=True,
                            help='path to CESM machines directory')

        parser.add_argument('--mpilib', nargs=1, required=False,
                            help='target mpi library')

        args = parser.parse_args()

        return args
