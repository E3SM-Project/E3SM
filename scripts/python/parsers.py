import argparse


class cliParser():

    def doCliParse(self, useMsg):
        """ Do CLI parse with one argument
        """
        parser = argparse.ArgumentParser(description=useMsg)

        parser.add_argument('compilerName', metavar='compiler', type=str,
                            help='enter the name of a '
                            'compiler for your pio build and ctests')

        args = parser.parse_args()

        return args
