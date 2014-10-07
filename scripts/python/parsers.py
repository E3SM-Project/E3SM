import argparse


class cliParser():

    def doCliParse(self, useMsg):
        """ Do CLI parse with one argument
        """
        parser = argparse.ArgumentParser(description=useMsg)

        parser.add_argument('platformName', metavar='platform', type=str,
                            help='enter the name of a '
                            'platform for your pio build')

        args = parser.parse_args()

        return args
