#!/usr/bin/python
import sys

import parsers
import builder


def runBuild(platform):
    """ run the build and configure.  Call a factory class to
        kick off the appropriate build.
    """
    print ("Configure and build for:: %s\n" % platform)
    bld = builder.platformBuilder.factory(platform)
    bld.cmakeCmd()
    bld.buildCmd()
    bld.testCmd()


def main(argv):
    """ everything starts here
    """
    useMsg = ("Build PIO and unit tests outside of the CESM framework."
              "Mostly for development and testing for faster turn-around."
              " >> python scripts/python/buildPio.py darwin")
    #
    # instantiate class to get directory name from CLI
    #
    cliP = parsers.cliParser()
    args = cliP.doCliParse(useMsg)

    #
    # cmake invocation and build take place in runBuild
    #
    runBuild(args.platformName)

if __name__ == "__main__":
    main(sys.argv[1:])
