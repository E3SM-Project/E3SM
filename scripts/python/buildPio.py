import os
import platform
import sys
# imports of NCAR scripts
import builder
import parsers
lib_path = os.path.join('scripts/python/contrib/unit_testing')
sys.path.append(lib_path)
from machine_setup import get_machine_name


def runBuild(compiler):
    """ run the build and configure.  Call a factory class to
        kick off the appropriate build.
    """
    platform = resolveName()

    print ("Configure, build and test for :: %s \n" %
           (platform + "_" + compiler))

    bld = builder.platformBuilder.factory(platform+"_"+compiler)
    bld.metaBuild()


def resolveName():
    """ osx requires a Darwin build (for laptops, desktops).
        bigger platforms will usually have a name (yellowstone,
        edison, etc...)
    """
    if platform.system() == "Darwin":
        name = "darwin"
    else:
        name = get_machine_name()
    return name


def main(argv):
    """ everything starts here
    """
    useMsg = ("Build PIO and unit tests outside of the CESM framework."
              "Mostly for development and testing for faster turn-around."
              " >> python scripts/python/buildPio.py <compiler>")

    cliP = parsers.cliParser()
    args = cliP.doCliParse(useMsg)

    runBuild(args.compilerName)

if __name__ == "__main__":
    main(sys.argv[1:])
