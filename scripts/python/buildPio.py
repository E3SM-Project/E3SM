#!/usr/bin/python
import sys
import platform

from machine_setup import get_machine_name
import parsers
import builder

def runBuild(compiler):
    """ run the build and configure.  Call a factory class to
        kick off the appropriate build.
    """
    platform=resolveName()
    
    print ("Try to configure and build for:: %s \n" % (platform+"_"+compiler) )
    
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
    #
    # instantiate class to get directory name from CLI
    #
    cliP = parsers.cliParser()
    args = cliP.doCliParse(useMsg)

    #
    # cmake invocation and build take place in runBuild
    #
    runBuild(args.compilerName)

if __name__ == "__main__":
    main(sys.argv[1:])
