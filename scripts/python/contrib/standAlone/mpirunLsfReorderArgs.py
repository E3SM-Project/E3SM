import sys
import platform
import argparse
import subprocess


def main(argv):
    """ Stand alone script.  Should only depend on python packages provided
        with original install.  Switches the order of -np <N> and the CLI
        executable so that we can automatically call mpirun.lsf from
        PIOs cmake/ctest framework.
    """
    useMsg = ("Provides a wrapper to mpirun.lsf that allows us to reroder "
              "arguments so we can run PIO's unit tests (via ctest) quickly "
              "and from CLI as we are developing PIO.")

    parser = argparse.ArgumentParser(description=useMsg)
    parser.add_argument('-np', metavar='np', type=str,
                        help='enter the name of a '
                        'numprocs to send to mpirun.lsf in the form '
                        '< -np numPes > ')
    parser.add_argument('exe', metavar='exe', type=str,
                        help='enter the name of a '
                        'command and any arguments to send to mpirun.lsf')
    args = parser.parse_args()

    exeCmd = ("mpirun.lsf " + args.exe + " -np " + args.np)

    print exeCmd

    p = subprocess.Popen(exeCmd, shell=True)
    p.wait()

if __name__ == "__main__":
    main(sys.argv[1:])
