.. _compass_clean_testcase:

clean\_testcase.py
==================

This script is used to clean one or more test cases that have already been
setup.

It will remove directories and driver scripts that were generated as part of
setting up a test case.

Command-line options::

    $ ./clean_testcase.py -h
    usage: clean_testcase.py [-h] [-o CORE] [-c CONFIG] [-r RES] [-t TEST]
                             [-n NUM] [-q] [-a] [--work_dir PATH]

    This script is used to clean one or more test cases that have already been
    setup.

    It will remove directories / driver scripts that were generated as part of
    setting up a test case.

    optional arguments:
      -h, --help            show this help message and exit
      -o CORE, --core CORE  Core that contains configurations to clean
      -c CONFIG, --configuration CONFIG
                            Configuration to clean
      -r RES, --resolution RES
                            Resolution of configuration to clean
      -t TEST, --test TEST  Test name within a resolution to clean
      -n NUM, --case_number NUM
                            Case number to clean, as listed from list_testcases.py. Can be a comma delimited list of case numbers.
      -q, --quiet           If set, script will not write a command_history file
      -a, --all             Is set, the script will clean all test cases in the work_dir.
      --work_dir PATH       If set, script will clean case directories in work_dir rather than the current directory.

