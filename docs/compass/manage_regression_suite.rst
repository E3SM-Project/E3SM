.. _compass_manage_regression_suite:

manage\_regression\_suite.py
============================

This script is used to manage regression suites. A regression suite is a set of
test cases that ensure one or more features in a model meet certain criteria.

Using this script one can setup or clean a regression suite.

When setting up a regression suite, this script will generate a script to run
all tests in the suite, and additionally setup each individual test case.

When cleaning a regression suite, this script will remove any generated files
for each individual test case, and the run script that runs all test cases.

Command-line options::

    $ ./manage_regression_suite.py -h
    usage: manage_regression_suite.py [-h] -t FILE [-f FILE] [-s] [-c] [-v]
                                      [-m FILE] [-b PATH] [--work_dir PATH]

    This script is used to manage regression suites. A regression suite is a set of
    test cases that ensure one or more features in a model meet certain criteria.

    Using this script one can setup or clean a regression suite.

    When setting up a regression suite, this script will generate a script to run
    all tests in the suite, and additionally setup each individual test case.

    When cleaning a regression suite, this script will remove any generated files
    for each individual test case, and the run script that runs all test cases.

    optional arguments:
      -h, --help            show this help message and exit
      -t FILE, --test_suite FILE
                            Path to file containing a test suite to setup
      -f FILE, --config_file FILE
                            Configuration file for test case setup
      -s, --setup           Option to determine if regression suite should be setup or not.
      -c, --clean           Option to determine if regression suite should be cleaned or not.
      -v, --verbose         Use verbose output from setup_testcase.py
      -m FILE, --model_runtime FILE
                            Definition of how to build model run commands on this machine
      -b PATH, --baseline_dir PATH
                            Location of baseslines that can be compared to
      --work_dir PATH       If set, script will setup the test suite in work_dir rather in this script's location.

