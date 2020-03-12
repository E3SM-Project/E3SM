.. _compass_list_testcases:

list\_testcases.py
==================

This script is used to list available test cases.

It iterates through the directory structure and prints out configuration
options to setup specific test cases. Additionally, the ``-o``, ``-c``, ``-r``, and ``-t``
flags can be used to narrow the information that this script prints. If any of
them are passed in, the script will only print test cases that match all
criteria.

Additionally, if ``-n`` is passed in to get information about a single test case,
it will only print the flags needed to setup that specific test case.

Command-line options::

    $ ./list_testcases.py -h
    usage: list_testcases.py [-h] [-o CORE] [-c CONFIG] [-r RES] [-t TEST]
                             [-n NUMBER]

    This script is used to list available test cases.

    It iterates through the directory structure and prints out configuration
    options to setup specific test cases. Additionally, the -o, -c, -r, and -t
    flags can be used to narrow the information that this script prints. If any of
    them are passed in, the script will only print test cases that match all
    criteria.

    Additionally, if -n is passed in to get information about a single test case,
    it will only print the flags needed to setup that specific test case.

    optional arguments:
      -h, --help            show this help message and exit
      -o CORE, --core CORE  Core to search for configurations within
      -c CONFIG, --configuration CONFIG
                            Configuration name to search for
      -r RES, --resolution RES
                            Resolution to search for
      -t TEST, --test TEST  Test name to search for
      -n NUMBER, --number NUMBER
                            If set, script will print the flags to use a the N'th configuration.

