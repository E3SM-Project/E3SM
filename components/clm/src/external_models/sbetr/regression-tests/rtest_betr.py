#!/usr/bin/env python
"""BeTR regression test driver

Author: Ben Andre <andre@ucar.edu>

"""

from __future__ import print_function

import sys

if sys.hexversion < 0x02070000:  # pragma: no coverage
    print(70 * "*")
    print("ERROR: {0} requires python >= 2.7.x. ".format(sys.argv[0]))
    print("It appears that you are running python {0}".format(
        ".".join(str(x) for x in sys.version_info[0:3])))
    print(70 * "*")
    sys.exit(1)

#
# built-in modules
#
import argparse
import copy
from datetime import datetime
import logging
import os
import os.path
import subprocess
import time
import traceback

if sys.version_info[0] == 2:  # pragma: no coverage
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser

#
# installed dependencies
#

#
# other modules in this package
#

#
# globals
#
BANNER = 70*'-'
SEPERATOR = '-+- {0}'.format(35*'-')


# -----------------------------------------------------------------------------
#
# User input
#
# -----------------------------------------------------------------------------
def commandline_options(ext_args=None):
    """Process the command line arguments.

    NOTE(bja, 201603): ext_args should only be used by the unit tests!
    The caller should pass a list of args. We default to none so that
    the arg parser will default back to sys.argv if nothing is
    specifed.

    """
    parser = argparse.ArgumentParser(
        description='FIXME: python program template.')

    parser.add_argument('--backtrace', action='store_true',
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument('--check-only', action='store_true',
                        help='check regression output from a previous run')

    parser.add_argument('--config', nargs='+', default=[],
                        help='path to config file(s)')

    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')

    parser.add_argument('--dry-run', action='store_true',
                        help='setup but don\'t run tests')

    parser.add_argument('--executable', nargs=1, required=True,
                        help='path to the executable')

    parser.add_argument('--timeout', nargs=1, default='90.0',
                        help='max runtime [seconds] before we timout a test.')

    parser.add_argument('--update-baseline', action='store_true',
                        help=('update the baseline regression file with '
                              'the results of the current run.'))

    options = parser.parse_args(ext_args)
    return options


def read_config_file(filename):
    """Read the configuration file and process

    """
    logging.info("Reading configuration file :")
    logging.info("    {0}".format(filename))

    cfg_file = os.path.abspath(filename)
    if not os.path.isfile(cfg_file):
        raise RuntimeError("Could not find file: {0}".format(cfg_file))

    config = config_parser()
    config.read(cfg_file)

    return config


def read_config_file_as_dict(filename):
    """Read the configuration file and return a dictionary

    """
    config = read_config_file(filename)
    config_dict = config_to_dict(config)
    return config_dict


# -----------------------------------------------------------------------------
#
# Classes
#
# -----------------------------------------------------------------------------
class RegressionTestSuite(object):
    """
    """

    def __init__(self, filename, conf, timeout):
        """
        """
        name = os.path.basename(filename)
        logging.info('Adding test suite from "{0}"'.format(name))
        self._name = name.split('.')[0]

        self._test_dir = os.path.abspath(os.path.dirname(filename))
        # default comparison object for this test suite
        # tests will get a copy and modify as needed.
        self._comparison = Comparison(self._name, conf)

        self._tests = []
        self._add_tests(conf, timeout)

    def _add_tests(self, conf, timeout):
        """
        """
        os.chdir(self._test_dir)

        for test_name in conf:
            comparison = copy.deepcopy(self._comparison)
            options = conf[test_name]
            try:
                test = RegressionTest(test_name, comparison, options, timeout)
                self._tests.append(test)
            except RuntimeError as e:
                msg = "{0} : {1}".format(self._name, e)
                raise RuntimeError(msg)

    def run_tests(self, executable, check_only, dry_run, update_baseline):
        """Run the tests in the test suite

        Note: we want to update the baselines from here rather than in
        a separate function call to the test suite so the update
        baseline occurs next to the test results in the global test
        log.

        """
        logging.info(BANNER)
        logging.info('Running tests for "{0}"'.format(self._name))
        os.chdir(self._test_dir)
        for test in self._tests:
            test.run(executable, check_only, dry_run)
            if update_baseline:
                test.update_baseline(dry_run)

    def num_tests(self):
        """Return the number of tests in the test suite
        """
        return len(self._tests)

    def status_summary(self):
        """Return the test results for the test suite.

        """
        status = StatusSummary()

        for test in self._tests:
            if test.status() == 'skip':
                status.add_skip()
            elif test.status() == 'fail':
                status.add_failure()
            elif test.status() == 'pass':
                status.add_pass()
            else:
                msg = ('development error: suite "{0}" : test "{1}" : '
                       'indeterminant status.'.format(self._name, test._name))
                raise RuntimeError(msg)
        if len(self._tests) != status.total():
            msg = ('development error: suite "{0}" : number of tests did not'
                   ' equal the number of valid statuses.'.format(self._name))
            raise RuntimeError(msg)
        return status


class RegressionTest(object):
    """
    """
    _REGRESSION = 'regression'
    _BASELINE = 'regression.baseline'

    def __init__(self, name, comparison, options, timeout):
        """Note: assume that we are in the test directory.

        We start with a None status because we want to explicitly se
        the pass/fail/skip status. If it is still unset at the end of
        the test, then there is a bug in the comparison code!

        """
        logging.info('  Adding test "{0}"'.format(name))

        # initialize to defaults
        self._name = name
        self._namelist_filename = None
        self._comparison = comparison
        self._status = None
        self._timeout = 60.0

        # setup the test
        self._comparison.set_name(self._name)
        self._process_options(options)
        self._set_timeout(timeout)
        self._check_namelist()

    def _process_options(self, options):
        """
        """
        for opt in options:
            name = opt
            value = options[opt]
            processed = False
            try:
                # check for non-tolerance options first, e.g. restart
                # FIXME(bja, 201603) should remove processed options
                # from the dict!
                if name == 'timeout':
                    self._set_timeout(value)
                    processed = True
                else:
                    # comparison will raise an exception if the name is unknown
                    self._comparison.update_from_name(name, value)
                    processed = True
            except RuntimeError as e:
                msg = '{0} : {1}'.format(self._name, e)
                raise RuntimeError(msg)

            if not processed:
                msg = '{0} : Unsupported test option "{1}"'.format(
                    self._name, name)
                raise RuntimeError(msg)

    def _set_timeout(self, timeout_data):
        """Set the timeout, maximum test run time, from the user specified
        value.

        precidence: default < test specific < command line
        """
        timeout = None
        units = None
        if timeout_data:
            data = timeout_data.split()
            if len(data) == 2:
                timeout = float(data[0])
                units = data[1]
            elif len(data) == 1:
                timeout = float(data[0])
            else:
                msg = ('timout data must be a single value in seconds '
                       'or a time followed by units.')
                raise RuntimeError(msg)
            if units:
                if units.lower()[0] == 's':
                    # seconds, do nothing
                    pass
                elif units.lower()[0] == 'm':
                    timeout = timeout * 60.0
                else:
                    msg = ('unknown timeout units "{0}". Valid units '
                           'are seconds and minutes'.format(units))
                    raise RuntimeError(msg)

            self._timeout = timeout

    def _check_namelist(self):
        """Check for the existance of the namelist file for the test.

        """
        namelist_filename = '{0}.namelist'.format(self._name)
        if not os.path.isfile(namelist_filename):
            logging.critical(
                'Missing namelist file for test "{0}"'.format(self._name))
            logging.critical(
                'expected "{0}" in directory : {1}'.format(namelist_filename,
                                                           os.getcwd()))
            self._status = 'skip'
        else:
            self._namelist_filename = namelist_filename

    def status(self):
        """
        """
        return self._status

    def run(self, executable, check_only, dry_run):
        """
        """
        if not check_only:
            self._run_test(executable, dry_run)
        if not self._status:
            self._check_test()

        if self._status == 'skip':
            print('s', end='')
        elif self._status == 'fail':
            print('F', end='')
        elif self._status == 'pass':
            print('.', end='')
        else:
            msg = ('  DEV_ERROR : {0} : unknown status type "{1}".'.format(
                self._name, self._status))
            logging.critical(msg)
            raise RuntimeError(msg)
        sys.stdout.flush()

    def _run_test(self, executable, dry_run):
        """Run the test
        """
        logging.info(SEPERATOR)
        logging.info('Running test "{0}"'.format(self._name))
        self._cleanup_previous_run()
        if self._status == 'skip':
            logging.critical('Skipping test "{0}"').format(self._name)
            return
        cmd = []
        cmd.append(executable)
        cmd.append(self._namelist_filename)

        logging.info("    cd {0}".format(os.getcwd()))
        logging.info("    {0}".format(" ".join(cmd)))
        if not dry_run:
            run_stdout = open(self._name + ".stdout", 'w')
            start = time.time()
            proc = subprocess.Popen(cmd,
                                    shell=False,
                                    stdout=run_stdout,
                                    stderr=subprocess.STDOUT)
            while proc.poll() is None:
                time.sleep(0.1)
                if time.time() - start > self._timeout:
                    proc.kill()
                    time.sleep(0.1)
                    msg = ('    FAILURE: "{0}" exceeded max run time '
                           '{1} seconds.'.format(self._name, self._timeout))
                    logging.critical(''.join(['\n', msg, '\n']))
                    self._status = 'fail'
            finish = time.time()
            logging.info("    {0} : run time : {1:.2f} seconds".format(
                self._name, finish - start))
            run_stdout.close()
            status = abs(proc.returncode)
            if status != 0:
                self._status = 'fail'
                logging.critical('    FAILURE: runtime error in "{0}". '
                                 'See {0}.stdout file for details.'.format(
                                     self._name))

    def _check_test(self):
        """Check the test results against the baseline
        """
        logging.info(SEPERATOR)
        logging.info('Checking test "{0}"'.format(self._name))

        # NOTE(bja, 201603) relying on implicit ordering. We check
        # baseline first, so if the current regression file is missing
        # we get a failure instead of a skip status!
        filename = '{0}.{1}'.format(self._name, self._BASELINE)
        baseline = None
        try:
            baseline = read_config_file_as_dict(filename)
        except RuntimeError as e:
            logging.critical(e)
            msg = ('  SKIP : Could not open baseline file for '
                   'test "{0}".'.format(self._name))
            logging.critical(msg)
            self._status = 'skip'

        filename = '{0}.{1}'.format(self._name, self._REGRESSION)
        regression = None
        try:
            regression = read_config_file_as_dict(filename)
        except RuntimeError as e:
            logging.critical(e)
            msg = ('  FAILURE : Could not open regression file for '
                   'test "{0}".'.format(self._name))
            logging.critical(msg)
            self._status = 'fail'

        if not regression or not baseline:
            return

        # TODO(bja, 201603) is it safe to directly assign status here?
        self._status = self._comparison.regression_to_baseline(
            regression, baseline)

    def update_baseline(self, dry_run):
        """Update the baseline regression file with the results from the
        current run.

        """
        current = '{0}.{1}'.format(self._name, self._REGRESSION)
        baseline = '{0}.{1}'.format(self._name, self._BASELINE)
        try:
            logging.info('  mv {0} {1}'.format(current, baseline))
            if not dry_run:
                os.rename(current, baseline)
                logging.critical('  {0} : updated baseline with current '
                                 'regression file.'.format(self._name))
        except OSError as e:
            msg = '    ERROR: {0} : {1} :\n  moving "{2}" to "{3}"'.format(
                self._name, e, current, baseline)
            logging.critical(msg)

    def _cleanup_previous_run(self):
        """Cleanup any output files that may be left over from a previous run
        so we are sure we test current results.

        """
        output_files = [
            '{0}.{1}'.format(self._name, self._REGRESSION),
            '{0}.stdout'.format(self._name),
            '{0}.output.nc'.format(self._name),
        ]
        for output in output_files:
            src = output
            dest = '{0}.bak'.format(src)
            if os.path.exists(src):
                os.rename(src, dest)


class Comparison(object):
    """Regression test class to isolate the comparison of simulation
    results to baselines. Doesn't have much data of it's own except
    the tolerances, but it make it much easier to see the comparison
    code and unit test if it is in it's own class.

    """

    def __init__(self, name, conf):
        """
        """
        self._name = name
        self._tolerances = Tolerances()
        self._tolerances.update(conf)
        self._status = None

    def update(self, conf):
        """Update from conf dict.

        Only comparison options are tolerances, so we just pass
        through.

        """
        self._tolerances.update(conf)

    def update_from_name(self, name, value):
        """Update object configuration.

        Only comparison options are tolerances, so we just pass
        through.

        """
        self._tolerances.update_from_name(name, value)

    def set_name(self, name):
        """Update the name of the test/suite this object is associated with.

        """
        self._name = name

    def regression_to_baseline(self, regression, baseline):
        """Comparing regression to baseline and baseline and set the status as
        pass if we didn't detect any errors..

        """
        self._compare_sections(baseline, 'baseline',
                               regression, 'regression')

        if not self._status:
            msg = ('  PASS : {0} '.format(self._name))
            logging.critical(msg)
            self._status = 'pass'
        return self._status

    def _compare_sections(self, a_data, a_name, b_data, b_name):
        """Compare the sections. Any sections in a or b but not both are a
        failure.

        """
        for section in a_data:
            if section not in b_data:
                msg = ('  FAILURE : section "{0}" present in {1} '
                       'but missing from {2}'.format(section, a_name, b_name))
                logging.critical(msg)
                self._status = 'fail'
            else:
                # in both a and b
                self._compare_options(
                    section,
                    a_data[section], a_name,
                    b_data[section], b_name)

        for section in b_data:
            # NOTE(bja, 201604) don't need the else clause here
            # because we already checked common sections.
            if section not in a_data:
                msg = ('  FAILURE : section "{0}" present in {1} '
                       'but missing from {2}'.format(section, b_name, a_name))
                logging.critical(msg)
                self._status = 'fail'


    def _compare_options(self, section, a_data, a_name, b_data, b_name):
        """
        """
        a_category = self._get_section_category(section, a_data)
        b_category = self._get_section_category(section, b_data)
        if a_category != b_category:
            self._status = 'fail'
            msg = ('  FAILURE : {0} : data category does not match:\n'
                   '      {1} = {2}\n'
                   '      {3} = {4}\n'.format(section, a_name, a_category,
                                              b_name, b_category))
            logging.critical(msg)

        for key in a_data:
            if key not in b_data:
                msg = ('  FAIURE : {0} : key "{1}" present in {2} '
                       'mising in {3}'.format(section, key, a_name, b_name))
                logging.critical(msg)
                self._status = 'fail'
            elif key == 'category':
                pass
            else:
                self._compare_values_with_tolerance(
                    a_category, section, key, a_data[key], b_data[key])

        for key in b_data:
            # NOTE(bja, 2016) don't need the else to call
            # compare_values again because we already compared
            # sections in both files!
            if key not in a_data:
                msg = ('  FAIURE : {0} : key "{1}" present in {2} '
                       'mising in {3}'.format(section, key, b_name, a_name))
                logging.critical(msg)
                self._status = 'fail'

    def _compare_values_with_tolerance(
            self, category, section, key, a_data, b_data):
        """
        """
        if category is "discrete":
            self._compare_integer_values_with_tolerance(
                category, section, key, a_data, b_data)
        else:
            self._compare_float_values_with_tolerance(
                category, section, key, a_data, b_data)

    def _compare_integer_values_with_tolerance(
            self, category, section, key, a_data, b_data):
        """
        """
        msg = ('  DEV_ERROR : {0} : {1} : comparison of "discrete" '
               'data has not been implemented!'.format(self._name, section))
        raise RuntimeError(msg)

    def _compare_float_values_with_tolerance(
            self, category, section, key, a_data, b_data):
        """
        """
        tol_type = self._tolerances.get(category, 'type')
        tol_value = self._tolerances.get(category, 'value')
        a_value = float(a_data)
        b_value = float(b_data)
        abs_diff = abs(a_value - b_value)
        denominator = self._set_denominator(a_value, b_value)

        if tol_type == Tolerances.ABSOLUTE:
            diff = abs_diff
        elif tol_type == Tolerances.RELATIVE:
            diff = abs_diff / denominator
        elif tol_type == Tolerances.PERCENT:
            diff = 100.0 * abs_diff / denominator
        else:
            # shouldn't be possible to get here if previous error
            # checking was good....
            msg = ('  DEV_ERROR : {0} : {1} : {2} : invalid tolerance '
                   'type "{3}".'.format(self._name, section, key, tol_type))
            raise RuntimeError(msg)

        pass_comparison = None
        if diff > tol_value:
            msg = ('  FAILURE : {0} : {1} : {2} : '
                   'tolerance failure : {3} > {4} [{5}]'.format(
                       self._name, section, key, diff, tol_value, tol_type))
            logging.critical(msg)
            # FIXME(bja, 201603) ugly side-effect, should be setting
            # based on return value.....
            self._status = 'fail'
            pass_comparison = False
        else:
            pass_comparison = True
        return pass_comparison

    def _set_denominator(self, a_value, b_value):
        """Set the denominator, accounting for zero values being a valid
        value, to avoid a floating point error.

        """
        if a_value == 0.0 and b_value != 0.0:
            denominator = b_value
        elif a_value != 0.0 and b_value == 0.0:
            denominator = a_value
        elif a_value == 0 and b_value == 0:
            # set to 1.0 to avoid a floating point exception without
            # special logic.
            denominator = 1.0
        else:  # a != 0 and b != 0
            denominator = a_value

        return denominator

    def _get_section_category(self, section, a):
        """Extract the 'category' value from the section.
        """
        category = None
        if 'category' not in a:
            msg = ('  DEV_ERROR : {0} : section "{1}" missing required '
                   'option "category"'.format(self._name, section))
            logging.critical(msg)
            raise RuntimeError(msg)

        category = a['category']
        self._tolerances.check_valid_category(category)
        return category


class Tolerances(object):
    """Class to manage tolerances for tests
    """
    # tolerance category names
    GENERAL = 'general'
    CONC = 'concentration'
    PRESSURE = 'pressure'
    FRACTION = 'pfraction'
    VELOCITY = 'velocity'
    # DISCRETE = 'discrete'
    _KNOWN_CATEGORIES = [
        GENERAL,
        CONC,
        PRESSURE,
        FRACTION,
        VELOCITY,
    ]

    # tolerance types
    ABSOLUTE = 'absolute'
    RELATIVE = 'relative'
    PERCENT = 'percent'

    _DEFAULT_EPSILON = 1.0e-16

    def __init__(self):
        """
        """
        self._tolerances = {}
        for tol in self._KNOWN_CATEGORIES:
            self._tolerances[tol] = {'value': self._DEFAULT_EPSILON,
                                     'type': self.ABSOLUTE,
                                     'min': 0.0,
                                     'max': sys.float_info.max,
                                     }

    def get(self, category, key):
        """
        """
        if category not in self._tolerances:
            msg = ('  DEV_ERROR : invalid tolerance category "{0}"'
                   ' requested after checking should be complete.')
            logging.critical(msg)
            raise RuntimeError(msg)
        if key not in self._tolerances[category]:
            msg = ('  DEV_ERROR : invalid tolerance data type "{0}"'
                   ' requested.')
            logging.critical(msg)
            raise RuntimeError(msg)

        return self._tolerances[category][key]

    def update(self, conf):
        """Add tolerances from a configuration parser object
        """
        if 'default_tolerances' not in conf:
            # It is legitimate not to specify default tolerances, so
            # we just return.
            return
        new_tolerances = conf['default_tolerances']
        for tol in new_tolerances:
            self.update_from_name(tol, new_tolerances[tol])

        # remove the default tolerances section so we can process the
        # rest of the dict as tests.
        del conf['default_tolerances']

    def update_from_name(self, category, data):
        """
        """
        if category in self._tolerances:
            value, tol_type = data.split()
            try:
                value = float(value)
            except ValueError as e:
                msg = ('tolerance value "{0} = {1}" could not be converted '
                       'to a floaitng point number.'.format(category, value))
                logging.critical(e)
                logging.critical(msg)
                raise RuntimeError(msg)

            if value < self._tolerances[category]['min']:
                msg = ('tolerance "{0} = {1} [{2}]" ; must be greater '
                       'than {3}'.format(
                           category, value, tol_type,
                           self._tolerances[category]['min']))
                logging.critical(msg)
                raise RuntimeError(msg)

            if tol_type == self.PERCENT:
                self._tolerances[category]['max'] = 100.0

            if value > self._tolerances[category]['max']:
                msg = ('tolerance "{0} = {1} [{2}]" ; must be less '
                       'than {3}'.format(
                           category, value, tol_type,
                           self._tolerances[category]['max']))
                logging.critical(msg)
                raise RuntimeError(msg)

            self._tolerances[category]['type'] = tol_type
            self._tolerances[category]['value'] = value
        else:
            msg = "Unknown tolerance type '{0}'.".format(category)
            raise RuntimeError(msg)

    def check_valid_category(self, category):
        """Check that an externally supplied category is one of the valid
        known types.

        """
        if category not in self._KNOWN_CATEGORIES:
            msg = ('  ERROR : invalid tolerance category "{0}". '
                   'Valid values are :\n    {1}'.format(
                       category, self._KNOWN_CATEGORIES))
            raise RuntimeError(msg)
        return True


class StatusSummary(object):
    """
    """

    def __init__(self):
        """
        """
        self._skips = 0
        self._failures = 0
        self._passes = 0
        self._total = 0

    def __add__(self, other):
        """Addition operator to combine status objects

        NOTE(bja, 201603) do NOT manually increment totals here. they
        are incremented automatically.

        """
        new = StatusSummary()
        new.add_skip(self.skips() + other.skips())
        new.add_failure(self.failures() + other.failures())
        new.add_pass(self.passes() + other.passes())
        return new

    def skips(self):
        """
        """
        return self._skips

    def add_skip(self, num=1):
        """
        """
        self._skips += num
        self._total += num

    def failures(self):
        """
        """
        return self._failures

    def add_failure(self, num=1):
        """
        """
        self._failures += num
        self._total += num

    def passes(self):
        """
        """
        return self._passes

    def add_pass(self, num=1):
        """
        """
        self._passes += num
        self._total += num

    def total(self):
        """
        """
        return self._total


# -----------------------------------------------------------------------------
#
# work functions
#
# -----------------------------------------------------------------------------
def find_all_config_files(current_dir):
    """walk the directory and look for cfg files that can be processed for
    tests

    """
    logging.info(SEPERATOR)
    logging.info('Searching for config files.')
    cfg_files = []
    for root, dirs, files in os.walk(current_dir):
        for name in files:
            if name.endswith('.cfg'):
                cfg_files.append(os.path.join(root, name))
    return cfg_files


def append_command_to_log(command, tempfile):
    """Run a shell command and append output to the log file
    """
    logging.info("$ {0}".format(" ".join(command)))
    with open(tempfile, "w") as tempinfo:
        subprocess.call(command, shell=False,
                        stdout=tempinfo,
                        stderr=subprocess.STDOUT)
        # NOTE(bja, 2013-06) : need a short sleep to ensure the
        # contents get written...?
        time.sleep(0.01)
    if os.path.isfile(tempfile):
        with open(tempfile, 'r') as tempinfo:
            for line in tempinfo.readlines():
                logging.info('    {0}'.format(line.strip()))


def setup_log(test_dir):
    """
    """
    time_stamp = datetime.now().strftime('%Y-%m-%d_%H%M%S')
    # log_filename = 'betr-tests-{0}.log'.format(time_stamp)
    log_filename = 'betr-tests.testlog'
    if os.path.isfile(log_filename):
        backup = '{0}.bak'.format(log_filename)
        os.rename(log_filename, backup)
    logging.basicConfig(filename=log_filename,
                        filemode='w',
                        level=logging.INFO,
                        format='%(message)s')
    logging.info("BeTR regression tests")
    logging.info("Run time : {0}".format(time_stamp))
    logging.info("System information :")
    logging.info("    {0}".format(os.uname()))
    logging.info("Test directory :")
    logging.info("    {0}".format(os.getcwd()))

    tempfile = "{0}/tmp-betr-test-info.txt".format(test_dir)

    logging.info(SEPERATOR)
    logging.info("BeTR repository status ")
    if os.path.isdir("{0}/../.git".format(test_dir)):
        cmd = ["git", "log", "-n", "1"]
        append_command_to_log(cmd, tempfile)
        cmd = ["git", "status"]
        append_command_to_log(cmd, tempfile)
    else:
        logging.info("    unknown")

    os.remove(tempfile)


def verify_executable(executable):
    """Check that the user specified executable exists and has the
    executable bit set.

    """
    logging.info(SEPERATOR)
    logging.info('Verifying executable:')
    logging.info('    {0}'.format(executable))

    if not os.path.isfile(executable):
        msg = 'Specified executable does not exist :\n    {0}'.format(
            executable)
        raise RuntimeError(msg)

    if not os.access(executable, os.X_OK):
        msg = 'Can not run specified executable :\n    {0}'.format(
            executable)
        raise RuntimeError(msg)

    logging.info('Executable exists and is can be run by this user.')


def config_to_dict(config):
    """Convert a config object into a dictionary.

    Note: we can probably just use config._sections, but:

      1) that is relying on a a private variable

      2) won't do cfg value interpolation, which we hopfully wouldn't
    use anyway....

    """
    config_dict = {}
    for section in config.sections():
        config_dict[section] = {}
        for option in config.options(section):
            value = config.get(section, option)
            config_dict[section][option] = value

    return config_dict


def convert_input_data(ncgen, input_dir, dry_run):
    """Check that text cdl-netcdf input files are available as binary
    netcdf by calling the external ncgen command.

    FIXME(bja, 201604) assumes that all input files are in a single
    hard coded directory. Should this be genralized....?

    """
    logging.info(SEPERATOR)
    logging.info('Converting input data from text cdl to binary nc.')
    all_files = os.listdir(input_dir)
    nc_files = [f for f in all_files if f.endswith('.cdl.nc')]
    cdl_files = [f for f in all_files if f.endswith('.cdl')]

    tempfile = "{0}/tmp-betr-test-info.txt".format(input_dir)

    for cdl in cdl_files:
        nc = cdl + '.nc'
        nc_filename = os.path.join(input_dir, nc)
        cdl_filename = os.path.join(input_dir, cdl)

        convert = False
        if nc not in nc_files:
            # nc file doesn't exist, need to convert
            convert = True
            logging.info('  {0}: nc doesn\'t exist.'.format(cdl))
        else:
            nc_timestamp = os.path.getmtime(nc_filename)
            cdl_timestamp = os.path.getmtime(cdl_filename)
            mod_time = cdl_timestamp - nc_timestamp
            if mod_time > 0:
                # cdl file in newer than nc file. need to convert
                convert = True
                logging.info(
                    '  {0}: cdl file is newer than nc file: {1}.'.format(
                        cdl, mod_time))

        if convert:
            msg = '  converting {0} to binary nc format'.format(cdl)
            logging.info(msg)
            command = [
                ncgen,
                '-o',
                nc_filename,
                cdl_filename,
            ]
            append_command_to_log(command, tempfile)
            os.remove(tempfile)
            logging.info('')


def find_executable(exe, paths):
    """Check if an executable exists is a list of paths and use it, or try
    to see if it is just in the user path.

    """
    logging.info('Searching for executable "{0}"'.format(exe))
    found = False
    full_path = None
    for path in paths:
        full_path = os.path.join(path, exe)
        try:
            verify_executable(full_path)
        except RuntimeError as e:
            # not useable, move on
            logging.info('    checked : "{0}"'.format(full_path))
            pass
        else:
            found = True
            break

    if not found:
        # check to see if it is just in the user path.
        try:
            command = [exe, '--help']
            subprocess.call(command, shell=False)
        except OSError:
            # not in path:
            pass
        else:
            full_path = exe
            found = True

    if not found:
        msg = ("ERROR: could not find a usable executable for '{0}' in:\n"
               "    {1}".format(exe, paths))
        raise RuntimeError(msg)
    logging.info('  Using "{0}" at :\n    {1}'.format(exe, full_path))
    return full_path

# -----------------------------------------------------------------------------
#
# main
#
# -----------------------------------------------------------------------------
def main(options):
    print(BANNER)
    print('BeTR regression test driver')
    start_time = time.time()

    # FIXME(bja, 201604) assuming that we are always being called from
    # the same directory. Is there a need to support calling from an
    # arbitrary location, e.g. integration with ctest?
    cwd = os.getcwd()
    setup_log(cwd)
    executable = os.path.abspath(options.executable[0])
    verify_executable(executable)
    exe_paths = [os.path.dirname(executable), '']

    if options.config:
        filenames = options.config
    else:
        test_root = os.path.join(cwd, 'tests')
        filenames = find_all_config_files(test_root)

    input_dir = os.path.join(cwd, 'input-data')
    if not os.path.isdir(input_dir):
        msg = ('  DEV_ERROR: input data directory does not exist:'
               '\n    {0}'.format(input_dir))
        raise RuntimeError(msg)

    dry_run = options.dry_run
    ncgen = find_executable('ncgen', exe_paths)
    convert_input_data(ncgen, input_dir, dry_run)

    print('Setting up tests.')
    test_suites = []
    for filename in filenames:
        os.chdir(cwd)
        logging.info(BANNER)
        config = read_config_file_as_dict(filename)
        suite = RegressionTestSuite(filename, config, options.timeout)
        if suite.num_tests() > 0:
            test_suites.append(suite)
        else:
            logging.info('   WARNING: No tests in suite {0}'.format(filename))
            del suite

    print('Running tests:')
    check_only = options.check_only
    update_baseline = options.update_baseline
    for suite in test_suites:
        suite.run_tests(executable, check_only, dry_run, update_baseline)

    print()
    print(BANNER)

    status = StatusSummary()
    for suite in test_suites:
        status += suite.status_summary()

    end_time = time.time()

    print("Status:")
    print("  total tests : {0}".format(status.total()))
    print("    skipped : {0}".format(status.skips()))
    print("    failed : {0}".format(status.failures()))
    print("    passed : {0}".format(status.passes()))
    print("  overall time : {0:5.2f} [s]".format(end_time - start_time))

    logging.shutdown()
    return status.total() - status.passes()


if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        if options.backtrace:
            traceback.print_exc()
        logging.shutdown()
        sys.exit(1)
