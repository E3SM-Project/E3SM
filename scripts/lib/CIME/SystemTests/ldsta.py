"""
CIME last date short term archiver test. This class inherits from SystemTestsCommon
It does a run without restarting, then runs the archiver with various last-date parameters
The test verifies the archive directory contains the expected files
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import expect
from CIME.case_st_archive import case_st_archive
from CIME.case_st_archive import _get_file_date

import datetime
import glob
import random

logger = logging.getLogger(__name__)

class LDSTA(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

    def run_phase(self):
        self.run_indv()
        # finished running, so all archive files should exist
        start_date = datetime.datetime(1, 1, 1)
        last_date = datetime.datetime(2, 1, 1)
        rest_dir = os.path.join(self._case.get_value('DOUT_S_ROOT'), 'rest')
        current_date = self._increment_month(start_date)
        next_datecheck = current_date
        while current_date < last_date:
            logger.info('Testing archiving with last date: {}'.format(current_date))
            current_date_str = '{:04}-{:02}-{:02}'.format(current_date.year,
                                                          current_date.month,
                                                          current_date.day)
            case_st_archive(self._case, last_date_str=current_date_str, copy_only=False)
            archive_dates = [_get_file_date(fname)
                             for fname in glob.glob(os.path.join(rest_dir, '*'))]
            while next_datecheck <= current_date:
                expect(next_datecheck in archive_dates,
                       'Not all dates generated and/or archived: '
                       + '{} is missing'.format(next_datecheck))
                next_datecheck = self._increment_month(next_datecheck)
            for date in archive_dates:
                expect(date <= current_date,
                       'Archived date greater than specified by last-date: '
                       + '{}'.format(date))
            num_months = random.randint(1, 4)
            for i in range(num_months):
                current_date = self._increment_month(current_date)

    @staticmethod
    def _increment_month(date):
        if date.month < 12:
            return datetime.datetime(date.year, date.month + 1, date.day)
        else:
            return datetime.datetime(date.year + 1, 1, date.day)
