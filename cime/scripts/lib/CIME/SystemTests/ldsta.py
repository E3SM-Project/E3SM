"""
CIME last date short term archiver test. This class inherits from SystemTestsCommon
It does a run without restarting, then runs the archiver with various last-date parameters
The test verifies the archive directory contains the expected files
"""

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import expect
from CIME.date import get_file_date

import datetime
import glob
import os
import random
import shutil

logger = logging.getLogger(__name__)

# datetime objects can't be used anywhere else
def _date_to_datetime(date_obj):
    return datetime.datetime(year = date_obj.year(),
                             month = date_obj.month(),
                             day = date_obj.day(),
                             hour = date_obj.hour(),
                             minute = date_obj.minute(),
                             second = date_obj.second())

class LDSTA(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        SystemTestsCommon.__init__(self, case)

    def run_phase(self):
        archive_dir = self._case.get_value('DOUT_S_ROOT')
        if os.path.isdir(archive_dir):
            shutil.rmtree(archive_dir)
        self.run_indv()
        # finished running, so all archive files should exist
        start_date = _date_to_datetime(get_file_date(self._case.get_value('RUN_STARTDATE')))
        rest_dir = os.path.join(archive_dir, 'rest')
        delta_day = datetime.timedelta(1)
        current_date = start_date + delta_day
        next_datecheck = current_date
        days_left = self._case.get_value('STOP_N')
        final_date = start_date + delta_day * days_left
        while current_date < final_date:
            logger.info('Testing archiving with last date: {}'.format(current_date))
            current_date_str = '{:04}-{:02}-{:02}'.format(current_date.year,
                                                          current_date.month,
                                                          current_date.day)
            self._case.case_st_archive(last_date_str=current_date_str, copy_only=False)
            archive_dates = [_date_to_datetime(get_file_date(fname))
                             for fname in glob.glob(os.path.join(rest_dir, '*'))]
            while next_datecheck <= current_date:
                expect(next_datecheck in archive_dates,
                       'Not all dates generated and/or archived: '
                       + '{} is missing'.format(next_datecheck))
                next_datecheck += delta_day
            for date in archive_dates:
                expect(date <= current_date,
                       'Archived date greater than specified by last-date: '
                       + '{}'.format(date))
            num_days = random.randint(1, min(3, days_left))
            days_left -= num_days
            current_date += num_days * delta_day
