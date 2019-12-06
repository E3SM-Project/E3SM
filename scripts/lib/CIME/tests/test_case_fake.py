#/usr/bin/env python

"""
This module contains unit tests of CaseFake
"""

import unittest
import tempfile
import os
import shutil

from CIME.tests.case_fake import CaseFake

class TestCaseFake(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_create_clone(self):
        # Setup
        old_caseroot = os.path.join(self.tempdir, 'oldcase')
        oldcase = CaseFake(old_caseroot)
        oldcase.set_value('foo', 'bar')

        # Exercise
        new_caseroot = os.path.join(self.tempdir, 'newcase')
        clone = oldcase.create_clone(new_caseroot)

        # Verify
        self.assertEqual('bar', clone.get_value('foo'))
        self.assertEqual('newcase', clone.get_value('CASE'))
        self.assertEqual('newcase', clone.get_value('CASEBASEID'))
        self.assertEqual(new_caseroot, clone.get_value('CASEROOT'))
        self.assertEqual(os.path.join(new_caseroot, 'run'),
                         clone.get_value('RUNDIR'))
