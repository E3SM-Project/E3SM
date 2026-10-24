"""Run with python3; validate CIME calendar configurations without model cases."""

import importlib.machinery
import importlib.util
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "cime"))
sys.path.insert(0, str(ROOT / "cime" / "CIME" / "Tools"))
os.environ.setdefault("CIME_MODEL", "e3sm")

from CIME.nmlgen import NamelistGenerator
from CIME.get_tests import get_test_data
from CIME.SystemTests.test_mods import find_test_mods
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import get_full_test_name, parse_test_name


class CalendarConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.drivers = {}
        for driver in ("driver-mct", "driver-moab"):
            path = ROOT / driver / "cime_config" / "buildnml"
            loader = importlib.machinery.SourceFileLoader(driver, str(path))
            spec = importlib.util.spec_from_loader(loader.name, loader)
            module = importlib.util.module_from_spec(spec)
            loader.exec_module(module)
            cls.drivers[driver] = module

    def generator(self, driver, calendar, profile=None):
        case = Mock()
        case.get_value.side_effect = {"DIN_LOC_ROOT": str(ROOT)}.get
        definition = ROOT / driver / "cime_config" / "namelist_definition_drv.xml"
        generator = NamelistGenerator(case, [str(definition)])
        infiles = []
        if profile is not None:
            infiles.append(str(ROOT / "cime_config" / "testmods_dirs" / "allactive" /
                               ("calendar_" + profile) / "user_nl_cpl"))
        generator.init_defaults(infiles, {"calendar": calendar}, skip_entry_loop=True)
        for name in ("orb_mode", "orb_iyear", "orb_iyear_align", "orb_use_equation_of_time"):
            generator.add_default(name)
        return generator

    def test_epoch_defaults_and_eot_default(self):
        for driver, module in self.drivers.items():
            for calendar in ("NO_LEAP", "GREGORIAN"):
                for token, expected in (("1850", "1850"), ("2010", "2010"), ("20TR", "1850")):
                    with self.subTest(driver=driver, calendar=calendar, token=token):
                        generator = self.generator(driver, calendar)
                        module._set_orbital_defaults(generator, token)
                        self.assertEqual(generator.get_value("orb_iyear"), expected)
                        self.assertEqual(generator.get_value("orb_iyear_align"), expected)
                        self.assertEqual(generator.get_value("orb_mode"), "fixed_year")
                        self.assertEqual(generator.get_value("orb_use_equation_of_time"), ".true.")

    def test_profiles(self):
        for driver, module in self.drivers.items():
            for profile in ("nf0", "nf1", "nv0", "nv1", "gf0", "gf1", "gv0", "gv1"):
                with self.subTest(driver=driver, profile=profile):
                    calendar = "GREGORIAN" if profile.startswith("g") else "NO_LEAP"
                    generator = self.generator(driver, calendar, profile)
                    module._set_orbital_defaults(generator, "1850")
                    variable = profile[1] == "v"
                    self.assertEqual(generator.get_value("orb_mode"),
                                     "variable_orbit_year" if variable else "fixed_year")
                    self.assertEqual(generator.get_value("orb_iyear"), "1999" if variable else "1850")
                    self.assertEqual(generator.get_value("orb_iyear_align"), "1999" if variable else "1850")
                    self.assertEqual(generator.get_value("orb_use_equation_of_time"),
                                     ".true." if profile.endswith("1") else ".false.")

    def test_explicit_epoch_and_legacy_alias(self):
        for driver, module in self.drivers.items():
            for calendar in ("NO_LEAP", "GREGORIAN"):
                with self.subTest(driver=driver, calendar=calendar):
                    generator = self.generator(driver, calendar)
                    generator.set_value("orb_mode", "variable_year")
                    generator.set_value("orb_iyear", "1950")
                    generator.set_value("orb_iyear_align", "2000")
                    module._set_orbital_defaults(generator, "1850")
                    self.assertEqual(generator.get_value("orb_mode"), "variable_year")
                    self.assertEqual(generator.get_value("orb_iyear"), "1950")
                    self.assertEqual(generator.get_value("orb_iyear_align"), "2000")

    def test_suite_modifications(self):
        for suffix, expected in (("smoke", 16), ("short", 128), ("cycle", 64)):
            tests = get_test_data("e3sm_calendar_" + suffix)[-1]
            self.assertEqual(len(tests), expected)
            self.assertEqual(len(set(tests)), expected)
            for test in tests:
                with self.subTest(suite=suffix, test=test):
                    test_name, modifiers = test.rsplit(".", 1)
                    full_name = get_full_test_name(
                        test_name, machine="pm-cpu", compiler="gnu", testmods_string=modifiers)
                    parsed = parse_test_name(full_name)
                    paths = find_test_mods("mct", parsed[6])
                    self.assertEqual(len(paths), 2)
                    for path in paths:
                        self.assertTrue((Path(path) / "shell_commands").is_file())

    def test_restart_pointers(self):
        for calendar in ("GREGORIAN", "NO_LEAP"):
            cycle_restart = ".2000-12-31-00000" if calendar == "GREGORIAN" else ".2001-01-01-00000"
            for start, stop, interval, expected in (
                    ("1999-01-01", 1461 if calendar == "GREGORIAN" else 1460, 730, cycle_restart),
                    ("2000-02-27", 6, 1, ".2000-02-28-00000"),
                    ("1999-12-29", 6, 1, ".1999-12-30-00000")):
                with self.subTest(calendar=calendar, start=start):
                    values = dict(STOP_N=stop, STOP_OPTION="ndays", COMP_ATM="datm",
                                  ATM_NCPL=24, COMP_INTERFACE="mct", NCPL_BASE_PERIOD="day",
                                  CALENDAR=calendar, START_TOD=0, RUN_STARTDATE=start)
                    case = Mock()
                    case.get_value.side_effect = values.get
                    case.get_values.return_value = ["CPL", "ATM"]
                    context = SimpleNamespace(
                        _case=case, _rest_n=interval,
                        _leap_year_correction=SystemTestsCommon._leap_year_correction)
                    SystemTestsCommon._set_restart_interval(context)
                    self.assertEqual(context._rest_time, expected)


if __name__ == "__main__":
    unittest.main()