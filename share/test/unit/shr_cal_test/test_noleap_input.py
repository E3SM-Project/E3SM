"""Local fixture/harness checks; production reader tests require model libraries."""

from datetime import datetime
from contextlib import redirect_stdout
import gzip
from io import StringIO
import json
from importlib.util import find_spec
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "cime"))
sys.path.insert(0, str(ROOT / "cime_config" / "SystemTests"))
sys.path.insert(0, str(ROOT / "cime_config" / "scripts"))
os.environ.setdefault("CIME_MODEL", "e3sm")

from noleap_input_cases import (
    CASES, expected_value, export_snapshot, fresh_logs, log_snapshot, require_expected_failure,
    selected_case, set_user_policy, stage_fixtures, validate_activation, validate_exports,
)


class NoleapInputHarnessTests(unittest.TestCase):
    def test_reverse_case_manifest_and_oracles(self):
        from noleap_input_cases import GREGORIAN_CASES

        for name in GREGORIAN_CASES:
            selected, specification = selected_case("allactive-calendar_gregorian_to_noleap_" + name)
            self.assertEqual(selected, name)
            self.assertEqual(specification["direction"], "gregorian_to_noleap")
            self.assertEqual(specification["model_calendar"], "NO_LEAP")
        records = [("2000-02-28T00:00:00", 10), ("2000-02-29T00:00:00", 100), ("2000-03-01T00:00:00", 14)]
        self.assertEqual(expected_value(records, datetime(2000, 2, 28, 12)), 55)
        self.assertEqual(expected_value(records, datetime(2000, 3, 1)), 14)
        self.assertEqual(expected_value([records[0], records[2]], datetime(2000, 2, 28, 12)), 11)
        subdaily = [("2000-02-28T18:00:00", 10), ("2000-03-01T00:00:00", 14)]
        self.assertAlmostEqual(expected_value(subdaily, datetime(2000, 2, 28, 21)), 10.4)

    def test_reverse_preflight_policy(self):
        from input_calendar_preflight import read_policy

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "user_nl_cpl"
            path.write_text("fudge_gregorian_to_noleap=.true.\n")
            policy = read_policy(path, True, "NO_LEAP")
            self.assertEqual(policy["direction"], "gregorian_to_noleap")
            self.assertTrue(policy["model_calendar_verified"])
            with self.assertRaisesRegex(ValueError, "requires CALENDAR=NO_LEAP"):
                read_policy(path, True, "GREGORIAN")
            with self.assertRaisesRegex(ValueError, "supply --model-calendar"):
                read_policy(path, True)
            path.write_text("fudge_gregorian_to_noleap=.true.\nfudge_noleap_for_gregorian=.true.\n")
            with self.assertRaisesRegex(ValueError, "mutually exclusive"):
                read_policy(path, True, "NO_LEAP")

    @unittest.skipUnless(find_spec("netCDF4"), "netCDF4 is required")
    def test_reverse_preflight_retains_leapday(self):
        import netCDF4
        from input_calendar_preflight import main

        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            policy = folder / "user_nl_cpl"
            policy.write_text("fudge_gregorian_to_noleap=.true.\n")
            path = folder / "forcing.nc"
            with netCDF4.Dataset(path, "w") as dataset:
                dataset.createDimension("time", 3)
                time = dataset.createVariable("time", "f8", ("time",))
                time.calendar = "gregorian"
                time.units = "days since 2000-01-01 00:00:00"
                time[:] = [58, 59, 60]
            output = StringIO()
            with redirect_stdout(output):
                status = main(["--policy-namelist", str(policy), "--groupless", "--model-calendar", "NO_LEAP",
                               "--input", "eamxx", str(path)])
            report = json.loads(output.getvalue())
            self.assertEqual(status, 0, report)
            self.assertEqual(report["inputs"][0]["leap_day_records"], 1)
            self.assertEqual(report["inputs"][0]["interpolation_calendar"], "GREGORIAN")

    @unittest.skipUnless(find_spec("netCDF4"), "netCDF4 is required")
    def test_calendar_preflight_configured_inventory(self):
        import netCDF4
        from input_calendar_preflight import main

        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            policy = folder / "user_nl_cpl"
            policy.write_text("fudge_noleap_for_gregorian=.true.\n"
                              "input_calendar_override_selectors='stream:forcing.txt'\n"
                              "input_calendar_override_calendars='365_day'\n")
            component = folder / "datm_in"
            component.write_text("&unrelated_nml\n unused=1\n/\n"
                                 "&shr_strdata_nml\n streams='forcing.txt 2000 2000 2000'\n unused=2\n/\n")
            (folder / "forcing.txt").write_text(
                "<fieldInfo><filePath>" + directory + "</filePath>"
                "<fileNames>forcing.nc</fileNames></fieldInfo>")
            with netCDF4.Dataset(folder / "forcing.nc", "w") as dataset:
                dataset.createDimension("time", 2)
                time = dataset.createVariable("time", "f8", ("time",))
                time.units = "days since 2000-02-28 00:00:00"
                time[:] = [0, 1]
            output = StringIO()
            with redirect_stdout(output):
                status = main(["--policy-namelist", str(policy), "--groupless",
                               "--rundir", directory, "--namelist", str(component)])
            report = json.loads(output.getvalue())
            self.assertEqual(status, 0, report)
            self.assertFalse(report["inventory_complete"])
            self.assertEqual(report["inputs"][0]["provenance"], "override")
            self.assertEqual(report["inputs"][0]["stream"], "forcing.txt")

    def test_calendar_preflight_resolution(self):
        from input_calendar_preflight import make_rules, normalize_calendar, resolve_calendar

        rules = make_rules(["file:/forcing/a.nc", "stream:forcing.txt"], ["365_day", "noleap"])
        self.assertEqual(normalize_calendar("\tNO_LEAP\0 "), "NO_LEAP")
        self.assertEqual(normalize_calendar("no leap"), "")
        self.assertEqual(resolve_calendar("", False, "/forcing/a.nc", "forcing.txt", rules),
                         ("NO_LEAP", "override", ["file:/forcing/a.nc", "stream:forcing.txt"]))
        with self.assertRaisesRegex(ValueError, "required"):
            resolve_calendar("", False, "a.nc", "", rules)
        with self.assertRaisesRegex(ValueError, "conflicts"):
            resolve_calendar("gregorian", True, "/forcing/a.nc", "", rules)
        with self.assertRaisesRegex(ValueError, "unsupported"):
            resolve_calendar("360_day", True, "/forcing/a.nc", "", rules)

    @unittest.skipUnless(find_spec("netCDF4"), "netCDF4 is required")
    def test_historical_calendar_flavors(self):
        import netCDF4
        from input_calendar_preflight import inspect_input, make_rules, resolve_calendar

        rules = make_rules(["file:a.nc"], ["proleptic_gregorian"])
        with self.assertRaisesRegex(ValueError, "conflicts"):
            resolve_calendar("standard", True, "a.nc", "", rules)
        self.assertNotEqual(rules, make_rules(["file:a.nc"], ["standard"]))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "historical.nc"
            entry = {"reader": "chemistry", "file": str(path)}
            for flavor, epoch, offsets, endpoints, accepted in (
                    ("standard", "0001-01-01", [730181, 730182], None, False),
                    ("gregorian", "2000-01-01", [0, 60], None, True),
                    ("proleptic_gregorian", "0001-01-01", [730181, 730182], None, True),
                    ("standard", "1582-10-15", [-1, 0], None, False),
                    ("standard", "1582-10-15", [0, 1], [[-0.1, 0.1], [0.9, 1.1]], False),
                    ("standard", "1582-10-10", [5, 6], None, False)):
                with self.subTest(flavor=flavor, epoch=epoch, endpoints=endpoints):
                    with netCDF4.Dataset(path, "w") as dataset:
                        dataset.createDimension("time", 2)
                        time = dataset.createVariable("time", "f8", ("time",))
                        time.calendar = flavor
                        time.units = "days since " + epoch + " 00:00:00"
                        time[:] = offsets
                        if endpoints is not None:
                            dataset.createDimension("nv", 2)
                            time.bounds = "time_bnds"
                            bounds = dataset.createVariable("time_bnds", "f8", ("time", "nv"))
                            bounds[:] = endpoints
                    report = inspect_input(entry, {}, Path(directory), "gregorian_to_noleap")
                    self.assertEqual(report["status"], "ok" if accepted else "error", report)
                    if not accepted:
                        self.assertIn("unsupported historical", report["error"])

    @unittest.skipUnless(find_spec("netCDF4"), "netCDF4 is required")
    def test_calendar_preflight_coordinates(self):
        import netCDF4
        from input_calendar_preflight import inspect_input, make_rules

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "forcing.nc"
            entry = {"file": str(path), "reader": "stream", "stream": "forcing.txt"}
            rules = make_rules(["stream:forcing.txt", "file:" + str(path)], ["noleap", "365_day"])
            with netCDF4.Dataset(path, "w") as dataset:
                dataset.createDimension("time", 2)
                time = dataset.createVariable("time", "f8", ("time",))
                time.units = "days since 2000-02-28 00:00:00"
                time[:] = [0, 1]
            report = inspect_input(entry, rules, Path(directory))
            self.assertEqual(report["status"], "ok", report)
            self.assertEqual(report["provenance"], "override")
            self.assertTrue(report["last"].startswith("2000-03-01"))
            with netCDF4.Dataset(path, "a") as dataset:
                dataset.variables["time"][:] = [1, 0]
            self.assertIn("strictly increasing", inspect_input(entry, rules, Path(directory))["error"])
            for offsets in ([0, 0], [0, float("nan")]):
                with netCDF4.Dataset(path, "a") as dataset:
                    dataset.variables["time"][:] = offsets
                self.assertEqual(inspect_input(entry, rules, Path(directory))["status"], "error")
            with netCDF4.Dataset(path, "a") as dataset:
                dataset.variables["time"][:] = [0, 1]
                quantity = dataset.createVariable("quantity", "f8", ("time",))
                quantity.cell_methods = "time: sum"
            chemistry = dict(entry, reader="chemistry")
            self.assertIn("cell_methods", inspect_input(chemistry, rules, Path(directory))["error"])
            with netCDF4.Dataset(path, "a") as dataset:
                dataset.variables["quantity"].cell_methods = "time: mean"
            self.assertIn("require mapped time bounds", inspect_input(chemistry, rules, Path(directory))["error"])
            with netCDF4.Dataset(path, "a") as dataset:
                dataset.createDimension("nv", 2)
                bounds = dataset.createVariable("time_bnds", "f8", ("time", "nv"))
                bounds[:] = [[-0.25, 0.25], [0.75, 1.25]]
            self.assertEqual(inspect_input(chemistry, rules, Path(directory))["status"], "ok")
            unsupported = dict(entry, reader="mpas")
            self.assertIn("unsupported reader", inspect_input(unsupported, rules, Path(directory))["error"])

    def test_positive_history_matches_atmosphere_cadence(self):
        from CIME.SystemTests.ers import ERS
        from nlf import NLF

        for frequency, hours in ((1, 24), (6, 4), (24, 1), (32, 3), (48, 1)):
            with self.subTest(frequency=frequency):
                changes = []
                test = object.__new__(NLF)
                values = {"ATM_NCPL": frequency, "NCPL_BASE_PERIOD": "day",
                          "RUNDIR": ".", "CASE": "calendar_feb"}
                test._case = SimpleNamespace(
                    get_value=values.get, set_value=lambda *args: changes.append(args),
                    create_namelists=lambda: None)
                test._fixture_days = 6
                test._specification = {"enabled": False}
                with patch("nlf.stage_fixtures", return_value={"enabled": False}), \
                        patch("nlf.log_snapshot", return_value={}), \
                        patch("nlf.export_snapshot", return_value={}), \
                        patch("nlf.fresh_logs", return_value=""), \
                        patch("nlf.validate_exports", return_value=(2, 1)), \
                        patch.object(ERS, "run_indv"):
                    test.run_indv()
                self.assertEqual(changes, [("HIST_OPTION", "nhours"), ("HIST_N", hours)])

    def test_each_stream_activation_required(self):
        manifest = {"enabled": True, "calendar": "noleap", "streams": ["first", "second"]}
        logs = ("WARNING NOLEAP_INPUT: fudge_noleap_for_gregorian is ENABLED\n"
                "WARNING NOLEAP_INPUT: stream 1 source=NO_LEAP model=GREGORIAN\n"
                "WARNING NOLEAP_INPUT: description=first 2000 2000 2000\n")
        with self.assertRaisesRegex(AssertionError, "activation for second"):
            validate_activation(manifest, logs)
        validate_activation(manifest, logs + "WARNING NOLEAP_INPUT: description=second\n")
        with self.assertRaisesRegex(AssertionError, "policy banner"):
            validate_activation({"enabled": False}, logs)

    def test_negative_staging_after_each_namelist_generation(self):
        from CIME.SystemTests.ers import ERS
        from nln import NLN

        events = []
        test = object.__new__(NLN)
        test._case = SimpleNamespace(create_namelists=lambda: events.append("namelists"))
        test._fixture_days = 6
        test._specification = {"enabled": True}

        def staged(case, specification, days):
            self.assertIs(case, test._case)
            self.assertEqual(days, 6)
            events.append("fixtures")

        def run(*args, **kwargs):
            self.assertTrue(test._skip_pnl)
            events.append(kwargs.get("suffix", "base"))

        with patch("nln.stage_fixtures", side_effect=staged), patch.object(ERS, "run_indv", side_effect=run):
            test._skip_pnl = False
            test.run_indv()
            test._skip_pnl = False
            test.run_indv(suffix="rest")
        self.assertEqual(events, ["namelists", "fixtures", "base", "namelists", "fixtures", "rest"])

    def test_restart_negative_restores_initial_policy_on_retry(self):
        from CIME.utils import CIMEError
        from nln import NLN

        for name in ("restart_on", "restart_off"):
            with self.subTest(fixture=name):
                test = object.__new__(NLN)
                test._case = SimpleNamespace(get_value={"RUNDIR": "."}.get)
                test._fixture_name, test._specification = selected_case("allactive-calendar_noleap_" + name)
                events = []
                def policy(case, enabled):
                    events.append(enabled)
                with patch("nln.set_user_policy", side_effect=policy), \
                        patch.object(test, "_ers_first_phase", side_effect=lambda: events.append("base")), \
                        patch.object(test, "_ers_second_phase", side_effect=CIMEError("expected abort")), \
                        patch("nln.log_snapshot", return_value={}), \
                        patch("nln.fresh_logs", return_value=test._specification["error"]):
                    test.run_phase()
                self.assertEqual(events, [test._specification["enabled"], "base",
                                          test._specification["restart_enabled"]])

    def test_daily_and_subdaily_oracles(self):
        query = datetime(2000, 2, 29)
        daily = [("2000-02-28T00:00:00", 10), ("2000-03-01T00:00:00", 14)]
        self.assertEqual(expected_value(daily, query), 12)
        self.assertEqual(expected_value(daily, datetime(2000, 2, 29, 12)), 13)
        self.assertEqual(expected_value(daily, query, enabled=False), 10)
        self.assertEqual(expected_value(daily, datetime(2000, 2, 28, 12), enabled=False), 12)
        self.assertEqual(expected_value(daily, datetime(2000, 2, 29, 12), enabled=False), 12)
        subdaily = [("2000-02-28T18:00:00", 10), ("2000-03-01T00:00:00", 14)]
        self.assertAlmostEqual(expected_value(subdaily, query), 10.8)
        self.assertEqual(expected_value(daily, datetime(2000, 3, 1)), 14)
        for outside in (datetime(2000, 2, 27), datetime(2000, 3, 2)):
            with self.assertRaises(ValueError):
                expected_value(daily, outside)

    def test_monthly_and_century_oracles(self):
        for year in (1900, 1999, 2000, 2004, 2100, 2400):
            leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
            records = [(f"{year}-02-01T00:00:00", 10), (f"{year}-03-01T00:00:00", 14)]
            query = datetime(year, 2, 29 if leap else 28)
            expected = 10 + 4 * (28 / 29 if leap else 27 / 28)
            self.assertAlmostEqual(expected_value(records, query), expected)

    def test_specific_rejection_required(self):
        error = RuntimeError("run aborted")
        require_expected_failure(error, "NOLEAP_INPUT: expected", "NOLEAP_INPUT: expected")
        for failure, logs in ((None, "NOLEAP_INPUT: expected"),
                              (error, "missing input data"), (error, "")):
            with self.assertRaises(AssertionError):
                require_expected_failure(failure, "NOLEAP_INPUT: expected", logs)

    def test_stale_and_appended_logs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cpl.log.old"
            path.write_text("NOLEAP_INPUT: old failure\n")
            before = log_snapshot(directory)
            self.assertEqual(fresh_logs(directory, before), "")
            path.write_text(path.read_text() + "new failure\n")
            self.assertEqual(fresh_logs(directory, before), "new failure\n")
            with gzip.open(Path(directory) / "cpl.log.new.gz", "wt") as stream:
                stream.write("current compressed failure\n")
            self.assertIn("current compressed failure", fresh_logs(directory, before))

    def test_case_selection(self):
        for name in CASES:
            full = "NLF_Ld6.grid.A.machine_gnu.allactive-calendar_gf1--allactive-calendar_noleap_" + name
            selected, specification = selected_case(full)
            self.assertEqual(selected, name)
            self.assertIn("enabled", specification)
        with self.assertRaises(ValueError):
            selected_case("ERS_Ld6.grid.A")

    def test_compressing_old_log_is_not_a_new_diagnostic(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cpl.log.old"
            path.write_text("NOLEAP_INPUT: stale failure\n")
            before = log_snapshot(directory)
            compressed = Path(str(path) + ".gz")
            with gzip.open(compressed, "wb") as stream:
                stream.write(path.read_bytes())
            path.unlink()
            self.assertEqual(fresh_logs(directory, before), "")
            before = log_snapshot(directory)
            with gzip.open(compressed, "at") as stream:
                stream.write("current unrelated failure\n")
            self.assertEqual(fresh_logs(directory, before), "current unrelated failure\n")

    def test_policy_update_preserves_other_settings(self):
        from CIME import namelist

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "user_nl_cpl"
            path.write_text("orb_use_equation_of_time = .false.\nfudge_noleap_for_gregorian = .true.\n")
            case = SimpleNamespace(get_value={"CASEROOT": directory}.get)
            set_user_policy(case, False)
            parsed = namelist.parse(in_file=str(path), groupless=True)
            self.assertEqual(parsed["fudge_noleap_for_gregorian"], [".false."])
            self.assertEqual(parsed["orb_use_equation_of_time"], [".false."])
            set_user_policy(case, True, "gregorian_to_noleap")
            parsed = namelist.parse(in_file=str(path), groupless=True)
            self.assertEqual(parsed["fudge_gregorian_to_noleap"], [".true."])
            self.assertEqual(parsed["fudge_noleap_for_gregorian"], [".false."])
            self.assertEqual(parsed["orb_use_equation_of_time"], [".false."])


@unittest.skipUnless(find_spec("netCDF4"), "netCDF4 is required for file-backed fixture checks")
class NoleapInputFixtureTests(unittest.TestCase):
    def setUp(self):
        import netCDF4

        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.rundir = Path(self.directory.name)
        with netCDF4.Dataset(self.rundir / "original.nc", "w") as dataset:
            dataset.createDimension("time", None)
            dataset.createDimension("ncol", 2)
            time = dataset.createVariable("time", "f8", ("time",))
            time[:] = [0]
            time.units = "days since 1999-01-01"
            for name, value in (("temperature", 285), ("humidity", 0.005)):
                variable = dataset.createVariable(name, "f8", ("time", "ncol"))
                variable[:] = value
                variable.cell_methods = "time: mean"
        stream = self.rundir / "datm.streams.txt.test"
        stream.write_text("<stream><fieldInfo><variableNames>\n"
                          "temperature tbot\nhumidity qbot\n</variableNames>"
                          "<filePath>" + str(self.rundir) + "</filePath>"
                          "<fileNames>original.nc</fileNames></fieldInfo></stream>")
        (self.rundir / "datm_in").write_text(
            "&shr_strdata_nml\nstreams='datm.streams.txt.test 1999 1999 2000'\n"
            "tintalgo='coszen'\ntaxmode='cycle'\n/\n")
        self.case = SimpleNamespace(get_value={
            "RUNDIR": str(self.rundir), "RUN_STARTDATE": "2000-02-27"}.get)

    def test_staging_uses_source_calendar_and_point_metadata(self):
        import netCDF4
        from CIME import namelist

        _, specification = selected_case("allactive-calendar_noleap_daily")
        manifest = stage_fixtures(self.case, specification, 6)
        self.assertEqual(len(manifest["streams"]), 1)
        self.assertEqual(expected_value(manifest["records"], datetime(2000, 2, 29)), 290)
        with netCDF4.Dataset(self.rundir / "calendar_noleap_fixtures/datm_in_0.nc") as dataset:
            self.assertEqual(dataset.variables["time"].calendar, "noleap")
            dates = netCDF4.num2date(dataset.variables["time"][:], dataset.variables["time"].units,
                                    calendar="noleap")
            self.assertFalse(any(date.month == 2 and date.day == 29 for date in dates))
            self.assertEqual(dataset.variables["temperature"].cell_methods, "time: point")
            self.assertTrue((dataset.variables["humidity"][:] == 0.005).all())
        configuration = namelist.parse(in_file=str(self.rundir / "datm_in"))
        self.assertEqual(configuration.get_value("tintalgo"), ["'linear'"])
        self.assertIn("2000 2000 2000", configuration.get_value("streams")[0])

    def test_subdaily_staging(self):
        _, specification = selected_case("allactive-calendar_noleap_subdaily")
        manifest = stage_fixtures(self.case, specification, 6)
        self.assertAlmostEqual(expected_value(manifest["records"], datetime(2000, 2, 29)), 291.2)

    def test_ocean_solar_export_checks(self):
        import netCDF4
        import numpy as np
        from noleap_input_cases import validate_solar_snapshot

        path = self.rundir / "solar.cpl.hi.fixture.nc"
        with netCDF4.Dataset(path, "w") as dataset:
            dataset.createDimension("time", 1)
            dataset.createDimension("ncol", 3)
            for name, values in (("xaoo_Faox_swdn", [100, 0, np.nan]),
                                 ("xaoo_Faox_swup", [-10, 0, np.nan]),
                                 ("xaoo_So_anidr", [0.1, 1, np.nan]), ("domo_mask", [1, 1, 0])):
                dataset.createVariable(name, "f8", ("time", "ncol"))[:] = [values]
        with netCDF4.Dataset(path) as dataset:
            validate_solar_snapshot(dataset, 0)
        for field, value in (("xaoo_Faox_swdn", np.nan), ("xaoo_Faox_swdn", 1.e20),
                             ("xaoo_Faox_swup", 10), ("xaoo_So_anidr", 1.1),
                             ("xaoo_Faox_swdn", netCDF4.default_fillvals["f8"])):
            with self.subTest(field=field, value=value):
                with netCDF4.Dataset(path, "a") as dataset:
                    original = dataset.variables[field][0, 0]
                    dataset.variables[field][0, 0] = value
                with netCDF4.Dataset(path) as dataset:
                    with self.assertRaises(AssertionError):
                        validate_solar_snapshot(dataset, 0)
                with netCDF4.Dataset(path, "a") as dataset:
                    dataset.variables[field][0, 0] = original

    def test_solar_stream_keeps_coszen_without_changing_temperature(self):
        import netCDF4
        from CIME import namelist

        with netCDF4.Dataset(self.rundir / "original.nc", "a") as dataset:
            dataset.createVariable("shortwave", "f8", ("time", "ncol"))[:] = 200
        solar = self.rundir / "datm.streams.txt.solar"
        solar.write_text((self.rundir / "datm.streams.txt.test").read_text().replace(
            "temperature tbot\nhumidity qbot", "shortwave swdn"))
        path = self.rundir / "datm_in"
        configuration = namelist.parse(in_file=str(path))
        configuration.set_variable_value("shr_strdata_nml", "streams", [
            "'datm.streams.txt.test 1999 1999 2000'", "'datm.streams.txt.solar 1999 1999 2000'"], var_size=2)
        configuration.write(str(path))
        _, specification = selected_case("allactive-calendar_noleap_daily")
        manifest = stage_fixtures(self.case, specification, 6)
        self.assertEqual(manifest["solar_streams"], [str(solar)])
        configuration = namelist.parse(in_file=str(path))
        self.assertEqual(configuration.get_value("tintalgo"), ["'linear'", "'coszen'"])

    def test_reverse_fixture_retains_leapday(self):
        import netCDF4

        _, specification = selected_case("allactive-calendar_gregorian_to_noleap_daily")
        manifest = stage_fixtures(self.case, specification, 6)
        self.assertEqual(expected_value(manifest["records"], datetime(2000, 2, 28, 12)), 304)
        with netCDF4.Dataset(manifest["files"][0]) as dataset:
            time = dataset.variables["time"]
            dates = netCDF4.num2date(time[:], time.units, calendar=time.calendar)
            self.assertEqual(sum(date.month == 2 and date.day == 29 for date in dates), 1)
            self.assertEqual(time.calendar, "gregorian")

    def test_reverse_cycle_and_export_calendar(self):
        import netCDF4
        import numpy as np

        self.case = SimpleNamespace(get_value={"RUNDIR": str(self.rundir), "RUN_STARTDATE": "1999-01-01"}.get)
        _, specification = selected_case("allactive-calendar_gregorian_to_noleap_cycling")
        manifest = stage_fixtures(self.case, specification, 1460)
        self.assertEqual(manifest["source_year"], 2000)
        self.assertEqual(manifest["records"][0][0], "2000-01-01T00:00:00")
        dates = netCDF4.num2date(np.arange(1460), "days since 1999-01-01", calendar="noleap")
        self.assertFalse(any(date.month == 2 and date.day == 29 for date in dates))
        self.assertEqual((dates[-1].year, dates[-1].month, dates[-1].day), (2002, 12, 31))
        path = self.rundir / "case.cpl.hi.reverse.nc"
        with netCDF4.Dataset(path, "w") as dataset:
            dataset.createDimension("time", 1)
            dataset.createDimension("ncol", 2)
            time = dataset.createVariable("time", "f8", ("time",))
            time.units = "days since 2001-02-28"
            time.calendar = "noleap"
            time[:] = [0.5]
            dataset.createVariable("a2x_Sa_tbot", "f8", ("time", "ncol"))[:] = 304
        self.assertEqual(validate_exports(self.rundir, manifest), (1, 0))
        with netCDF4.Dataset(path, "a") as dataset:
            dataset.variables["time"].calendar = "gregorian"
        with self.assertRaisesRegex(AssertionError, "must declare a NO_LEAP"):
            validate_exports(self.rundir, manifest)

    def test_reverse_negative_fixture_staging(self):
        import netCDF4
        from CIME import namelist

        (self.rundir / "drv_in").write_text("&seq_infodata_inparm\n restart_pfile='rpointer.cpl'\n/\n")
        _, specification = selected_case("allactive-calendar_gregorian_to_noleap_override_conflict")
        manifest = stage_fixtures(self.case, specification, 6)
        generated = namelist.parse(in_file=str(self.rundir / "drv_in"))
        self.assertIn("file:" + manifest["files"][0], generated.get_value("input_calendar_override_selectors")[0])
        self.assertEqual(generated.get_value("input_calendar_override_calendars"), ["'noleap'"])
        self.setUp()
        _, specification = selected_case("allactive-calendar_gregorian_to_noleap_mixed_calendar")
        manifest = stage_fixtures(self.case, specification, 6)
        self.assertEqual(len(manifest["files"]), 2)
        for filename, calendar in zip(manifest["files"], ("gregorian", "noleap")):
            with netCDF4.Dataset(filename) as dataset:
                self.assertEqual(dataset.variables["time"].calendar, calendar)

    def test_reverse_restart_marker_mutations(self):
        import netCDF4
        from noleap_input_cases import mutate_restart_marker

        (self.rundir / "drv_in").write_text("&seq_infodata_inparm\n restart_pfile='rpointer.cpl'\n/\n")
        (self.rundir / "rpointer.cpl").write_text("restart.nc\n")
        path = self.rundir / "restart.nc"
        for action in ("corrupt", "missing"):
            with netCDF4.Dataset(path, "w") as dataset:
                dataset.createVariable("seq_infodata_gregorian_input", "i4")[...] = 1
            mutate_restart_marker(self.case, action)
            with netCDF4.Dataset(path) as dataset:
                if action == "missing":
                    self.assertNotIn("seq_infodata_gregorian_input", dataset.variables)
                else:
                    self.assertEqual(dataset.variables["seq_infodata_gregorian_input"][...], 99)

    def test_packed_forcing_keeps_decoded_physical_values(self):
        import netCDF4

        with netCDF4.Dataset(self.rundir / "original.nc", "a") as dataset:
            dataset.renameVariable("humidity", "unused_humidity")
            variable = dataset.createVariable("humidity", "i2", ("time", "ncol"))
            variable.scale_factor = 1e-4
            variable[:] = 0.005
        _, specification = selected_case("allactive-calendar_noleap_daily")
        stage_fixtures(self.case, specification, 6)
        with netCDF4.Dataset(self.rundir / "calendar_noleap_fixtures/datm_in_0.nc") as dataset:
            self.assertTrue((abs(dataset.variables["humidity"][:] - 0.005) < 1e-12).all())
            self.assertNotIn("scale_factor", dataset.variables["humidity"].ncattrs())
            self.assertEqual(dataset.variables["temperature"].units, "K")

    def test_multiple_stream_arrays_are_preserved(self):
        from CIME import namelist

        first = self.rundir / "datm.streams.txt.test"
        second = self.rundir / "datm.streams.txt.second"
        second.write_text(first.read_text().replace("temperature tbot\n", ""))
        path = self.rundir / "datm_in"
        configuration = namelist.parse(in_file=str(path))
        configuration.set_variable_value("shr_strdata_nml", "streams", [
            "'datm.streams.txt.test 1999 1999 2000'",
            "'datm.streams.txt.second 1999 1999 2000'"], var_size=2)
        configuration.write(str(path))
        _, specification = selected_case("allactive-calendar_noleap_daily")
        self.assertEqual(len(stage_fixtures(self.case, specification, 6)["streams"]), 2)
        configuration = namelist.parse(in_file=str(path))
        for name in ("streams", "tintalgo", "taxmode", "offset", "dtlimit"):
            self.assertEqual(len(configuration.get_value(name)), 2, name)

    def test_missing_calendar_fixture(self):
        import netCDF4

        _, specification = selected_case("allactive-calendar_noleap_missing_calendar")
        stage_fixtures(self.case, specification, 6)
        with netCDF4.Dataset(self.rundir / "calendar_noleap_fixtures/datm_in_0.nc") as dataset:
            self.assertNotIn("calendar", dataset.variables["time"].ncattrs())

    def test_bounded_and_declared_mean_negative_fixtures(self):
        import netCDF4

        for name in ("mean", "ambiguous_bounds", "empty_calendar"):
            with self.subTest(fixture=name):
                self.setUp()
                _, specification = selected_case("allactive-calendar_noleap_" + name)
                stage_fixtures(self.case, specification, 6)
                with netCDF4.Dataset(self.rundir / "calendar_noleap_fixtures/datm_in_0.nc") as dataset:
                    if name == "mean":
                        self.assertEqual(dataset.variables["temperature"].cell_methods, "time: mean")
                    elif name == "ambiguous_bounds":
                        self.assertEqual(dataset.variables["time"].bounds, "time_bnds")
                        self.assertEqual(dataset.variables["time_bnds"].shape[1], 2)
                        self.assertNotIn("cell_methods", dataset.variables["temperature"].ncattrs())
                    else:
                        self.assertEqual(dataset.variables["time"].calendar, "")

    def test_export_validation_detects_wrong_values(self):
        import netCDF4

        _, specification = selected_case("allactive-calendar_noleap_daily")
        manifest = stage_fixtures(self.case, specification, 6)
        path = self.rundir / "case.cpl.hi.2000-02-29-00000.nc"
        with netCDF4.Dataset(path, "w") as dataset:
            dataset.createDimension("time", 1)
            dataset.createDimension("ncol", 2)
            time = dataset.createVariable("time", "f8", ("time",))
            time.units = "days since 2000-01-01"
            time.calendar = "gregorian"
            time[:] = [59]
            field = dataset.createVariable("a2x_Sa_tbot", "f8", ("time", "ncol"))
            field[:] = 290
        self.assertEqual(validate_exports(self.rundir, manifest), (1, 1))
        with self.assertRaisesRegex(AssertionError, "No instantaneous"):
            validate_exports(self.rundir, manifest, export_snapshot(self.rundir))
        with netCDF4.Dataset(path, "a") as dataset:
            dataset.variables["a2x_Sa_tbot"][0, 1] = 288
        with self.assertRaisesRegex(AssertionError, "Wrong forcing export"):
            validate_exports(self.rundir, manifest)

    def test_core2_export_excludes_downstream_polar_corrections(self):
        import netCDF4

        manifest = {"records": [("2000-02-28T00:00:00", 280), ("2000-03-01T00:00:00", 284)],
                    "enabled": True, "calendar": "noleap", "datamode": "CORE2_NYF"}
        path = self.rundir / "case.cpl.hi.2000-02-29-00000.nc"
        with netCDF4.Dataset(path, "w") as dataset:
            dataset.createDimension("time", 1)
            dataset.createDimension("ncol", 3)
            time = dataset.createVariable("time", "f8", ("time",))
            time.units = "days since 2000-02-27 00:00:00"
            time.calendar = "gregorian"
            time[:] = [2]
            dataset.createVariable("doma_lat", "f8", ("time", "ncol"))[:] = [[-80, 0, 80]]
            dataset.createVariable("a2x_Sa_tbot", "f8", ("time", "ncol"))[:] = [[285, 282, 279]]
        self.assertEqual(validate_exports(self.rundir, manifest), (1, 1))
        with netCDF4.Dataset(path, "a") as dataset:
            dataset.variables["a2x_Sa_tbot"][0, 1] = 281
        with self.assertRaisesRegex(AssertionError, "Wrong forcing export"):
            validate_exports(self.rundir, manifest)


if __name__ == "__main__":
    unittest.main()