"""Prescribed no-leap forcing values and exact restart equivalence."""

from math import gcd

from CIME.SystemTests.ers import ERS
from noleap_input_cases import (
    export_snapshot, fresh_logs, log_snapshot, selected_case, stage_fixtures,
    validate_activation, validate_exports,
)


class NLF(ERS):
    def __init__(self, case, **kwargs):
        super().__init__(case, **kwargs)
        self._fixture_days = case.get_value("STOP_N")
        self._fixture_name, self._specification = selected_case(case.get_value("CASE"))
        if "error" in self._specification:
            raise ValueError("Negative fixture requires NLN, not NLF")

    def run_indv(self, *args, **kwargs):
        frequency = self._case.get_value("ATM_NCPL")
        if self._case.get_value("NCPL_BASE_PERIOD") != "day" or frequency <= 0:
            raise ValueError("NLF requires positive day-based atmosphere coupling frequency")
        self._case.set_value("HIST_OPTION", "nhours")
        self._case.set_value("HIST_N", 24 // gcd(24, frequency))
        self._case.create_namelists()
        manifest = stage_fixtures(self._case, self._specification, self._fixture_days)
        self._skip_pnl = True
        rundir = self._case.get_value("RUNDIR")
        before = log_snapshot(rundir)
        previous_exports = export_snapshot(rundir)
        super().run_indv(*args, **kwargs)
        samples, leap_samples = validate_exports(rundir, manifest, previous_exports)
        logs = fresh_logs(rundir, before)
        validate_activation(manifest, logs)
        reverse = self._specification.get("direction") == "gregorian_to_noleap"
        if reverse and leap_samples != 0:
            raise AssertionError("NO_LEAP history contains February 29")
        if not reverse and "calendar_feb" in self._case.get_value("CASE") and leap_samples == 0:
            raise AssertionError("No February 29 forcing samples were checked")
        if samples < 2:
            raise AssertionError("Insufficient forcing-export samples")