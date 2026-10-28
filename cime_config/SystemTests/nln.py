"""Require a specific no-leap input rejection, never an unrelated run failure."""

from CIME.SystemTests.ers import ERS
from CIME.utils import CIMEError
from noleap_input_cases import (
    fresh_logs, log_snapshot, require_expected_failure, selected_case,
    set_user_policy, stage_fixtures,
    mutate_restart_marker,
)


class NLN(ERS):
    def __init__(self, case, **kwargs):
        super().__init__(case, **kwargs)
        self._fixture_days = case.get_value("STOP_N")
        self._fixture_name, self._specification = selected_case(case.get_value("CASE"))

    def run_indv(self, *args, **kwargs):
        self._case.create_namelists()
        stage_fixtures(self._case, self._specification, self._fixture_days)
        if "restart_marker" in self._specification and self._case.get_value("CONTINUE_RUN"):
            mutate_restart_marker(self._case, self._specification["restart_marker"])
        self._skip_pnl = True
        return super().run_indv(*args, **kwargs)

    def run_phase(self):
        name, specification = self._fixture_name, self._specification
        if "error" not in specification:
            raise ValueError("Positive fixture requires NLF, not NLN: " + name)
        if "restart_enabled" in specification:
            def policy(enabled):
                if specification.get("direction") == "gregorian_to_noleap":
                    set_user_policy(self._case, enabled, "gregorian_to_noleap")
                else:
                    set_user_policy(self._case, enabled)
            policy(specification["enabled"])
            self._ers_first_phase()
            policy(specification["restart_enabled"])
            run = self._ers_second_phase
        elif "restart_override" in specification or "restart_marker" in specification:
            self._specification = dict(specification)
            self._specification.pop("override", None)
            self._ers_first_phase()
            if "restart_override" in specification:
                self._specification["override"] = "gregorian"
            run = self._ers_second_phase
        else:
            run = self.run_indv
        rundir = self._case.get_value("RUNDIR")
        before = log_snapshot(rundir)
        failure = None
        try:
            run()
        except CIMEError as error:
            failure = error
        require_expected_failure(failure, specification["error"], fresh_logs(rundir, before))