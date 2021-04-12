from .core_parameter import CoreParameter


class TimeSeriesParameter(CoreParameter):
    def __init__(self):
        super(TimeSeriesParameter, self).__init__()

    def check_values(self):
        # Do not call super -- we want to overwrite check_values from CoreParameter
        test_ref_start_yr_both_set = hasattr(self, "test_start_yr") and hasattr(
            self, "ref_start_yr"
        )
        if hasattr(self, "start_yr"):
            # Use `start_yr` as a default value for other parameters.
            if not hasattr(self, "test_start_yr"):
                self.test_start_yr: str = self.start_yr
            if not hasattr(self, "ref_start_yr"):
                self.ref_start_yr: str = self.start_yr
        elif test_ref_start_yr_both_set and self.test_start_yr == self.ref_start_yr:
            # Derive the value of self.start_yr
            self.start_yr: str = self.test_start_yr

        test_ref_end_yr_both_set = hasattr(self, "test_end_yr") and hasattr(
            self, "ref_end_yr"
        )
        if hasattr(self, "end_yr"):
            # Use `end_yr` as a default value for other parameters.
            if not hasattr(self, "test_end_yr"):
                self.test_end_yr: str = self.end_yr
            if not hasattr(self, "ref_end_yr"):
                self.ref_end_yr: str = self.end_yr
        elif test_ref_end_yr_both_set and self.test_end_yr == self.ref_end_yr:
            # Derive the value of self.end_yr
            self.end_yr: str = self.test_end_yr

        if hasattr(self, "start_yr"):
            # We need to re-evaluate this variable, since these attributes could have been set.
            test_ref_end_yr_both_set = hasattr(self, "test_end_yr") and hasattr(
                self, "ref_end_yr"
            )
            if not (hasattr(self, "end_yr") or test_ref_end_yr_both_set):
                msg = "To use 'start_yr' you need to also define 'end_yr' or both 'test_end_yr' and 'ref_end_yr'."
                raise RuntimeError(msg)

        if hasattr(self, "end_yr"):
            # We need to re-evaluate this variable, since these attributes could have been set.
            test_ref_start_yr_both_set = hasattr(self, "test_start_yr") and hasattr(
                self, "ref_start_yr"
            )
            if not (hasattr(self, "start_yr") or test_ref_start_yr_both_set):
                msg = "To use 'end_yr' you need to also define 'start_yr' or both 'test_start_yr' and 'ref_start_yr'."
                raise RuntimeError(msg)
