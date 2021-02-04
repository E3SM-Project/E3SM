from .core_parameter import CoreParameter


class StreamflowParameter(CoreParameter):
    def __init__(self):
        super(StreamflowParameter, self).__init__()
        self.gauges_path = None
        self.main_title_seasonality_map = "Seasonality Map"
        self.main_title_annual_map = "Mean Annual Streamflow Map"
        self.main_title_annual_scatter = "Mean Annual Streamflow Scatter Plot"
        self.max_num_gauges = None
        self.output_file_seasonality_map = "seasonality_map"
        self.output_file_annual_map = "annual_map"
        self.output_file_annual_scatter = "annual_scatter"
        self.print_statements = False
        self.ref_timeseries_input = True
        self.test_timeseries_input = True
        self.granulate.remove("seasons")

    # TODO: Other Parameter classes have the same check_values method. Move to CoreParameter.
    # FIXME: start_yr and end_yr attributes never seem to be instantiated
    def check_values(self):
        test_ref_start_yr_both_set = hasattr(
            self, "test_start_yr"
        ) and hasattr(self, "ref_start_yr")
        if hasattr(self, "start_yr"):
            # Use `start_yr` as a default value for other parameters.
            if not hasattr(self, "test_start_yr"):
                # FIXME: error: Cannot determine type of 'start_yr'
                self.test_start_yr = self.start_yr  # type: ignore
            if not hasattr(self, "ref_start_yr"):
                # FIXME: error: Cannot determine type of 'start_yr'
                self.ref_start_yr = self.start_yr  # type: ignore
        elif (
            test_ref_start_yr_both_set
            and self.test_start_yr == self.ref_start_yr
        ):
            # Derive the value of self.start_yr
            self.start_yr = self.test_start_yr

        test_ref_end_yr_both_set = hasattr(self, "test_end_yr") and hasattr(
            self, "ref_end_yr"
        )
        if hasattr(self, "end_yr"):
            # Use `end_yr` as a default value for other parameters.
            if not hasattr(self, "test_end_yr"):
                # FIXME: error: Cannot determine type of 'end_yr'
                self.test_end_yr = self.end_yr  # type: ignore
            if not hasattr(self, "ref_end_yr"):
                # FIXME: error: Cannot determine type of 'end_yr'
                self.ref_end_yr = self.end_yr  # type: ignore
        elif test_ref_end_yr_both_set and self.test_end_yr == self.ref_end_yr:
            # Derive the value of self.end_yr
            self.end_yr = self.test_end_yr

        if hasattr(self, "start_yr"):
            # We need to re-evaluate this variable, since these attributes could have been set.
            test_ref_end_yr_both_set = hasattr(
                self, "test_end_yr"
            ) and hasattr(self, "ref_end_yr")
            if not (hasattr(self, "end_yr") or test_ref_end_yr_both_set):
                msg = "To use 'start_yr' you need to also define 'end_yr' or both 'test_end_yr' and 'ref_end_yr'."
                raise RuntimeError(msg)

        if hasattr(self, "end_yr"):
            # We need to re-evaluate this variable, since these attributes could have been set.
            test_ref_start_yr_both_set = hasattr(
                self, "test_start_yr"
            ) and hasattr(self, "ref_start_yr")
            if not (hasattr(self, "start_yr") or test_ref_start_yr_both_set):
                msg = "To use 'end_yr' you need to also define 'start_yr' or both 'test_start_yr' and 'ref_start_yr'."
                raise RuntimeError(msg)
