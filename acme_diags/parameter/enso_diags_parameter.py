from .core_parameter import CoreParameter


class EnsoDiagsParameter(CoreParameter):
    def __init__(self):
        super(EnsoDiagsParameter, self).__init__()
        self.granulate.remove("seasons")
        self.nino_region = "NINO34"
        self.plot_type = "map"
        self.print_statements = False
        self.ref_timeseries_input = True
        self.test_timeseries_input = True

    # TODO: Other Parameter classes have the same check_values method. Move to CoreParameter. EnsoDiagsParameter has some variation (e.g. valid_nino_regions), just extend the method
    # FIXME: start_yr and end_yr attributes never seem to be instantiated
    def check_values(self):
        valid_nino_regions = ["NINO3", "NINO34", "NINO4"]
        if self.nino_region not in valid_nino_regions:
            msg = "nino_region={} not in {}".format(
                self.nino_region, valid_nino_regions
            )
            raise RuntimeError(msg)

        valid_plot_types = ["map", "scatter"]
        if self.plot_type not in valid_plot_types:
            msg = "plot_type={} not in {}".format(self.plot_type, valid_plot_types)
            raise RuntimeError(msg)

        # NOTE: Everything below this is duplicated with other Parameter classes
        test_ref_start_yr_both_set = hasattr(self, "test_start_yr") and hasattr(
            self, "ref_start_yr"
        )
        if hasattr(self, "start_yr"):
            # Use `start_yr` as a default value for other parameters.
            if not hasattr(self, "test_start_yr"):
                # FIXME: error: Cannot determine type of 'start_yr'
                self.test_start_yr = self.start_yr  # type: ignore
            if not hasattr(self, "ref_start_yr"):
                # FIXME: error: Cannot determine type of 'start_yr'
                self.ref_start_yr = self.start_yr  # type: ignore
        elif test_ref_start_yr_both_set and self.test_start_yr == self.ref_start_yr:
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
