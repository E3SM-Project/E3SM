import unittest
import collections
from acme_diags.run import Run
from acme_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from acme_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter


class TestRun(unittest.TestCase):
    """
    Tests running e3sm_diags via the API.
    """
    def setUp(self):
        self.runner = Run()
        self.core_param = CoreParameter()
        
        # These parameters don't need to map to real files since we're
        # just testing if the parameters are created correctly.
        self.core_param.reference_data_path = '/this/can/be/whatever/for/'
        self.core_param.reference_data_path += 'testing/obs_for_acme_diags/'
        self.core_param.test_data_path = '/this/can//also/be/whatever/for/testing/'
        self.core_param.test_data_path += 'test_model_data_for_acme_diags/climatology/'
        self.core_param.test_name = 'SomeTimeStamp.beta?.SomeModelThing.ne30_ne30.somemachine'

        self.core_param.results_dir = 'results'

    def test_lat_lon_ann(self):
        self.core_param.seasons = ['ANN']
        self.runner.sets_to_run = ['lat_lon']
        parameters = self.runner.get_final_parameters([self.core_param])

        for param in parameters:
            bad_seasons = ['DJF', 'MAM', 'JJA', 'SON']
            for season in bad_seasons:
                if season in param.seasons:
                    msg = '{} shouldn\'t be a season in this parameter.'
                    self.fail(msg.format(season))

    def test_all_sets_and_all_seasons(self):
        ts_param = AreaMeanTimeSeriesParameter()
        ts_param.start_yr = '2000'
        ts_param.end_yr = '2004'

        enso_param = EnsoDiagsParameter()
        enso_param.start_yr = '2000'
        enso_param.end_yr = '2004'
        
        parameters = self.runner.get_final_parameters([self.core_param, ts_param, enso_param])
        # Counts the number of each set and each seasons to run the diags on.
        set_counter, season_counter = collections.Counter(), collections.Counter()
        for param in parameters:
            for set_name in param.sets:
                set_counter[set_name] += 1
            for season in param.seasons:
                season_counter[season] += 1
        
        for set_name in set_counter:
            count = set_counter[set_name]
            if count <= 0:
                msg = 'Count for {} is invalid: {}'
                self.fail(msg.format(set_name, count))
        
        all_season_counts = list(season_counter.values())
        # enso_diags only runs ANN, no seasons
        # So, reduce the ANN count by the number of times enso_diags appears
        all_season_counts[0] -= set_counter['enso_diags']
        if not all(all_season_counts[0] == count for count in all_season_counts):
            self.fail('Counts for the seasons don\'t match: {}'.format(all_season_counts))

    def test_zonal_mean_2d(self):
        # Running zonal_mean_2d with the core param only.
        self.runner.sets_to_run = ['zonal_mean_2d']
        core_only_results = self.runner.get_final_parameters([self.core_param])

        # Running zonal_mean_2d with a set-specific param.
        # We pass in both the core and this parameter.
        zonal_mean_2d_param = ZonalMean2dParameter()
        zonal_mean_2d_param.plevs = [10.0, 20.0, 30.0]
        both_results = self.runner.get_final_parameters([self.core_param, zonal_mean_2d_param])
        
        # Check that the plevs value is actually changed.
        for param in both_results:
            if param.plevs != zonal_mean_2d_param.plevs:
                msg = 'plevs are {} when they should be {}.'
                self.fail(msg.format(param.plevs, zonal_mean_2d_param.plevs))

        # Running zonal_mean_2d without a CoreParameter.
        # When doing so, we still must have valid values for
        # certain parameters.
        another_zonal_mean_2d_param = ZonalMean2dParameter()
        another_zonal_mean_2d_param.reference_data_path = '/something'
        another_zonal_mean_2d_param.test_data_path = '/something/else'
        another_zonal_mean_2d_param.results_dir = '/something/else/too'
        zm_2d_only_results = self.runner.get_final_parameters([another_zonal_mean_2d_param])

        # Now checking that all three configs have the same output.
        len1, len2, len3 = len(core_only_results), len(both_results), len(zm_2d_only_results)
        if len1 == 0 or not(len1 == len2 == len3):
            msg = 'The lengths are either 0 or not equal: {} {} {}'
            self.fail(msg.format(len1, len2, len3))

if __name__ == '__main__':
    unittest.main()
