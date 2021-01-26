import unittest
import os
from PIL import Image, ImageChops
from acme_diags.run import runner
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
from acme_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from acme_diags.parameter.qbo_parameter import QboParameter

# Due to the large amount of data required to run, this test will be run manually on Anvil
# (rather than as part of the CI tests).

class TestCompleteRun(unittest.TestCase):
    def test_complete_run(self):
        # Anvil
        # Run `source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified.sh`
        test_data_prefix = '/lcrc/group/e3sm/public_html/e3sm_diags_test_data'
        ref_data_prefix = '/lcrc/group/e3sm/public_html/diagnostics/observations/Atm'
        html_prefix = '.'

        param = CoreParameter()

        param.reference_data_path = os.path.join(ref_data_prefix, 'climatology')
        param.test_data_path = os.path.join(test_data_prefix, 'climatology/')
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN","JJA"]  # Default setting: seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

        param.results_dir = os.path.join(html_prefix, 'tutorial_2020_all_sets')
        param.multiprocessing = True
        param.num_workers = 30

        # Additional parameters:
        #param.short_test_name = 'beta0.FC5COSP.ne30'
        #param.run_type = 'model_vs_model'
        #param.diff_title = 'Difference'
        #param.output_format = ['png']
        #param.output_format_subplot = ['pdf']
        #param.save_netcdf = True


        # Set specific parameters for new sets
        enso_param = EnsoDiagsParameter()
        enso_param.reference_data_path = os.path.join(ref_data_prefix, 'time-series/')
        enso_param.test_data_path = os.path.join(test_data_prefix, 'time-series/E3SM_v1/')
        enso_param.test_name = 'e3sm_v1'
        enso_param.start_yr = '1990'
        enso_param.end_yr = '1999'

        qbo_param = QboParameter()
        qbo_param.reference_data_path = os.path.join(ref_data_prefix, 'time-series/')
        qbo_param.test_data_path = os.path.join(test_data_prefix, 'time-series/E3SM_v1/')
        qbo_param.test_name = 'e3sm_v1'
        qbo_param.start_yr = '1990'
        qbo_param.end_yr = '1999'

        ts_param = AreaMeanTimeSeriesParameter()
        ts_param.reference_data_path = os.path.join(ref_data_prefix, 'time-series/')
        ts_param.test_data_path = os.path.join(test_data_prefix, 'time-series/E3SM_v1/')
        ts_param.test_name = 'e3sm_v1'
        ts_param.start_yr = '1990'
        ts_param.end_yr = '1999'

        runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d', 'enso_diags', 'qbo','area_mean_time_series']
        runner.run_diags([param, enso_param, qbo_param, ts_param])

        actual_images_dir = param.results_dir
        # The expected_images_file lists all 475 images we expect to compare.
        # It was generated with the following steps:
        # cd /lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/tutorial_2020_all_sets
        # find . -type f -name '*.png' > ../expected_images_complete_run.txt
        expected_images_file = '/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/expected_images_complete_run.txt'
        expected_images_dir = '/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/tutorial_2020_all_sets'

        mismatched_images = []

        with open(expected_images_file) as f:
            for line in f:
                image = line.strip('./').strip('\n')
                path_to_actual_png = os.path.join(actual_images_dir, image)
                path_to_expected_png = os.path.join(expected_images_dir, image)

                actual_png = Image.open(path_to_actual_png).convert('RGB')
                expected_png = Image.open(path_to_expected_png).convert('RGB')
                diff = ImageChops.difference(actual_png, expected_png)

                bbox = diff.getbbox()
                if not bbox:
                    # If `diff.getbbox()` is None, then the images are in theory equal
                    self.assertIsNone(diff.getbbox())
                else:
                    # Sometimes, a few pixels will differ, but the two images appear identical.
                    # https://codereview.stackexchange.com/questions/55902/fastest-way-to-count-non-zero-pixels-using-python-and-pillow
                    nonzero_pixels = diff.crop(bbox).point(lambda x: 255 if x else 0).convert("L").point(bool).getdata()
                    num_nonzero_pixels = sum(nonzero_pixels)
                    print('\npath_to_actual_png={}'.format(path_to_actual_png))
                    print('path_to_expected_png={}'.format(path_to_expected_png))
                    print('diff has {} nonzero pixels.'.format(num_nonzero_pixels))
                    width, height = expected_png.size
                    num_pixels = width * height
                    print('total number of pixels={}'.format(num_pixels))
                    fraction = num_nonzero_pixels/num_pixels
                    print('num_nonzero_pixels/num_pixels fraction={}'.format(fraction))
                    # Fraction of mismatched pixels should be less than 0.02%
                    if fraction >= 0.0002:
                        mismatched_images.append(image)

        self.assertEqual(mismatched_images, [])

if __name__ == '__main__':
    unittest.main()
