import unittest
import os
import re
import subprocess
from PIL import Image, ImageChops


# Run these tetsts on Cori by doing the following:
# cd tests/system
# module load python/2.7-anaconda-4.4
# source activate e3sm_diags_env_dev
# If code in acme_diags has been changed:
# pip install /global/homes/f/<username>/e3sm_diags/
# python test_diags.py


# Set to True to place the results directory on Cori's web server
# Set to False to place the results directory in tests/system
CORI_WEB = False


def get_results_dir(output_list):
    """Given output from acme_diags_driver, extract the path to results_dir."""
    for line in output_list:
        match = re.search('Viewer HTML generated at (.*)viewer.*.html', line)
        if match:
            results_dir = match.group(1)
            return results_dir
    message = 'No viewer directory listed in output: {}'.format(output_list)
    raise RuntimeError(message)


def move_to_web(machine_path_re_str, html_prefix_format_str, results_dir):
    command = 'git rev-parse --show-toplevel'
    top_level = subprocess.check_output(command.split()).decode('utf-8').splitlines()[0]
    match = re.search(machine_path_re_str, top_level)
    if match:
        username = match.group(1)
    else:
        message = 'Username could not be extracted from top_level={}'.format(top_level)
        raise RuntimeError(message)
    html_prefix = html_prefix_format_str.format(username)
    print('html_prefix={}'.format(html_prefix))
    new_results_dir = '{}/{}'.format(html_prefix, results_dir)
    print('new_results_dir={}'.format(new_results_dir))
    if os.path.exists(new_results_dir):
        command = 'rm -r {}'.format(new_results_dir)
        subprocess.check_output(command.split())
    command = 'mv {} {}'.format(results_dir, new_results_dir)
    subprocess.check_output(command.split())
    command = 'chmod -R 755 {}'.format(new_results_dir)
    subprocess.check_output(command.split())
    return new_results_dir


def count_images(directory):
    images = []
    for _,_,filenames in os.walk(directory):
        for file in filenames:
            if file.endswith('.png'):
                images.append(file)
    return len(images)


def compare_images(test, mismatched_images, image_name, path_to_actual_png, path_to_expected_png):
    # https://stackoverflow.com/questions/35176639/compare-images-python-pil

    actual_png = Image.open(path_to_actual_png).convert('RGB')
    expected_png = Image.open(path_to_expected_png).convert('RGB')
    diff = ImageChops.difference(actual_png, expected_png)

    bbox = diff.getbbox()
    if not bbox:
        # If `diff.getbbox()` is None, then the images are in theory equal
        test.assertIsNone(diff.getbbox())
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
        fraction = num_nonzero_pixels / num_pixels
        print('num_nonzero_pixels/num_pixels fraction={}'.format(fraction))
        # Fraction of mismatched pixels should be less than 0.02%
        if fraction >= 0.0002:
            mismatched_images.append(image_name)


class TestAllSets(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        command = 'python all_sets.py -d all_sets.cfg'
        output_list = subprocess.check_output(command.split()).decode('utf-8').splitlines()
        TestAllSets.results_dir = get_results_dir(output_list)
        print('TestAllSets.results_dir={}'.format(TestAllSets.results_dir))
        if CORI_WEB:
            TestAllSets.results_dir = move_to_web(
                '/global/u1/f/(.*)/e3sm_diags', '/global/cfs/cdirs/acme/www/{}',
                TestAllSets.results_dir)

    def check_html_image(self, html_path, png_path, full_png_path):
        img_src = None
        option_value = None
        href = None
        with open(html_path, 'r') as html:
            for line in html:
                # If `img_src` is not defined yet:
                if not img_src:
                    re_str = '<img src="../../../../{}">'.format(png_path)
                    img_src = re.search(re_str, line)
                # If `option_value` is not defined yet:
                if not option_value:
                    re_str = '<option value="../../../../{}">'.format(png_path)
                    option_value = re.search(re_str, line)
                # If `href` is not defined yet:
                if not href:
                    re_str = 'href="../../../../{}">'.format(png_path)
                    href = re.search(re_str, line)
        self.assertIsNotNone(img_src)
        self.assertIsNotNone(option_value)
        self.assertIsNotNone(href)

        image_name = os.path.split(png_path)[-1]
        path_to_actual_png = full_png_path
        path_to_expected_png = '../unit_test_images/{}'.format(png_path)

        mismatched_images = []
        compare_images(self, mismatched_images, image_name, path_to_actual_png, path_to_expected_png)
        self.assertEqual(mismatched_images, [])

    def check_plots_generic(self, set_name, case_id, ref_name, variables, region, plev=None):
        case_id_lower = case_id.lower()
        ref_name_lower = ref_name.lower()
        region_lower = region.lower()
        seasons = ['ANN']
        for variable in variables:
            variable_lower = variable.lower()
            for season in seasons:
                season_lower = season.lower()
                if plev:
                    # 200.9 would just show up as 200 in the file paths.
                    plev_str = '%.0f' % plev
                    plev_png_str = '{}-'.format(plev_str)
                    plev_html_str = '{}mb-'.format(plev_str)
                else:
                    plev_png_str = ''
                    plev_html_str = ''
                png_path = '{}/{}/{}-{}-{}{}-{}.png'.format(
                    set_name, case_id, ref_name, variable, plev_png_str, season, region)
                full_png_path = '{}{}'.format(TestAllSets.results_dir, png_path)
                self.assertTrue(os.path.exists(full_png_path))
                html_path = '{}viewer/{}/{}/{}-{}{}-{}/{}.html'.format(
                    TestAllSets.results_dir, set_name, case_id_lower, variable_lower,
                    plev_html_str, region_lower, ref_name_lower, season_lower)
                self.check_html_image(html_path, png_path, full_png_path)

    def check_plots_2d(self, set_name):
        self.check_plots_generic(
            set_name=set_name,
            case_id='ERA-Interim',
            ref_name='ERA-Interim',
            variables=['T'],
            region='global'
        )

    def check_plots_plevs(self, set_name, region, plevs):
        for plev in plevs:
            self.check_plots_generic(
                set_name=set_name,
                case_id='ERA-Interim',
                ref_name='ERA-Interim',
                variables=['T'],
                region=region,
                plev=plev
            )

    def check_enso_map_plots(self, case_id):
        case_id_lower = case_id.lower()
        nino_region_lower = 'NINO34'.lower()
        set_name = 'enso_diags'
        variables = ['TREFHT']
        for variable in variables:
            variable_lower = variable.lower()
            png_path = '{}/{}/regression-coefficient-{}-over-{}.png'.format(
                set_name, case_id, variable_lower, nino_region_lower)
            full_png_path = '{}{}'.format(TestAllSets.results_dir, png_path)
            self.assertTrue(os.path.exists(full_png_path))
            html_path = '{}viewer/{}/map/{}/plot.html'.format(
                TestAllSets.results_dir, set_name, case_id_lower)
            self.check_html_image(html_path, png_path, full_png_path)

    def check_enso_scatter_plots(self, case_id):
        case_id_lower = case_id.lower()
        set_name = 'enso_diags'
        variables = ['TREFHT']
        for variable in variables:
            region = 'NINO3'
            png_path = '{}/{}/feedback-{}-{}-TS-NINO3.png'.format(
                set_name, case_id, variable, region)
            full_png_path = '{}{}'.format(TestAllSets.results_dir, png_path)
            self.assertTrue(os.path.exists(full_png_path))
            html_path = '{}viewer/{}/scatter/{}/plot.html'.format(
                TestAllSets.results_dir, set_name, case_id_lower)
            self.check_html_image(html_path, png_path, full_png_path)

    def check_streamflow_plots(self):
        case_id = 'RIVER_DISCHARGE_OVER_LAND_LIQ_GSIM'
        case_id_lower = case_id.lower()
        set_name = 'streamflow'
        variables = ['RIVER_DISCHARGE_OVER_LAND_LIQ']
        for variable in variables:
            for plot_type in ['seasonality_map', 'annual_map', 'annual_scatter']:
                png_path = '{}/{}/{}.png'.format(
                    set_name, case_id, plot_type)
                expected = 'streamflow/RIVER_DISCHARGE_OVER_LAND_LIQ_GSIM/{}.png'.format(plot_type)
                self.assertEqual(png_path, expected)
                full_png_path = '{}{}'.format(TestAllSets.results_dir, png_path)
                self.assertTrue(os.path.exists(full_png_path))
                if plot_type == 'seasonality_map':
                    plot_label = 'seasonality-map'
                elif plot_type == 'annual_map':
                    plot_label = 'mean-annual-streamflow-map'
                elif plot_type == 'annual_scatter':
                    plot_label = 'mean-annual-streamflow-scatter-plot'
                else:
                    raise Exception('Invalid plot_type={}'.format(plot_type))
                html_path = 'viewer/{}/{}/{}-{}/plot.html'.format(
                    set_name, plot_label, case_id_lower, plot_type)
                expected = 'viewer/streamflow/{}/river_discharge_over_land_liq_gsim-{}/plot.html'.format(
                    plot_label, plot_type)
                self.assertEqual(html_path, expected)
                full_html_path = '{}{}'.format(TestAllSets.results_dir, html_path)
                self.check_html_image(full_html_path, png_path, full_png_path)

    # Test results_dir
    def test_results_dir(self):
        self.assertTrue(TestAllSets.results_dir.endswith('all_sets_results_test/'))

    # Test the image count
    def test_image_count(self):
        actual_num_images = count_images('all_sets_results_test')
        expected_num_images = count_images('../unit_test_images')
        self.assertEqual(actual_num_images, expected_num_images)

    # Test sets
    def test_area_mean_time_series(self):
        set_name = 'area_mean_time_series'
        variables = ['TREFHT']
        for variable in variables:
            variable_lower = variable.lower()
            png_path = '{}/{}.png'.format(
                set_name, variable)
            full_png_path = '{}{}'.format(TestAllSets.results_dir, png_path)
            self.assertTrue(os.path.exists(full_png_path))
            html_path = '{}viewer/{}/variable/{}/plot.html'.format(
                TestAllSets.results_dir, set_name, variable_lower)
            self.check_html_image(html_path, png_path, full_png_path)

    def test_cosp_histogram(self):
        self.check_plots_generic(
            set_name='cosp_histogram',
            case_id='MISR-COSP',
            ref_name='MISRCOSP',
            variables=['COSP_HISTOGRAM_MISR'],
            region='global'
        )

    def test_enso_diags_map(self):
        case_id = 'TREFHT-response-map'
        self.check_enso_map_plots(case_id)

    def test_enso_diags_map_start_yrs(self):
        case_id = 'TREFHT-response-map-start-yrs'
        self.check_enso_map_plots(case_id)

    def test_enso_diags_map_test_ref_yrs(self):
        case_id = 'TREFHT-response-map-test-ref-yrs'
        self.check_enso_map_plots(case_id)

    def test_enso_diags_scatter(self):
        case_id = 'TREFHT-response-scatter'
        self.check_enso_scatter_plots(case_id)

    def test_enso_diags_scatter_start_yrs(self):
        case_id = 'TREFHT-response-scatter-start-yrs'
        self.check_enso_scatter_plots(case_id)

    def test_enso_diags_scatter_test_ref_yrs(self):
        case_id = 'TREFHT-response-scatter-test-ref-yrs'
        self.check_enso_scatter_plots(case_id)

    def test_lat_lon(self):
        self.check_plots_plevs('lat_lon', 'global', [850.0])

    def test_lat_lon_regional(self):
        self.check_plots_plevs('lat_lon', 'CONUS_RRM', [850.0])

    def test_meridional_mean_2d(self):
        self.check_plots_2d('meridional_mean_2d')

    def test_polar(self):
        self.check_plots_plevs('polar', 'polar_S', [850.0])

    def test_qbo(self):
        case_id = 'qbo-test'
        case_id_lower = case_id.lower()
        set_name = 'qbo'
        png_path = '{}/{}/qbo_diags.png'.format(set_name, case_id)
        full_png_path = '{}{}'.format(TestAllSets.results_dir, png_path)
        print(full_png_path)
        self.assertTrue(os.path.exists(full_png_path))
        # viewer/qbo/variable/era-interim/plot.html
        html_path = '{}viewer/{}/variable/{}/plot.html'.format(
            TestAllSets.results_dir, set_name, case_id_lower)
        self.check_html_image(html_path, png_path, full_png_path)

    def test_streamflow(self):
        self.check_streamflow_plots()

    def test_zonal_mean_2d(self):
        self.check_plots_2d('zonal_mean_2d')

    def test_zonal_mean_xy(self):
        self.check_plots_plevs('zonal_mean_xy', 'global', [200.0])


if __name__ == '__main__':
    unittest.main()
