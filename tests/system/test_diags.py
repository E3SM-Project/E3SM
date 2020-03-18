import unittest
import os
import re
import subprocess


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


class TestAllSets(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        command = 'python all_sets.py -d all_sets.cfg'
        output_list = subprocess.check_output(command.split()).decode('utf-8').splitlines()
        TestAllSets.results_dir = get_results_dir(output_list)
        print('TestAllSets.results_dir={}'.format(TestAllSets.results_dir))
        if CORI_WEB:
            TestAllSets.results_dir = move_to_web(
                '/global/u1/f/(.*)/e3sm_diags', '/global/cfs/cdirs/acme/www/{}',
                TestAllSets.results_dir)

    def check_html_image(self, html_path, png_path):
        img_src = None
        option_value = None
        href = None
        with open(html_path, 'r') as html:
            for line in html:
                if not img_src:
                    re_str = '<img src="../../../../{}">'.format(png_path)
                    img_src = re.search(re_str, line)
                if not option_value:
                    re_str = '<option value="../../../../{}">'.format(png_path)
                    option_value = re.search(re_str, line)
                if not href:
                    re_str = 'href="../../../../{}">'.format(png_path)
                    href = re.search(re_str, line)
        self.assertIsNotNone(img_src)
        self.assertIsNotNone(option_value)
        self.assertIsNotNone(href)

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
                self.check_html_image(html_path, png_path)

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

    # Test results_dir
    def test_results_dir(self):
        self.assertTrue(TestAllSets.results_dir.endswith('all_sets_results_test/'))

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
            self.check_html_image(html_path, png_path)

    def test_cosp_histogram(self):
        self.check_plots_generic(
            set_name='cosp_histogram',
            case_id='MISR-COSP',
            ref_name='MISRCOSP',
            variables=['COSP_HISTOGRAM_MISR'],
            region='global'
        )

    def test_enso_diags_map(self):
        case_id = 'TREFHT-response'
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
            self.check_html_image(html_path, png_path)

    def test_enso_diags_scatter(self):
        case_id = 'TREFHT'
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
            self.check_html_image(html_path, png_path)

    def test_lat_lon(self):
        self.check_plots_plevs('lat_lon', 'global', [850.0])

    def test_meridional_mean_2d(self):
        self.check_plots_2d('meridional_mean_2d')

    def test_polar(self):
        self.check_plots_plevs('polar', 'polar_S', [850.0])

    def test_zonal_mean_2d(self):
        self.check_plots_2d('zonal_mean_2d')

    def test_zonal_mean_xy(self):
        self.check_plots_plevs('zonal_mean_xy', 'global', [200.0])


if __name__ == '__main__':
    unittest.main()
