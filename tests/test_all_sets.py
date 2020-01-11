import unittest
import os
import re
import shlex
import shutil
import subprocess

def count_images(directory, file_type='png'):
    """Count the number of images of type file_type in directory"""
    count = 0
    for _, __, files in os.walk(directory):
        for f in files:
            if f.endswith(file_type):
                count += 1
    return count

class TestAllSets(unittest.TestCase):

    def get_results_dir(self, output):
        """Given output from acme_diags_driver, extract the path to results_dir."""
        for line in output:
            match = re.search('Viewer HTML generated at (.*)viewer.*.html', line)
            if match:
                results_dir = match.group(1)
                return results_dir
        self.fail('No viewer directory listed in output: {}'.format(output))

    def run_test(self, backend):
        pth = os.path.dirname(__file__)
        pth = os.path.join(pth, 'system')
        cfg_pth = os.path.join(pth, 'all_sets_modified.cfg')
        cfg_pth = os.path.abspath(cfg_pth)

        if backend == 'mpl':
            backend_option = ''
            expected_num_diags = 8
        elif backend == 'vcs':
            backend_option = ' --backend vcs'
            expected_num_diags = 7
        else:
            raise RuntimeError('Invalid backend: {}'.format(backend))
        # *_data_path needs to be added b/c the tests runs the diags from a different location
        cmd = 'acme_diags_driver.py -d {}{} --reference_data_path {} --test_data_path {}'
        cmd = cmd.format(cfg_pth, backend_option, pth, pth)
        # This raises a CalledProcessError if cmd has a non-zero return code.
        out = subprocess.check_output(cmd.split()).decode('utf-8').splitlines()

        # count the number of pngs in viewer_dir
        results_dir = self.get_results_dir(out)
        count = count_images(results_dir)
        # -1 is needed because of the E3SM logo in the viewer html
        self.assertEqual(count - 1, expected_num_diags)

        shutil.rmtree(results_dir) # remove all generated results from the diags

    def test_all_sets_mpl(self):
        self.run_test('mpl')

    def test_all_sets_vcs(self):
        self.run_test('vcs')

if __name__ == '__main__':
    unittest.main()
