import unittest
import subprocess
import shlex
import os
import shutil


class TestAllSets(unittest.TestCase):

    def get_viewer_dir(self, output):
        """Given output from acme_diags_driver,
        extract the path to results_dir."""
        for line in output.splitlines():
            if 'viewer' in line.lower() and 'html' in line.lower():
                viewer_dir = line
                # get the /path/to/results_dir/viewer/index.html
                viewer_dir = viewer_dir.split()[-1]
                # get the /path/to/results_dir/
                viewer_dir = viewer_dir.split('viewer')[0]
                return viewer_dir
        raise RuntimeError('Could not find the viewer directory from: {}'.format(output))

    def count_images(self, directory, file_type='png'):
        """Count the number of images of type file_type in directory"""
        count = 0
        for _, __, files in os.walk(directory):
            for f in files:
                if f.endswith(file_type):
                    count += 1
        return count

    def setUp(self):
        # The number of diag runs in all_sets.cfg
        self.NUM_OF_DIAGS = 5

    def test_all_sets_mpl(self):
        pth = os.path.dirname(__file__)
        cfg_pth = os.path.join(pth, 'all_sets.cfg')
        cfg_pth = os.path.abspath(cfg_pth)

        # *_data_path needs to be added b/c the tests runs the diags from a different location
        cmd = 'acme_diags_driver.py -d {} --reference_data_path {} --test_data_path {}'
        cmd = cmd.format(cfg_pth, pth, pth)
        process = subprocess.Popen(shlex.split(cmd),
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        process.wait()
        out = process.stdout.read()

        # count the number of pngs in viewer_dir
        viewer_dir = self.get_viewer_dir(out)
        count = self.count_images(viewer_dir)
        # -1 is needed because of the E3SM logo in the viewer html
        self.assertEqual(count - 1, self.NUM_OF_DIAGS)

        shutil.rmtree(viewer_dir)  # remove all generated results from the diags

    def test_all_sets_vcs(self):
        pth = os.path.dirname(__file__)
        cfg_pth = os.path.join(pth, 'all_sets.cfg')
        cfg_pth = os.path.abspath(cfg_pth)

        # *_data_path needs to be added b/c the tests runs the diags from a different location
        cmd = 'acme_diags_driver.py -d {} --backend vcs --reference_data_path {} --test_data_path {}'
        cmd = cmd.format(cfg_pth, pth, pth)
        process = subprocess.Popen(shlex.split(cmd),
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        process.wait()
        out = process.stdout.read()

        # count the number of pngs in viewer_dir
        viewer_dir = self.get_viewer_dir(out)
        count = self.count_images(viewer_dir)
        # -1 is needed because of the E3SM logo in the viewer html
        self.assertEqual(count - 1, self.NUM_OF_DIAGS)

        shutil.rmtree(viewer_dir)  # remove all generated results from the diags


if __name__ == '__main__':
    unittest.main()
