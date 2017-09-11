import unittest
import os
import subprocess
import shlex


class TestFlake8(unittest.TestCase):

    def test_flake8(self):
        # E402 is ignore b/c we do logic between import
        # statements so numpy can run without X
        pth = os.path.dirname(__file__)
        pth = os.path.join(pth, '..')
        pth = os.path.abspath(pth)
        cmd = 'flake8 {} --ignore=E402 --max-line-length=140'.format(pth)

        process = subprocess.Popen(shlex.split(cmd),
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        out = process.stdout.read()
        if out != "":
            print(out)
        self.assertEqual(out, "")


if __name__ == '__main__':
    unittest.main()
