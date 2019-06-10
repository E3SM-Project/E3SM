import unittest
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.acme_diags_driver import _get_default_diags


class TestACMEDriver(unittest.TestCase):
    def setUp(self):
        self.parameter = CoreParameter()

    def test_get_default_diags(self):
        diags_pth = _get_default_diags('5', 'AMWG')
        self.assertTrue('lat_lon_AMWG' in diags_pth)


if __name__ == '__main__':
    unittest.main()
