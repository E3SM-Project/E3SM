import unittest
import cdms2
from acme_diags.derivations import acme


class TestACMEDerivations(unittest.TestCase):

    def test_convert_units(self):
        precc_file = cdms2.open('precc.nc')
        var = precc_file('PRECC')
        precc_file.close()

        new_var = acme._convert_units(var, 'mm/day')

        self.assertEqual(new_var.units, 'mm/day')

if __name__ == '__main__':
    unittest.main()
