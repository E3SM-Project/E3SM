import unittest
import cdms2
from acme_diags.derivations import acme


class TestACMEDerivations(unittest.TestCase):

    def test_convert_units(self):
        precc_file = cdms2.open('precc.nc')
        var1 = precc_file('PRECC')
        precc_file.close()

        precl_file = cdms2.open('precl.nc')
        var2 = precl_file('PRECL')
        precl_file.close()

        new_var1, new_var2 = acme._convert_units(var1, var2, 'mm/day')

        self.assertEqual(new_var1.units, new_var2.units)

if __name__ == '__main__':
    unittest.main()
