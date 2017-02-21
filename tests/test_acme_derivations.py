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

    def test_process_derived_var_passes(self):
        derived_var = {
            'PRECT': [
                (['pr'],  lambda x: x),
                (['PRECC'], lambda x: 'this is some function')
            ],
        }

        precc_file = cdms2.open('precc.nc')
        acme.process_derived_var('PRECT', derived_var, precc_file)
        precc_file.close()

    def test_process_derived_var_with_wrong_dict(self):
        # pr, nothing, and nothing2 are not variables in the file we open
        wrong_derived_var = {
            'PRECT': [
                (['pr'],  lambda x: x),
                (['nothing','nothing2'], lambda x, y: 'this is some function')
            ],
        }

        precc_file = cdms2.open('precc.nc')
        with self.assertRaises(RuntimeError):
            acme.process_derived_var('PRECT', wrong_derived_var, precc_file)
        precc_file.close()


if __name__ == '__main__':
    unittest.main()
