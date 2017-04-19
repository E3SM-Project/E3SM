import unittest
import os
import cdms2
from acme_diags.derivations import acme
from acme_diags.acme_parameter import ACMEParameter

def get_abs_file_path(relative_path):
    return os.path.dirname(os.path.abspath(__file__)) + '/' + relative_path

class TestACMEDerivations(unittest.TestCase):
    def setUp(self):
        self.parameter = ACMEParameter()

    def test_convert_units(self):
        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
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

        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        acme.process_derived_var('PRECT', derived_var, precc_file, self.parameter)
        precc_file.close()

    def test_process_derived_var_with_wrong_dict(self):
        # pr, nothing, and nothing2 are not variables in the file we open
        wrong_derived_var = {
            'PRECT': [
                (['pr'],  lambda x: x),
                (['nothing','nothing2'], lambda x, y: 'this is some function')
            ],
        }

        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        with self.assertRaises(RuntimeError):
            acme.process_derived_var('PRECT', wrong_derived_var, precc_file, self.parameter)
        precc_file.close()

    def test_that_process_derived_var_adds_to_dict(self):
        # the one that's usually in the parameters file
        derived_var_dict = {
            'PRECT': [(['test'], lambda x: x)]
        }
        # use this instead of the acme.derived_variables one
        default_derived_vars = {
            'PRECT': [
                (['pr'],  lambda x: x),
                (['PRECC'], lambda x: 'this is some function')
            ],
        }

        # add derived_var_dict to default_derived_vars
        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        self.parameter.derived_variables = derived_var_dict
        acme.process_derived_var('PRECT', default_derived_vars, precc_file, self.parameter)
        precc_file.close()

        # check that (['test'], lambda x: x) was inserted correctly
        for derived_var in default_derived_vars['PRECT']:
            if ['test'] == derived_var[0]:
                return  # everything is okay
        self.fail("Failed to insert test derived variable")

    def test_mask_by_passes(self):
        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        prcc = precc_file('PRECC')
        precl_file = cdms2.open(get_abs_file_path('precl.nc'))
        prcl = precl_file('PRECL')
        
        acme.mask_by(prcc, prcl, low_limit = 2.0)
        precc_file.close()
        precl_file.close()

if __name__ == '__main__':
    unittest.main()
