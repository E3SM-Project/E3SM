import unittest
import os
import cdms2
from collections import OrderedDict
from acme_diags.derivations import acme as acme_derivations
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.driver.utils.dataset import Dataset

def get_abs_file_path(relative_path):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), relative_path)


class TestDataset(unittest.TestCase):
    def setUp(self):
        self.parameter = CoreParameter()

    def test_convert_units(self):
        with cdms2.open(get_abs_file_path('precc.nc')) as precc_file:
            var = precc_file('PRECC')

        new_var = acme_derivations.convert_units(var, 'mm/day')
        self.assertEqual(new_var.units, 'mm/day')

    def test_add_user_derived_vars(self):
        my_vars = {
            'A_NEW_VAR': {
                ('v1', 'v2'): lambda v1, v2: v1+v2,
                ('v3', 'v4'): lambda v3, v4: v3-v4
            },
            'PRECT': {
                ('MY_PRECT',): lambda my_prect: my_prect 
            }
        }
        self.parameter.derived_variables = my_vars
        data = Dataset(self.parameter, test=True)
        self.assertTrue('A_NEW_VAR' in data.derived_vars)

        # In the default my_vars, each entry
        # ('PRECT', 'A_NEW_VAR', etc) is an OrderedDict.
        # We must check that what the user inserted is
        # first, so it's used first.
        self.assertTrue(list(data.derived_vars['PRECT'].keys())[0] == ('MY_PRECT',))

    def test_is_timeseries(self):
        self.parameter.ref_timeseries_input = True
        data = Dataset(self.parameter, ref=True)
        self.assertTrue(data.is_timeseries())

        self.parameter.test_timeseries_input = True
        data = Dataset(self.parameter, test=True)
        self.assertTrue(data.is_timeseries())

        self.parameter.ref_timeseries_input = False
        data = Dataset(self.parameter, ref=True)
        self.assertFalse(data.is_timeseries())

        self.parameter.test_timeseries_input = False
        data = Dataset(self.parameter, test=True)
        self.assertFalse(data.is_timeseries())

    def test_is_climo(self):
        self.parameter.ref_timeseries_input = True
        data = Dataset(self.parameter, ref=True)
        self.assertFalse(data.is_climo())

        self.parameter.test_timeseries_input = True
        data = Dataset(self.parameter, test=True)
        self.assertFalse(data.is_climo())

        self.parameter.ref_timeseries_input = False
        data = Dataset(self.parameter, ref=True)
        self.assertTrue(data.is_climo())

        self.parameter.test_timeseries_input = False
        data = Dataset(self.parameter, test=True)
        self.assertTrue(data.is_climo())

    def test_get_attr_from_climo(self):
        # We pass in the path to a file, so the input directory
        # to the tests doesn't need to be like how it is for when e3sm_diags
        # is ran wit a bunch of diags.
        self.parameter.reference_data_path = './system'
        self.parameter.ref_file = "ta_ERA-Interim_ANN_198001_201401_climo.nc"        
        data = Dataset(self.parameter, ref=True)
        self.assertEqual(data.get_attr_from_climo('Conventions', 'ANN'), 'CF-1.0')



    '''
    def test_process_derived_var_passes(self):
        derived_var = {
            'PRECT': {
                ('pr'): lambda x: x,
                ('PRECC'): lambda x: 'this is some function'
            }
        }

        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        acme_derivations.process_derived_var('PRECT', derived_var,
                                 precc_file, self.parameter)
        precc_file.close()

    def test_process_derived_var_with_wrong_dict(self):
        # pr, nothing, and nothing2 are not variables in the file we open
        wrong_derived_var = {
            'PRECT': {
                ('pr'): lambda x: x,
                ('nothing1', 'nothing2'): lambda x, y: 'this is some function'
            }
        }

        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        with self.assertRaises(RuntimeError):
            acme_derivations.process_derived_var(
                'PRECT', wrong_derived_var, precc_file, self.parameter)
        precc_file.close()

    def test_process_derived_var_adds_to_dict(self):
        # the one that's usually in the parameters file
        derived_var_dict = {
            'PRECT': {('test'): lambda x: x}
        }
        # use this instead of the acme.derived_variables one
        default_derived_vars = {
            'PRECT': {
                ('pr'): lambda x: x,
                ('PRECC'): lambda x: 'this is some function'
            }
        }

        # add derived_var_dict to default_derived_vars
        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        self.parameter.derived_variables = derived_var_dict
        acme_derivations.process_derived_var(
            'PRECT', default_derived_vars, precc_file, self.parameter)
        precc_file.close()

        if 'test' not in default_derived_vars['PRECT']:
            self.fail("Failed to insert test derived variable")

        # make sure that test is the first one
        if 'test' != default_derived_vars['PRECT'].keys()[0]:
            self.fail(
                "Failed to insert test derived variable before default derived vars")

    def test_process_derived_var_adds_duplicate_to_dict(self):
        # the one that's usually in the parameters file
        # the function for PRECC below is different than the one in
        # default_derived_vars
        derived_var_dict = {
            'PRECT': {('PRECC'): lambda x: 'PRECC'}
        }
        # use this instead of the acme_derivations.derived_variables one
        default_derived_vars = {
            'PRECT': {
                ('pr'): lambda x: x,
                ('PRECC'): lambda x: 'this is some function'
            }
        }

        # add derived_var_dict to default_derived_vars
        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        self.parameter.derived_variables = derived_var_dict
        msg = acme_derivations.process_derived_var(
            'PRECT', default_derived_vars, precc_file, self.parameter)
        precc_file.close()
        if msg != 'PRECC':
            self.fail("Failed to insert a duplicate test derived variable")

    def test_process_derived_var_works_with_ordereddict(self):
        derived_var_dict = {
            'PRECT': OrderedDict([
                (('something'), lambda x: 'something')
            ])
        }

        default_derived_vars = {
            'PRECT': OrderedDict([
                (('pr'), lambda x: x),
                (('PRECC'), lambda x: 'this is some function')
            ])
        }

        precc_file = cdms2.open(get_abs_file_path('precc.nc'))
        self.parameter.derived_variables = derived_var_dict
        acme_derivations.process_derived_var(
            'PRECT', default_derived_vars, precc_file, self.parameter)
        precc_file.close()
        # Check that 'something' was inserted first
        self.assertEqual(['something', 'pr', 'PRECC'],
                         default_derived_vars['PRECT'].keys())

    def test_mask_by(self):
        with cdms2.open(get_abs_file_path('precc.nc')) as precc_file:
            prcc = precc_file('PRECC')
        with cdms2.open(get_abs_file_path('precl.nc')) as precl_file:
            prcl = precl_file('PRECL')

        acme_derivations.mask_by(prcc, prcl, low_limit=2.0)
    '''

if __name__ == '__main__':
    unittest.main()
