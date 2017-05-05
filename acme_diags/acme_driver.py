import os
import sys
import json
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_parameter import ACMEParameter
from acme_diags.acme_viewer import create_viewer

def _get_default_diags(set_num):
    """Get the data from the json corresponding to set_num"""
    folder = 'set{}'.format(set_num)
    fnm = 'set{}_diags_AMWG_default.json'.format(set_num)
    pth = os.path.join(sys.prefix, 'share', 'acme_diags', folder, fnm)
    with open(pth) as json_file:
        json_data = json.loads(json_file.read())
    return json_data

def make_parameters(original_parameter):
    """Create multiple parameters given a list of
    parameters in a json and an original parameter"""

    if hasattr(original_parameter, 'custom_diags'):
        with open(original_parameter.custom_diags) as json_file:
            json_data = json.loads(json_file.read())
    else:
        json_data = {'': []}  # the first key doesn't hold any value, so it's ''
        for set_num in original_parameter.set:
            default_set_runs = _get_default_diags(set_num)
            for _, set_runs in default_set_runs.iteritems():
                for single_run in set_runs:
                    json_data[''].append(single_run)
    parameters = []
    for key in json_data:
        for single_run in json_data[key]:
            p = ACMEParameter()
            for attr_name in single_run:
                setattr(p, attr_name, single_run[attr_name])

            # Add attributes of original_parameter to p
            p.__dict__.update(original_parameter.__dict__)
            p.check_values()
            parameters.append(p)
    return parameters

if __name__ == '__main__':
    parser = ACMEParser()
    original_parameter = parser.get_parameter(default_vars=False)
    if not hasattr(original_parameter, 'results_dir'):
        import datetime
        dt = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        original_parameter.results_dir = '{}-{}'.format('acme_diags_results', dt)

    if not os.path.exists(original_parameter.results_dir):
        os.makedirs(original_parameter.results_dir)
    parameters = make_parameters(original_parameter)

    for parameter in parameters:
        if not hasattr(parameter, 'set'):
            raise RuntimeError('Parameter needs to have an attribute set')
        for pset in parameter.set:
            pset = str(pset)
            if pset == '5':
                from acme_diags.driver.set5_driver import compute
            elif pset == '7':
                from acme_diags.driver.set7_driver import compute
            else:
                print('Plot set {} is not supported yet. Please give us time.'.format(pset))
                quit()
            compute(parameter)

    create_viewer(original_parameter.results_dir, parameters, 'png')
