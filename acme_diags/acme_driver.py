from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_parameter import ACMEParameter
from acme_diags.acme_viewer import create_viewer

def make_parameters(orginal_parameter):
    """ Create multiple parameters given a list of
    parameters in a json and an original parameter """

    if hasattr(original_parameter, 'custom_diags'):
        f_data = open(original_parameter.custom_diags).read()
    else:
        pth = os.path.join(sys.prefix, 'share', 'acme_diags', 'set5', 'set5_diags_AMWG_default.json')
        f_data = open(pth).read()
    json_file = json.loads(f_data)

    parameters = []
    for key in json_file:
        for single_run in json_file[key]:
            p = ACMEParameter()
            for attr_name in single_run:
                setattr(p, attr_name, single_run[attr_name])

            # Add attributes of original_parameter to p
            p.__dict__.update(orginal_parameter.__dict__)
            parameters.append(p)
    return parameters


parser = ACMEParser()
original_parameter = parser.get_parameter(default_vars=False)
if not hasattr(original_parameter, 'results_dir'):
    import datetime
    dt = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    original_parameter.results_dir = '{}-{}'.format('acme_diags_results', dt)

parameters = make_parameters(original_parameter)

for parameter in parameters:
    if not hasattr(parameter, 'set'):
        raise RuntimeError('Parameter needs to have an attribute set')
    for pset in parameter.set:
        if pset == '5':
            from acme_diags.driver.set5_driver import compute
        elif pset == '7':
            from acme_diags.driver.set7_driver import compute
        compute(parameter)

create_viewer(original_parameter.results_dir, parameters, 'png')
