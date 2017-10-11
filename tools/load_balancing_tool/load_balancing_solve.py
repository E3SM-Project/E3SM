#!/usr/bin/env python
"""
Reads timing data created with load_balancing_submit.py (or otherwise,
see --timing_files option) and solves an mixed integer optimization problem
using these timings. The default layout (IceLndAtmOcn) minimizes the cost per
model day assuming the layout:
              ____________________
             | ICE  |  LND  |     |
             |______|_______|     |
             |              | OCN |
             |    ATM       |     |
             |______________|_____|

It is possible to extend this tool to solve for other layouts.
"""
import re
import json

try:
    from Tools.standard_script_setup import *
except ImportError, e:
    print "Error importing Tools.standard_script_setup"
    print "May need to add cime/scripts to PYTHONPATH\n"
    raise ImportError(e)

from CIME.utils import expect
from CIME.XML.machines import Machines
logger = logging.getLogger(__name__)

# These values can be overridden on the command line
DEFAULT_TESTID = "lbt"
DEFAULT_BLOCKSIZE = 1
DEFAULT_LAYOUT = "IceLndAtmOcn"
COMPONENT_LIST = ['ATM', 'ICE', 'CPL', 'LND', 'WAV', 'ROF', 'OCN', 'GLC', 'ESP']

###############################################################################
def parse_command_line(args, description):
###############################################################################
    help_str = """
    Solve a Mixed Integer Linear Program to find a PE layout that minimizes
    the wall-clock time per model day.
    """
    parser = argparse.ArgumentParser(usage=help_str,
                                     description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument('--test-id', default=DEFAULT_TESTID,
                        help='test-id to use for all timing runs')

    parser.add_argument("-r", "--test-root",
                        help="Where test cases were created."
                        " Will default to output root as defined in the config_machines file")

    parser.add_argument('--timing-dir', help='alternative to using casename '
                        'to find timing data, instead read all files in'
                        ' this directory')

    parser.add_argument('--blocksize',
                        help='default minimum size of blocks to assign to all '
                        'components. Components can be assigned different '
                        'blocksizes using --blocksize_XXX. Default 1', type=int)

    for c in COMPONENT_LIST:
        parser.add_argument('--blocksize-%s' % c.lower(),
                            help='minimum blocksize for component %s, if '
                            'different from --blocksize', type=int)

    parser.add_argument('--total-tasks', type=int,
                        help='Number of pes available for assignment')

    parser.add_argument("--layout",
                        help="name of layout to solve (default selected internally)")

    parser.add_argument("--graph-models", action="store_true",
                        help="plot cost v. ntasks models. requires matplotlib")

    parser.add_argument("--print-models", action="store_true",
                        help="print all costs and ntasks")

    parser.add_argument("--pe-output", help="write pe layout to file")

    parser.add_argument('--json-output', help="write MILP data to .json file")

    parser.add_argument('--json-input', help="solve using data from .json file")

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args,
                                                                     parser)
    if args.total_tasks is None and args.json_input is None:
        expect(args.total_tasks is not None or args.json_input is not None,
               "--total-tasks or --json-input option must be set")

    blocksizes = {}
    for c in COMPONENT_LIST:
        attrib = 'blocksize_%s' % c.lower()
        if getattr(args, attrib) is not None:
            blocksizes[c] = getattr(args, attrib)
        elif args.blocksize is not None:
            blocksizes[c] = args.blocksize
    test_root = args.test_root
    if test_root is None:
        machobj = Machines()
        test_root = machobj.get_value("CIME_OUTPUT_ROOT")

    return (args.test_id, test_root, args.timing_dir, blocksizes,
            args.total_tasks, args.layout, args.graph_models,
            args.print_models, args.pe_output, args.json_output,
            args.json_input)


def _locate_timing_files(test_root, test_id, timing_dir):
    """
    Find all possible directories for timing files
    """
    timing_files = []
    timing_cases_tmp = []
    timing_dirs = []

    # Add command-line timing directory if it exists
    if timing_dir is not None:
        logger.info('found directory ' + timing_dir)
        timing_dirs.append(timing_dir)
    else:
        # Add script_dir/casename_prefix_*/timing
        for fn in os.listdir(test_root):
            if fn.endswith(test_id):
                fn = os.path.join(test_root, fn, "timing")
                if os.path.isdir(fn):
                    print "found {}".format(fn)
                    timing_cases_tmp.append(fn)
        timing_dirs = sorted(timing_cases_tmp)

    # Now add all non-.gz files in the directories to be read in
    for td in timing_dirs:
        full_fn = None
        for fn in os.listdir(td):
            full_fn = os.path.join(td, fn)
            if full_fn.find('.gz') < 0:
                timing_files.append(full_fn)
        if full_fn is None:
            logger.warning("WARNING: no timing files found in directory %s", (td))
    return timing_files

def _parse_timing_files(timing_files):
    """
    Parse every file in list for timing information and return data dict
    """
    data = {}
    for timing_file in timing_files:
        timing = _read_timing_file(timing_file)
        logger.debug('ntasks: %s' % "; ".join([str(k) + ":" +
                                               str(timing[k]['ntasks'])
                                               for k in timing.keys()]))
        logger.debug('cost: %s' % "; ".join([str(k) + ":" +
                                             str(timing[k]['cost'])
                                             for k in timing.keys()]))
        for key in timing:
            if key not in data:
                data[key] = {'cost':[], 'ntasks':[], 'nthrds':[]}

            if timing[key]['ntasks'] in data[key]['ntasks']:
                logger.warning('WARNING: duplicate timing run data in %s '
                               'for %s ntasks=%d.', timing_file, key,
                               timing[key]['ntasks'])
                index = data[key]['ntasks'].index(timing[key]['ntasks'])
                logger.warning('Existing value: cost=%s. Ignoring new value: '
                               'cost=%s', data[key]['cost'][index],
                               timing[key]['cost'])
            elif 'name' in data[key] and data[key]['name'] != timing[key]['name']:
                expect(False, "Timing files have inconsistant model components {} has {} vs {}"
                       .format(key, data[key]['name'], timing[key]['name']))
            else:
                data[key]['name'] = timing[key]['name']
                data[key]['cost'].append(timing[key]['cost'])
                data[key]['ntasks'].append(timing[key]['ntasks'])
                data[key]['nthrds'].append(timing[key]['nthrds'])
    return data

def _set_blocksizes(data, blocksizes):
    """
    Set blocksizes according to command line arguments.
    Specific command line arguments override current data, but
    do not set to default if it already exists
    """
    for key in COMPONENT_LIST:
        if key in data:
            if key in blocksizes:
                data[key]['blocksize'] = blocksizes[key]
            elif 'blocksize' not in data[key]:
                data[key]['blocksize'] = DEFAULT_BLOCKSIZE

def _read_timing_file(filename):
    """
    Read in timing files to get the costs (time/mday) for each test

    return model dictionaries. Example
    {'ICE':{'ntasks':8,'nthrds':1,'cost':40.6},
     'ATM':{'ntasks':8,'nthrds':1,'cost':120.4},
     ...
    }
    """

    logger.info('Reading timing file %s', filename)
    try:
        timing_file = open(filename, "r")
        timing_lines = timing_file.readlines()
        timing_file.close()
    except Exception, e:
        logger.critical("Unable to open file %s", filename)
        raise e
    models = {}
    for line in timing_lines:
        # Get number of tasks and thrds
        #  atm = xatm       8      0         8      x     1    1  (1 )
        #(\w+) = (\w+) \s+ \d+ \s+ \d+ \s+ (\d+)\s+ x\s+(\d+)
        m = re.search(r"(\w+) = (\w+)\s+\d+\s+\d+\s+(\d+)\s+x\s+(\d+)", line)
        if m:
            component = m.groups()[0].upper()
            name = m.groups()[1].upper()
            ntasks = int(m.groups()[2])
            nthrds = int(m.groups()[3])
            if component in models:
                models[component]['ntasks'] = ntasks
                models[component]['nthrds'] = nthrds
                models[component]['name'] = name
            else:
                models[component] = {'name':name,'ntasks':ntasks, 'nthrds':nthrds}
            continue

        # get cost
        # ATM Run Time:      17.433 seconds        1.743 seconds/mday
        #(\w+)Run Time: \s  \d+.\d+ seconds \s+(\d+.\d+) seconds/mday
        m = re.search(r"(\w+) Run Time:\s+(\d+\.\d+) seconds \s+(\d+\.\d+)"
                      " seconds/mday", line)
        if m:
            component = m.groups()[0]
            cost = float(m.groups()[1])
            if component != "TOT":
                if component in models:
                    models[component]['cost'] = cost
                else:
                    models[component] = {'cost':cost}
    return models

################################################################################
def load_balancing_solve(test_id, test_root, timing_dir, blocksizes, total_tasks, layout, graph_models, print_models, pe_output, json_output, json_input):
################################################################################
    if json_input is not None:
        # All data is read from given json file
        with open(json_input, "r") as jsonfile:
            try:
                data = json.load(jsonfile)
            except ValueError, e:
                logger.critical("Unable to parse json file %s", jsonfile)
                raise e
        # layout, totaltasks, blocksizes may already be set by json file
        # but can be overriden by options
        if layout is not None:
            data['layout'] = layout
        if total_tasks is not None:
            data['totaltasks'] = total_tasks

    else:
        # find and parse timing files
        timing_files = _locate_timing_files(test_root,
                                            test_id,
                                            timing_dir)

        expect(len(timing_files) > 0, "No timing data found")

        data = _parse_timing_files(timing_files)

        data['totaltasks'] = total_tasks
        if layout is None:
            # try to determine layout automatically
            if 'ATM' in data and 'OCN' in data and 'WAV' in data:
                aname = data['ATM']['name']
                oname = data['OCN']['name']
                wname = data['WAV']['name']
                if aname not in ('DATM', 'XATM', 'SATM') and \
                   oname not in ('DOCN', 'XOCN', 'SOCN'):
                    if wname in ('DWAV', 'XWAV', 'SWAV'):
                        data['layout'] = "IceLndAtmOcn"
                    else:
                        data['layout'] = "IceLndWavAtmOcn"

                    logger.info("Using layout  = {}".format(data['layout']))
                else:
                    expect(False, "Could not automatically determine layout")
        else:
            data['layout'] = layout

    _set_blocksizes(data, blocksizes)

    # Allow dumping to json file before trying to load optimization
    if json_output is not None:
        logger.info("Writing MILP data to %s", json_output)
        with open(json_output, "w") as outfile:
            json.dump(data, outfile, indent=4)

    import optimize_model

    # Use atm-lnd-ocn-ice linear program
    opt = optimize_model.solver_factory(data)
    opt.optimize()
    if graph_models:
        opt.graph_costs()
    if print_models:
        opt.write_timings(fd=None, level=logging.INFO)
    else:
        opt.write_timings(fd=None, level=logging.DEBUG)

    logger.info("Solving Mixed Integer Linear Program using PuLP interface to "
                "COIN-CBC")

    status = opt.optimize()
    logger.info("PuLP solver status: " + opt.get_state_string(status))
    solution = opt.get_solution()
    for k in sorted(solution):
        if k[0] == 'N':
            logger.info("%s = %d", k, solution[k])
        else:
            logger.info("%s = %f", k, solution[k])

    if pe_output:
        opt.write_pe_file(pe_output)

    return 0

###############################################################################
def _main_func(description):
###############################################################################
    test_id, test_root, timing_dir, blocksizes, total_tasks, layout, graph_models, print_models, pe_output, json_output, json_input = parse_command_line(sys.argv, description)

    sys.exit(load_balancing_solve(test_id, test_root, timing_dir, blocksizes, total_tasks, layout, graph_models, print_models, pe_output, json_output, json_input))

###############################################################################

if __name__ == "__main__":
    _main_func(__doc__)
