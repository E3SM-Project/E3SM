#!/usr/bin/env python
"""
Script to submit a series of ACME runs to get data for
time vs nprocessors model. This data will be used to generate
a processor layout that achieves high efficiency
"""
from xml.etree.ElementTree import ParseError
import shutil

try:
    from Tools.standard_script_setup import *
except ImportError, e:
    print 'Error importing Tools.standard_script_setup'
    print 'May need to add cime/scripts to PYTHONPATH\n'
    raise ImportError(e)

from CIME.utils import expect, get_full_test_name
from CIME.case import Case
from CIME.XML.pes import Pes
from CIME.XML.machines import Machines
from CIME.test_scheduler import TestScheduler

logger = logging.getLogger(__name__)


# Default CIME variables, these can be overridden using the
# --extra-options-file option
CIME_DEFAULTS = {
    'STOP_OPTION':'ndays',
    'STOP_N':'10',
    'REST_OPTION':'never',
    'DOUT_S':'FALSE',
    'COMP_RUN_BARRIERS':'TRUE',
    'TIMER_LEVEL':'9'
}

DEFAULT_TESTID = 'lbt'

###############################################################################
def parse_command_line(args, description):
###############################################################################
    help_str = """
Requires a pes xml file listing the timing runs you will submit and
their corresponding pe layouts. Use the 'pesize' tag to name each run.

After running this submission tool, run load_balancing_solve.py to
solve the mixed integer linear program minimizing the simulation time.

example_pes.xml:
<config_pes>
  <grid name="any">
    <mach name="any">
      <pes compset="any" pesize="0">
        <comment>none</comment>
        <ntasks>
          <ntasks_atm>8</ntasks_atm>
          <ntasks_lnd>8</ntasks_lnd>
          <ntasks_rof>8</ntasks_rof>
          <ntasks_ice>8</ntasks_ice>
          <ntasks_ocn>8</ntasks_ocn>
          <ntasks_glc>8</ntasks_glc>
          <ntasks_wav>8</ntasks_wav>
          <ntasks_cpl>8</ntasks_cpl>
        </ntasks>
        <nthrds>
          <nthrds_atm>1</nthrds_atm>
          <nthrds_lnd>1</nthrds_lnd>
          <nthrds_rof>1</nthrds_rof>
          <nthrds_ice>1</nthrds_ice>
          <nthrds_ocn>1</nthrds_ocn>
          <nthrds_glc>1</nthrds_glc>
          <nthrds_wav>1</nthrds_wav>
          <nthrds_cpl>1</nthrds_cpl>
        </nthrds>
        <rootpe>
          <rootpe_atm>0</rootpe_atm>
          <rootpe_lnd>0</rootpe_lnd>
          <rootpe_rof>0</rootpe_rof>
          <rootpe_ice>0</rootpe_ice>
          <rootpe_ocn>0</rootpe_ocn>
          <rootpe_glc>0</rootpe_glc>
          <rootpe_wav>0</rootpe_wav>
          <rootpe_cpl>0</rootpe_cpl>
        </rootpe>
      </pes>
    </mach>
  </grid>

  <grid name="any">
    <mach name="any">
      <pes compset="any" pesize="1">
        <comment>none</comment>
        <ntasks>
          <ntasks_atm>32</ntasks_atm>
          <ntasks_lnd>32</ntasks_lnd>
          <ntasks_rof>32</ntasks_rof>
          <ntasks_ice>32</ntasks_ice>
          <ntasks_ocn>32</ntasks_ocn>
          <ntasks_glc>32</ntasks_glc>
          <ntasks_wav>32</ntasks_wav>
          <ntasks_cpl>32</ntasks_cpl>
        </ntasks>
        <nthrds>
          <nthrds_atm>1</nthrds_atm>
          <nthrds_lnd>1</nthrds_lnd>
          <nthrds_rof>1</nthrds_rof>
          <nthrds_ice>1</nthrds_ice>
          <nthrds_ocn>1</nthrds_ocn>
          <nthrds_glc>1</nthrds_glc>
          <nthrds_wav>1</nthrds_wav>
          <nthrds_cpl>1</nthrds_cpl>
        </nthrds>
        <rootpe>
          <rootpe_atm>0</rootpe_atm>
          <rootpe_lnd>0</rootpe_lnd>
          <rootpe_rof>0</rootpe_rof>
          <rootpe_ice>0</rootpe_ice>
          <rootpe_ocn>0</rootpe_ocn>
          <rootpe_glc>0</rootpe_glc>
          <rootpe_wav>0</rootpe_wav>
          <rootpe_cpl>0</rootpe_cpl>
        </rootpe>
      </pes>
    </mach>
  </grid>

  <grid name="any">
    <mach name="any">
      <pes compset="any" pesize="2">
        <comment>none</comment>
        <ntasks>
          <ntasks_atm>128</ntasks_atm>
          <ntasks_lnd>128</ntasks_lnd>
          <ntasks_rof>128</ntasks_rof>
          <ntasks_ice>128</ntasks_ice>
          <ntasks_ocn>128</ntasks_ocn>
          <ntasks_glc>128</ntasks_glc>
          <ntasks_wav>128</ntasks_wav>
          <ntasks_cpl>128</ntasks_cpl>
        </ntasks>
        <nthrds>
          <nthrds_atm>1</nthrds_atm>
          <nthrds_lnd>1</nthrds_lnd>
          <nthrds_rof>1</nthrds_rof>
          <nthrds_ice>1</nthrds_ice>
          <nthrds_ocn>1</nthrds_ocn>
          <nthrds_glc>1</nthrds_glc>
          <nthrds_wav>1</nthrds_wav>
          <nthrds_cpl>1</nthrds_cpl>
        </nthrds>
        <rootpe>
          <rootpe_atm>0</rootpe_atm>
          <rootpe_lnd>0</rootpe_lnd>
          <rootpe_rof>0</rootpe_rof>
          <rootpe_ice>0</rootpe_ice>
          <rootpe_ocn>0</rootpe_ocn>
          <rootpe_glc>0</rootpe_glc>
          <rootpe_wav>0</rootpe_wav>
          <rootpe_cpl>0</rootpe_cpl>
        </rootpe>
      </pes>
    </mach>
  </grid>
</config_pes>
"""
    parser = argparse.ArgumentParser(usage=help_str,
                                     description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    CIME.utils.setup_standard_logging_options(parser)
    # Required arguments
    parser.add_argument('--compset',
                        help='Specify compset', required=True)
    parser.add_argument('--res',
                        help='Specify resolution', required=True)
    parser.add_argument('--pesfile', required=True)

    # Optional pass-through arguments to create_newcase
    parser.add_argument('--compiler', help='Choose compiler to build with')

    parser.add_argument('--project', help='Specify project id')

    parser.add_argument('--machine', help='machine name')

    parser.add_argument('--mpilib', help='mpi library name')

    parser.add_argument("-r", "--test-root",
                        help="Where test cases will be created."
                        " Will default to output root as defined in the config_machines file")

    parser.add_argument('--extra-options-file',
                        help='file listing options to be run using xmlchange')
    parser.add_argument('--test-id', default=DEFAULT_TESTID,
                        help='test-id to use for all timing runs')
    parser.add_argument('--force-purge', action='store_true')

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args, parser)

    return (args.compset, args.res, args.pesfile, args.mpilib,
            args.compiler, args.project, args.machine, args.extra_options_file,
            args.test_id, args.force_purge, args.test_root)

################################################################################
def load_balancing_submit(compset, res, pesfile, mpilib, compiler, project, machine,
                          extra_options_file, test_id, force_purge, test_root):
################################################################################
    # Read in list of pes from given file
    expect(os.access(pesfile, os.R_OK), 'ERROR: File %s not found', pesfile)

    logger.info('Reading XML file %s. Searching for pesize entries:', pesfile)
    try:
        pesobj = Pes(pesfile)
    except ParseError:
        expect(False, 'ERROR: File %s not parseable', pesfile)

    pesize_list = []
    grid_nodes = pesobj.get_children("grid")
    for gnode in grid_nodes:
        mach_nodes = pesobj.get_children("mach", root=gnode)
        for mnode in mach_nodes:
            pes_nodes = pesobj.get_children("pes", root=mnode)
            for pnode in pes_nodes:
                pesize = pesobj.get(pnode, 'pesize')
                if not pesize:
                    logger.critical('No pesize for pes node in file %s', pesfile)
                if pesize in pesize_list:
                    logger.critical('pesize %s duplicated in file %s', pesize, pesfile)
                pesize_list.append(pesize)

    expect(pesize_list, 'ERROR: No grid entries found in pes file {}'.format(pesfile))

    machobj = Machines(machine=machine)
    if test_root is None:
        test_root = machobj.get_value("CIME_OUTPUT_ROOT")
    if machine is None:
        machine = machobj.get_machine_name()
        print "machine is {}".format(machine)
    if compiler is None:
        compiler = machobj.get_default_compiler()
        print "compiler is {}".format(compiler)
    if mpilib is None:
        mpilib = machobj.get_default_MPIlib({"compiler":compiler})




    test_names = []
    for i in xrange(len(pesize_list)):
        test_names.append(get_full_test_name("PFS_I{}".format(i),grid=res, compset=compset,
                                             machine=machine, compiler=compiler))
        casedir = os.path.join(test_root, test_names[-1] + "." + test_id)
        print "casedir is {}".format(casedir)
        if os.path.isdir(casedir):
            if force_purge:
                logger.info('Removing directory %s', casedir)
                shutil.rmtree(casedir)
            else:
                expect(False,
                       "casedir {} already exists, use the --force-purge option, --test-root or"
                       " --test-id options".format(casedir))

    tests = TestScheduler(test_names, no_setup = True,
                          compiler=compiler, machine_name=machine, mpilib=mpilib,
                          test_root=test_root, test_id=test_id, project=project)
    success = tests.run_tests(wait=True)
    expect(success, "Error in creating cases")
    testnames = []
    for test in tests.get_testnames():
        testname =  os.path.join(test_root, test + "." + test_id)
        testnames.append( testname)
        logger.info("test is {}".format(testname))
        with Case(testname) as case:
            pes_ntasks, pes_nthrds, pes_rootpe, _, _, _ = \
                                                    pesobj.find_pes_layout('any', 'any', 'any', pesize_opts=pesize_list.pop(0))
            for key in pes_ntasks:
                case.set_value(key, pes_ntasks[key])
            for key in pes_nthrds:
                case.set_value(key, pes_nthrds[key])
            for key in pes_rootpe:
                case.set_value(key, pes_rootpe[key])

            if extra_options_file is not None:
                try:
                    extras = open(extra_options_file, 'r')
                    for line in extras.readlines():
                        split = line.split('=')
                        if len(split) == 2:
                            logger.info('setting %s=%s', split[0], split[1])
                            case.set_value(split[0], split[1])
                        else:
                            logger.debug('ignoring line in {}: {}'.format(
                                extra_options_file, line))
                    extras.close()
                except IOError:
                    expect(False, "ERROR: Could not read file {}".format(extra_options_file))


    tests = TestScheduler(test_names, use_existing=True, test_root=test_root, test_id=test_id)
    success = tests.run_tests(wait=False)
    expect(success, "Error in running cases")

    # need to fix
    logger.info('Timing jobs submitted. After jobs completed, run to optimize '
                'pe layout:\n  load_balancing_solve --test-id {} --test-root {}'.
                format(test_id, test_root))

###############################################################################
def _main_func(description):
###############################################################################
    compset, res, pesfile, mpilib, compiler, project, machine, extra_options_file, casename_prefix,  \
        force_purge, test_root = parse_command_line(sys.argv, description)

    sys.exit(load_balancing_submit(compset, res, pesfile, mpilib,
                                   compiler, project, machine,
                                   extra_options_file, casename_prefix,
                                   force_purge, test_root))

###############################################################################

if __name__ == '__main__':
    _main_func(__doc__)
