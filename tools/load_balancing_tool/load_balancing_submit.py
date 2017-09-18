#!/usr/bin/env python
"""
Script to submit a series of ACME runs to get data for
time vs nprocessors model. This data will be used to generate
a processor layout that achieves high efficiency
"""
from xml.etree.ElementTree import ParseError
import argparse
import shutil

try:
    from Tools.standard_script_setup import *
except ImportError, e:
    print 'Error importing Tools.standard_script_setup'
    print 'May need to add cime/scripts to PYTHONPATH\n'
    raise ImportError(e)

from CIME.utils import run_cmd_no_fail
from CIME.XML import pes

logger = logging.getLogger(__name__)

DEFAULT_CASENAME_PREFIX = 'lbt_timing_run_'

# Default CIME variables, these can be overridden using the
# --extra_options_file option
CIME_DEFAULTS = {
    'STOP_OPTION':'ndays',
    'STOP_N':'10',
    'REST_OPTION':'never',
    'DOUT_S':'FALSE',
    'COMP_RUN_BARRIERS':'TRUE',
    'TIMER_LEVEL':'9'
}

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
      <pes compset="any" pesize="1">
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
      <pes compset="any" pesize="2">
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
      <pes compset="any" pesize="3">
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

    parser.add_argument('--extra_options_file',
                        help='file listing options to be run using xmlchange')
    parser.add_argument('--casename_prefix', default=DEFAULT_CASENAME_PREFIX,
                        help='casename prefix to use for all timing runs')
    parser.add_argument('--force_purge', action='store_true')

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args, parser)

    return (args.compset, args.res, args.pesfile,
            args.compiler, args.project, args.machine, args.extra_options_file,
            args.casename_prefix, args.force_purge)

def _set_xml_val(option, value, casedir):
    """
    Call xmlchange from command line to set option=value
    """
    cmd = './xmlchange %s=%s' % (option, value)
    logger.info(cmd)
    return run_cmd_no_fail(cmd, from_dir=casedir)

################################################################################
def load_balancing_submit(compset, res, pesfile, compiler, project, machine,
                          extra_options_file, casename_prefix, force_purge):
################################################################################
    create_newcase_flags = ' --compset %s --res %s --handle-preexisting-dirs r ' % (compset, res)
    if machine is not None:
        create_newcase_flags += ' --machine ' + machine
    if project is not None:
        create_newcase_flags += ' --project ' + project
    if compiler is not None:
        create_newcase_flags += ' --compiler ' + compiler

    # Read in list of pes from given file
    if not os.access(pesfile, os.R_OK):
        logger.critical('ERROR: File %s not found', pesfile)
        raise SystemExit(1)
    logger.info('Reading XML file %s. Searching for pesize entries:', pesfile)
    try:
        pesobj = CIME.XML.pes.Pes(pesfile)
        logger.info(str(pesobj))
    except ParseError:
        logger.critical('ERROR: File %s not parseable', pesfile)
        raise SystemExit(1)

    pesize_list = []
    for node in pesobj.get_nodes('pes'):
        pesize = node.get('pesize')
        if not pesize:
            logger.critical('No pesize for pes node in file %s', pesfile)
        if pesize in pesize_list:
            logger.critical('pesize %s duplicated in file %s', pesize, pesfile)
        logger.info('  ' + pesize)
        pesize_list.append(pesize)

    if not pesize_list:
        logger.critical('ERROR: No grid entries found in pes file %s', pesfile)
        raise SystemExit(1)

    # Submit job for each entry in pesfile
    script_dir = CIME.utils.get_scripts_root()
    for pesize in pesize_list:
        casename = casename_prefix + pesize
        casedir = os.path.join(script_dir, casename)
        logger.info('Case %s...', casename)
        if os.path.exists(casedir):
            if force_purge:
                logger.info('Removing directory %s', casedir)
                shutil.rmtree(casedir)
            elif os.path.exists(os.path.join(casedir, 'timing')):
                logger.info('Skipping timing run for %s, because case '
                            'already exists.\n   To force rerun of all tests, '
                            'use --force_purge option.\n   To force rerun of '
                            'this test only, remove directory \n      %s',
                            casename, casedir)
                continue
            else:
                logger.critical("ERROR: directory %s exists,\n  but no timing "
                                "information available. Either run with "
                                "--force_purge\n  (remove all timing runs) or "
                                "remove this directory and try again", casedir)
                sys.exit(1)

        cmd = '%s/create_newcase --case %s --output-root %s ' % \
              (script_dir, casename, casedir)
        cmd += create_newcase_flags
        logger.info(cmd)
        run_cmd_no_fail(cmd, from_dir=script_dir)
        for var in CIME_DEFAULTS:
            _set_xml_val(var, CIME_DEFAULTS[var], casedir)

        if extra_options_file is not None:
            try:
                extras = open(extra_options_file, 'r')
                for line in extras.readlines():
                    split = line.split('=')
                    if len(split) == 2:
                        logger.info('setting %s=%s', split[0], split[1])
                        _set_xml_val(split[0], split[1], casedir)
                    else:
                        logger.debug('ignoring line in %s: %s',
                                     extra_options_file, line)
                extras.close()
            except IOError, e:
                logger.critical("ERROR: Could not read file %s",
                                extra_options_file)
                raise SystemExit(1)

        # There should be a better way to do this
        pes_ntasks, pes_nthrds, pes_rootpe, ignore = \
            pesobj.find_pes_layout('any', 'any', 'any', pesize_opts=pesize)

        for key in pes_ntasks:
            _set_xml_val(key, pes_ntasks[key], casedir)
        for key in pes_nthrds:
            _set_xml_val(key, pes_nthrds[key], casedir)
        for key in pes_rootpe:
            _set_xml_val(key, pes_rootpe[key], casedir)

        cmd = './case.setup'
        logger.info(cmd)
        run_cmd_no_fail(cmd, from_dir=casedir)

        cmd = './case.build'
        logger.info(cmd)
        run_cmd_no_fail(cmd, from_dir=casedir)

        cmd = './case.submit'
        logger.info(cmd)
        run_cmd_no_fail(cmd, from_dir=casedir)

    logger.info('Timing jobs submitted. After jobs completed, run to optimize '
                'pe layout:\n  load_balancing_solve --casename_prefix %s',
                (casename_prefix))

###############################################################################
def _main_func(description):
###############################################################################
    compset, res, pesfile, compiler, project, machine, extra_options_file, casename_prefix, force_purge = parse_command_line(sys.argv, description)

    sys.exit(load_balancing_submit(compset, res, pesfile,
                                   compiler, project, machine,
                                   extra_options_file, casename_prefix,
                                   force_purge))

###############################################################################

if __name__ == '__main__':
    _main_func(__doc__)
