#!/usr/bin/env python

"""
Generates a cylc workflow file for the case.  See https://cylc.github.io for details about cylc
"""

from standard_script_setup import *

from CIME.case              import Case
from CIME.utils             import expect, transform_vars

import argparse, re
logger = logging.getLogger(__name__)

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter)

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument("caseroot", nargs="?", default=os.getcwd(),
                        help="Case directory for which namelists are generated.\n"
                        "Default is current directory.")

    parser.add_argument('--cycles', default=1,
                        help="The number of cycles to run, default is RESUBMIT")

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args, parser)

    return args.caseroot, args.cycles

def cylc_batch_job_template(job, caseroot, batch_system_type):
    return """
  [[{job}]]
    script = cd {caseroot}; ./case.submit --job {job}
    [[[job]]]
      batch_system = {batch_system_type}
    [[[directives]]]
""".format(job=job, caseroot=caseroot, batch_system_type=batch_system_type) + "{{ batchdirectives }}\n"


def cylc_script_job_template(job, caseroot):
    return """
  [[{job}]]
    script = cd {caseroot}; ./case.submit --job {job}
""".format(job=job, caseroot=caseroot)

###############################################################################
def _main_func(description):
###############################################################################
    caseroot, cycles = parse_command_line(sys.argv, description)
    
    expect(os.path.isfile(os.path.join(caseroot, "CaseStatus")),
           "case.setup must be run prior to running {}".format(__file__))
    with Case(caseroot, read_only=True) as case:
        if cycles == 1:
            cycles = max(1, case.get_value('RESUBMIT'))
        env_batch = case.get_env('batch')
        env_workflow = case.get_env('workflow')
        jobs = env_workflow.get_jobs()
        casename = case.get_value('CASE')
        input_template = os.path.join(case.get_value("MACHDIR"),"cylc_suite.rc.template")
        batch_system_type = env_batch.get_batch_system_type()
        overrides = {"cycles":cycles, 
                     'casename':casename}
        input_text = open(input_template).read()
        for job in jobs:
            jobname = job
            if job == 'case.st_archive':
                continue
            if job == 'case.run':
                jobname = 'run'
                overrides.update(env_batch.get_job_overrides(job, case))
                overrides.update({'job_id':'run.'+casename})
                input_text = input_text + cylc_batch_job_template(jobname, caseroot, batch_system_type)
            else:
                depends_on = env_workflow.get_value('dependency', subgroup=job)
                if depends_on.startswith('case.'):
                    depends_on = depends_on[5:]
                input_text = input_text.replace(' => '+depends_on,' => '+depends_on+' => '+job)
                    

                overrides.update(env_batch.get_job_overrides(job, case))
                overrides.update({'job_id':job+'.'+casename})
                if 'total_tasks' in overrides and overrides['total_tasks'] > 1:
                    input_text = input_text + cylc_batch_job_template(jobname, caseroot, batch_system_type)
                else:
                    input_text = input_text + cylc_script_job_template(jobname, caseroot)


            overrides.update({'batchdirectives':env_batch.get_batch_directives(case,job, overrides=overrides, 
                                                                           output_format='cylc')})


            input_text = transform_vars(input_text, case=case, subgroup=job, overrides=overrides)
        with open("cylc_suite.rc", "w") as f:
                    f.write(input_text)
                


if (__name__ == "__main__"):
    _main_func(__doc__)
