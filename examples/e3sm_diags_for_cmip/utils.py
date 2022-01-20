import os
import shlex
import time
from subprocess import Popen, PIPE

# -----------------------------------------------------------------------------
def submitScript(scriptFile, export='ALL'):

    command = "sbatch --export=%s %s" % (export,scriptFile)

    # Actual submission
    p1 = Popen(shlex.split(command),stdout=PIPE,stderr=PIPE)
    (stdout, stderr) = p1.communicate()
    status = p1.returncode
    out = stdout.decode().strip()
    print('...%s' % (out))
    if status != 0 or not out.startswith("Submitted batch job"):
        print('Problem submitting script %s' % (scriptFile))
        print(stderr)
        raise Exception
    jobid = int( out.split()[-1] )

    # Small pause to avoid overloading queueing system
    time.sleep(0.5)

    return jobid
