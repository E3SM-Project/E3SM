To enable multiprocessing rather than running in serial, the program will need to be run in an
**interactive session** on compute nodes, or as a **batch job**. 


Interactive session on compute nodes
'''''''''''''''''''''''''''''''''''''

First, request an interactive session with a single node
(32 cores with Cori Haswell, 68 cores with Cori KNL)
for one hour (running this example should take much less than this).
If obtaining a session takes too long, try to use the ``debug`` partition.
Note that the maximum time allowed for that partition is ``00:30:00``.

    ::

        salloc --nodes=1 --partition=regular --time=01:00:00 -C haswell


Once the session is available, launch E3SM Diagnostics, to activate ``e3sm_unified``:

    ::

        source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh
        python run_e3sm_diags.py --multiprocessing --num_workers=32


We could have also set these multiprocessing parameters in the ``run_e3sm_diags.py`` as well
but we're showing that you can still submit parameters via the command line.

Batch job
'''''''''

Alternatively, you can also create a script and submit it to the batch system.
Copy and paste the code below into a file named ``diags.bash``.

    .. code:: bash
    
        #!/bin/bash -l
        #SBATCH --job-name=diags
        #SBATCH --output=diags.o%j
        #SBATCH --partition=regular
        #SBATCH --account=acme
        #SBATCH --nodes=1
        #SBATCH --time=01:00:00
        #SBATCH -C haswell

        source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh
        python run_e3sm_diags.py --multiprocessing --num_workers=32

And then submit it:

    ::

        sbatch diags.bash

View the status of your job with ``squeue -u <username>``.
Here's the meaning of some values under the State (``ST``) column:

* ``PD``: Pending
* ``R``: Running
* ``CA``: Cancelled
* ``CD``: Completed
* ``F``: Failed
* ``TO``: Timeout
* ``NF``: Node Failure
