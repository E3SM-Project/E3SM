To enable multiprocessing rather than running in serial, the program will need to be run in an
**interactive session** on compute nodes, or as a **batch job**.

Here are some hardware details for `Chrysalis`:
   * 64 cores/node
   * 256 GB DRAM/node
   * 32768 total cores


Interactive session on compute nodes
'''''''''''''''''''''''''''''''''''''

First, request an interactive session with a single node
for one hour (running this example should take much less than this).

    ::

        srun --pty --nodes=1 --time=01:00:00 /bin/bash



Once the session is available, launch E3SM Diagnostics, to activate ``e3sm_unified``:

    ::

        source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
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
        #SBATCH --account=e3sm
        #SBATCH --nodes=1
        #SBATCH --time=01:00:00

        source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
        python run_e3sm_diags.py --multiprocessing --num_workers=32

And then submit it:

    ::

        sbatch diags.bash
