Creating a New Environment from Scratch
=======================================

Every once in a while (fairly often), the yml environmental files
`in the conda folder <https://github.com/E3SM-Project/e3sm_diags/tree/master/conda>`__
break. It's either because:

- A dependency is removed.
- A dependency of a dependency is removed.
- ???

When this happens, you need to create a new environment from scratch.
Here's how you do it.

1. On a Linux machine, use the command below to create an environment: ::

    conda create -n e3sm_diags_env e3sm_diags python=3 mesalib -c e3sm -c conda-forge -c cdat

It gets the latest Python 3 version of ``e3sm_diags``.


2. Activate your new environment and ``cd`` to the ``tests/system/`` folder.


3. Run the command below, checking that the output is correct. ::

    python all_sets.py -d all_sets.cfg

**If you're creating an env for a version below v2.0.0, use the command below:** ::

    e3sm_diags -p all_sets.py -d all_sets.cfg

4. If you get an error while running, contact the CDAT team.
It might be a CDAT dependency that's causing the issue.
They might make you use Anaconda to force install the correct versions of things.


5. Once you have a working, environment, dump it into a file called ``e3sm_diags_env.yml``. ::

    conda env export > e3sm_diags_env.yml


6. Edit the file and remove most of the build numbers.

  * The build numbers are the numbers after the last '=' sign.
    Ex: It's ``py36h0aa2c8f_1004`` in ``cartopy=0.17.0=py36h0aa2c8f_1004``.
    So the line would be ``cartopy=0.17.0``, without the final '=' sign.
  * **If the build number doesn't have a hash and seems to have some build-specific**
    **feature, do not remove the build number.**
    In the example below, we will **not remove the build numbers**, because they look important.

    * ``openblas`` in ``blas=1.1=openblas``
    * ``nompi`` in ``hdf5=1.10.5=nompi_h3c11f04_1100``
    * ``mpich`` in ``mpi=1.0=mpich``
    * ``blas_openblas`` in ``numpy=1.16.2=py36_blas_openblash1522bff_0``
    * ``blas_openblas`` in ``py36_blas_openblash1522bff_0``
    * ``mesalib`` in ``py36_mesalibh24c825c_0``


7. Also, in the ``yml`` file, update any versions of any dependencies you need.
Ex: Updating ``cdp=1.4.2`` to ``cdp=1.6.0``.


8. To be thorough, create an environment with this new file again to
see if everything is working after all of the changes. ::

    conda env create -f e3sm_diags_env.yml --name e3sm_diags_env_test


9. Run steps 2 and 3 again.


10. Then, ``cp`` this file into a ``e3sm_diags_env_dev.yml``.
Open this new file and do the following:

    * Change ``name: e3sm_diags_env`` to ``name: e3sm_diags_env_dev``.
    * Under channels, remove: ::

        - e3sm
    * Under ``dependencies``, remove: ::

      - e3sm_diags=1.7.0
    * Change ``prefix: /global/homes/z/zshaheen/anaconda2/envs/e3sm_diags_env`` 
      to ``prefix: /global/homes/z/zshaheen/anaconda2/envs/e3sm_diags_env_dev``.


11. Commit and push both ``e3sm_diags_env.yml`` and ``e3sm_diags_env_dev.yml``.
