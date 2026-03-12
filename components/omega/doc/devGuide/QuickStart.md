(omega-dev-quick-start)=

# Quick Start for Developers

## Getting set up

### Creating an E3SM fork

We ask developers to make their own fork of the E3SM to use for development
branches.  To do this, go to the page for the
[E3SM-Project/E3SM](https://github.com/E3SM-Project/E3SM) repository and
click on the `fork` button near the upper right corner.  Set "owner" to your
GitHub username (this should be the default) and click "Create fork".

There is no need to have a separate fork for Omega, since the Omega repository
is also a fork of E3SM.  Your E3SM fork will serve for both E3SM and Omega
development.  In fact, GitHub will not let you make an Omega fork if you
already have an E3SM fork (or vice versa).

### Check out the code

Make sure you have added the SSH key for the machine you are using to your
GitHub account, see
[Adding a new SSH key to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).

Clone the develop branch of the Omega repo:
```sh
git clone git@github.com:E3SM-Project/Omega.git develop
```

You may wish to rename `origin` to `E3SM-Project/Omega` for clarity:
```sh
git remote rename origin E3SM-Project/Omega
```

Consider adding some other remotes, e.g.:
```sh
git remote add E3SM-Project/E3SM git@github.com:E3SM-Project/E3SM.git
git remote add <github_username>/E3SM git@github.com:<github_username>/E3SM.git
```
where `<github_username>` is your username or that of a collaborator.

To sync your local clone with all the remotes, run:
```sh
git fetch --all -p
```

### Create a conda environment

If you do not have conda set up in your own space yet, follow
{ref}`omega-dev-install-miniforge3` for instructions on how to install it.


Activate the `base` conda environment.  Go to the base of the Omega branch you
plan to develop and create a conda environment for
linting your code and building the documentation:
```sh
cd components/omega/
conda create -n omega_dev --file dev-conda.txt
conda activate omega_dev
```
(You can reuse `omega_dev` for other branches as long as `dev-conda.txt` has
not changed between branches.)

The first time you set up the environment, you need to initiate `pre-commit`:
```sh
pre-commit install
```
This needs to be done once per environment, so re-run it as necessary if you
re-create the development environment.

Please activate this environment each time you commit code.  This will ensure
that `pre-commit` and the associated linting utilities are available and that
they check your code as it is committed (rather than requiring fix-up commits
later on).

(omega-dev-quick-start-build-test)=

## Building and testing Omega

If you are using Polaris, you generally should not need to manually build
Omega. Instead, rely on Polaris auto-build during `polaris setup` or
`polaris suite`:

- Polaris quick start (developers):
    <https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html>
- Building E3SM components:
    <https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#building-e3sm-components>

For Omega development, you are expected to validate changes with the Polaris
`omega_pr` suite before non-trivial PRs (see {ref}`omega-dev-testing`):

- `omega_pr` suite:
    <https://docs.e3sm.org/polaris/main/users_guide/ocean/suites.html#omega-pr-suite>

Both Polaris-based and standalone/manual workflows are supported.

(omega-dev-quick-start-ctest-util)=

### Polaris CTest Utility

If you are using Polaris, you may wish to use its
[Omega CTest utility](https://github.com/E3SM-Project/polaris/tree/main/utils/omega/ctest)
to build and test Omega. The utility automates many of the steps below.

(omega-dev-quick-start-build)=

### Building Omega

If you are not using Polaris, follow the standalone/manual CMake workflow in
{ref}`omega-dev-cmake-build`.

At a high level:

1. Update required submodules.
2. Configure CMake in an out-of-source build directory.
3. Build with `./omega_build.sh`.
4. Run unit tests with `./omega_ctest.sh`.

(omega-dev-quick-start-getting-meshes)=
### Getting test meshes

Some tests require a valid Omega mesh file. Different tests require different
meshes. At the moment, mesh files need to be copied or linked to specifically
named files under the `test` directory. Appropriate mesh files can be
downloaded from:
- [Ocean Mesh](https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/oQU240/ocean.QU.240km.151209.nc)
- [Global Mesh](https://web.lcrc.anl.gov/public/e3sm/polaris/ocean/polaris_cache/global_convergence/icos/cosine_bell/Icos480/init/initial_state.230220.nc)
- [Planar Mesh](https://gist.github.com/mwarusz/f8caf260398dbe140d2102ec46a41268/raw/e3c29afbadc835797604369114321d93fd69886d/PlanarPeriodic48x48.nc)
```sh
wget -O ocean_test_mesh.nc https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/oQU240/ocean.QU.240km.151209.nc
wget -O global_test_mesh.nc https://web.lcrc.anl.gov/public/e3sm/polaris/ocean/polaris_cache/global_convergence/icos/cosine_bell/Icos480/init/initial_state.230220.nc
wget -O planar_test_mesh.nc https://gist.github.com/mwarusz/f8caf260398dbe140d2102ec46a41268/raw/e3c29afbadc835797604369114321d93fd69886d/PlanarPeriodic48x48.nc
cd test
ln -sf  ../ocean_test_mesh.nc OmegaMesh.nc
ln -sf  ../global_test_mesh.nc OmegaSphereMesh.nc
ln -sf  ../planar_test_mesh.nc OmegaPlanarMesh.nc
cd ..
```

### Running CTests

Omega includes several unit tests that run through CTest. The unit tests need
to be run on a compute node.

To run the tests:
```sh
./omega_ctest.sh
```

The results should look something like:
```
Test project /gpfs/fs1/home/ac.xylar/e3sm_work/polaris/add-omega-ctest-util/build_omega/build_chrysalis_intel
    Start 1: DATA_TYPES_TEST
1/9 Test #1: DATA_TYPES_TEST ..................   Passed    0.38 sec
    Start 2: MACHINE_ENV_TEST
2/9 Test #2: MACHINE_ENV_TEST .................   Passed    0.98 sec
    Start 3: BROADCAST_TEST
3/9 Test #3: BROADCAST_TEST ...................   Passed    1.13 sec
    Start 4: LOGGING_TEST
4/9 Test #4: LOGGING_TEST .....................   Passed    0.03 sec
    Start 5: DECOMP_TEST
5/9 Test #5: DECOMP_TEST ......................   Passed    1.20 sec
    Start 6: HALO_TEST
6/9 Test #6: HALO_TEST ........................   Passed    1.08 sec
    Start 7: IO_TEST
7/9 Test #7: IO_TEST ..........................   Passed    2.94 sec
    Start 8: CONFIG_TEST
8/9 Test #8: CONFIG_TEST ......................   Passed    1.01 sec
    Start 9: KOKKOS_TEST
9/9 Test #9: KOKKOS_TEST ........................   Passed    0.03 sec

100% tests passed, 0 tests failed out of 9

Total Test time (real) =   8.91 sec
```

### Debugging tips

If Omega CTests are failing or simulations are crashing, setting
`OMEGA_BUILD_TYPE` to `Debug` can be helpful for debugging purposes. If you
need to identify which test has failed, it may be useful to examine the CTest
log file located at `$BUILD_DIR/Testing/Temporary/LastTest.log`.

(omega-dev-parmetis-libs)=

### Metis and Parmetis libraries

The following table shows locations for Metis and Parmetis libraries on
supported E3SM machines. The pattern is:
```
<polaris_base>/<machine>/spack/dev_polaris_0_10_0_<compiler>_<mpi>/var/spack/environments/dev_polaris_0_10_0_<compiler>_<mpi>/.spack-env/view
```

```{eval-rst}
+--------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Machine      | Compiler         | Parmetis path                                                                                                                                                               |
+==============+==================+=============================================================================================================================================================================+
| aurora       | oneapi-ifx       | /lus/flare/projects/E3SM_Dec/soft/polaris/aurora/spack/dev_polaris_0_10_0_oneapi-ifx_mpich/var/spack/environments/dev_polaris_0_10_0_oneapi-ifx_mpich/.spack-env/view       |
+--------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| chicoma-cpu  | gnu              | /usr/projects/e3sm/polaris/chicoma-cpu/spack/dev_polaris_0_10_0_gnu_mpich/var/spack/environments/dev_polaris_0_10_0_gnu_mpich/.spack-env/view                               |
+--------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| chrysalis    | intel            | /lcrc/soft/climate/polaris/chrysalis/spack/dev_polaris_0_10_0_intel_openmpi/var/spack/environments/dev_polaris_0_10_0_intel_openmpi/.spack-env/view                         |
|              +------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|              | gnu              | /lcrc/soft/climate/polaris/chrysalis/spack/dev_polaris_0_10_0_gnu_openmpi/var/spack/environments/dev_polaris_0_10_0_gnu_openmpi/.spack-env/view                             |
+--------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| frontier     | craygnu          | /ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_10_0_craygnu_mpich/var/spack/environments/dev_polaris_0_10_0_craygnu_mpich/.spack-env/view                   |
|              +------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|              | craygnu-mphipcc  | /ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_10_0_craygnu-mphipcc_mpich/var/spack/environments/dev_polaris_0_10_0_craygnu-mphipcc_mpich/.spack-env/view   |
|              +------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|              | craycray         | /ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_10_0_craycray_mpich/var/spack/environments/dev_polaris_0_10_0_craycray_mpich/.spack-env/view                 |
|              +------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|              | craycray-mphipcc | /ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_10_0_craycray-mphipcc_mpich/var/spack/environments/dev_polaris_0_10_0_craycray-mphipcc_mpich/.spack-env/view |
+              +------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|              | crayamd          | /ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_10_0_crayamd_mpich/var/spack/environments/dev_polaris_0_10_0_crayamd_mpich/.spack-env/view                   |
|              +------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|              | crayamd-mphipcc  | /ccs/proj/cli115/software/polaris/frontier/spack/dev_polaris_0_10_0_crayamd-mphipcc_mpich/var/spack/environments/dev_polaris_0_10_0_crayamd-mphipcc_mpich/.spack-env/view   |
+--------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| pm-cpu       | gnu              | /global/cfs/cdirs/e3sm/software/polaris/pm-cpu/spack/dev_polaris_0_10_0_gnu_mpich/var/spack/environments/dev_polaris_0_10_0_gnu_mpich/.spack-env/view                       |
|              +------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|              | intel            | /global/cfs/cdirs/e3sm/software/polaris/pm-cpu/spack/dev_polaris_0_10_0_intel_mpich/var/spack/environments/dev_polaris_0_10_0_intel_mpich/.spack-env/view                   |
+--------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| pm-gpu       | gnugpu           | /global/cfs/cdirs/e3sm/software/polaris/pm-gpu/spack/dev_polaris_0_10_0_gnugpu_mpich/var/spack/environments/dev_polaris_0_10_0_gnugpu_mpich/.spack-env/view                 |
+--------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
```

## Code development

### Code style

We require that the C++ code in Omega adhere to the
[LLVM code style](https://llvm.org/docs/CodingStandards.html). These
conventions are enforced using the linting tools described in
{ref}`omega-dev-linting`.  This style may take new developers some time to get
used to but the hope is that it leads to a coherent code style in Omega.

### VS Code

You may wish to consider using an integrated development environment (IDE) to
develop your code.  A convenient option for developing on HPC is
[Visual Studio Code (VS Code)](https://code.visualstudio.com/).  It has plugins
for [C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools),
[CMake](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools),
[Python](https://marketplace.visualstudio.com/items?itemName=ms-python.python),
and so on.  If configured correctly, it should help to enforce the
[code formatting style](https://code.visualstudio.com/docs/cpp/cpp-ide#_code-formatting)
by recognizing Omega's `.clang-format` file.

A convenient feature is the ability to connect and edit code directly on HPC
systems from your laptop or desktop
[using an SSH connection](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh).

VS Code also provides a convenient way to preview the Markdown files used in
the Omega documentation as you are writing them.
