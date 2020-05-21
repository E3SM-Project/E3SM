# Directory Structure

* `fortran`: Original fortran SAM (will not change)
* `cpp`: Progressively ported C++ version of SAM
* `test`: Testing framework

# Simplifying SAM

We are only porting the 1-moment microphysics, turbulence, and dynamical core to C++. Therefore, we stripped out CLUBB, M2005 microphysics, MAML, ESMT, and ECPP. In the end, we only had roughly 8K lines of code left to port. The resulting code is in `SCAMPAM/dynamics/sam++/fortran`.

# Testing Approach

We use a single monolithic standalone test driver to avoid the overhead of writing unit tests for each function since the porting effort needs to be fairly speedy. This test exercises all of the code in the model using sample CRM instances output from a realistic model run. The tests are located in `dynamics/sam++/test`

We create a baseline `sam++/fortran` directory that is run every time as well as a `sam++/cpp` directory containing the progressively ported hybrid C++ - Fortran code. The `fortran` directory's code will stay the same to serve as a reliable baseline for testing to make sure we are getting the same answer as we port portions of the code to C++.

# Making the Fortran Code C Compatible

In order to easily enable `bind(C)` on each of the routines for easier testing, we changed all `integer` types to `integer(crm_iknd)`, and we changed all `logical` types to `logical(crm_lknd)`. This was done with a `sed` script, and we then had to set `crm_rknd = c_double`, `crm_iknd = c_int`, and `crm_lknd = c_bool`. Note that `c_bool` does not work with the PGI compiler, but we're using GNU, and the goal is not to run the Fortran code outright but to use it for test-driven development for the C++ code.

Next, we placed `bind(C)` on all of the scalar module-level variables, and we mirrored the Fortran `parameter`s in those files with C++ `constexpr` scalars in `const.h`. next, we created a list of `extern` variables in `vars.h` for C++ code to access the module-level scalars declared with `bind(C)` in the Fortran code. Finally, we declare all module-level arrays as unmanaged Kokkos `View`s with `LayoutRight` memory (fully contiguous with the right-most index varying the fastest) in the global C++ namespace. Then, we pass all of the module-level arrays to a `wrap_arrays(...)` C++ function that initializes the Kokkos Views with pointers from Fortran. In the `iso_c_binding` interface, each array is accepted as a single-dimensional assumed-size array for ease: e.g.,

```fortran
real(crm_rknd), dimension(*) :: arr1, arr2, ...
```

Once these Fortran pointers are wrapped by Kokkos Views, we can access them in the global namespace, allowing us to keep the code similar to the existing Fortran code, which does not pass much data by parameter but rather accesses it via `use` statements from modules.

By keeping the Fortran interfaces to C++ routines simpler with already created View-wrapped Fortran arrays in the C++ code, we are able to do the C++ porting work in a fast and readable manner.

# C++ Porting Workflow

Regarding actually porting routines to C++, since we've placed pre-wrapped Views of all Fortran module data in the C++ global namespace already, the parameter lists do not need to change. Any parameters that must be passed follow the interface given in the next section's table. Each source file must have its own corresponding header file. In the `sam++/cpp` directory, you can use `abcoefs.*`, `adams.*`, `vars.*`, and `const.h` to aid you in what porting a routine will look like. The following are things to make sure you pay attention to:

* Since we're sharing indirect indexing integers with Fortran, you will need to subtract 1 from them in the C++ code
  * E.g., `dudt(nc-1,k,j,i,icrm)`
* You must transpose the dimension in all multi-dimensional arrays
  * `rho(icrm,k)`  -->  `rho(k,icrm)`
* You'll need to use offsets for any Fortran array that has non-1 lower bounds, and these are declared in `const.h`, there's an example of their use in `adams.cpp`, and all arrays that need offsets (and which dimensions need offsets) are given at the **top** of `vars.h`.
  * E.g., `u(k,j+offy_u,i+offx_u,icrm)`
* For loops should loop over values that are 1 less than the Fortran loops
  * E.g., `do i = 1 , nx`  -->  `for (int i=0; i<nx; i++) {`
* Don't forget to put `extern "C"` before all C++ function definitions. This tells the C++ compiler to give it a simple name in the object file that is compatible with the `iso_c_binding`
* Don't forget to add the interface to `cpp_interface_mod.F90`, to `use` that module in the context you're calling the routine, and to remove the Fortran reference to the routine.
  * E.g. `use adams_mod, only: adams` --> `use cpp_interface_mod`

## An Example Conversion

The `adams.F90` file is below:

```fortran
module adams_mod
  use params, only: asyncid
  implicit none
contains
  subroutine adams(ncrms)
    !       Adams-Bashforth scheme
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) dtdx, dtdy, dtdz, rhox, rhoy, rhoz
    integer i,j,k,icrm

    dtdx = dtn/dx
    dtdy = dtn/dy

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dtdz = dtn/dz(icrm)
            rhox = rho(icrm,k)*dtdx
            rhoy = rho(icrm,k)*dtdy
            rhoz = rhow(icrm,k)*dtdz
            dudt(icrm,i,j,k,nc) = u(icrm,i,j,k) + dt3(na) *(at*dudt(icrm,i,j,k,na)+bt*dudt(icrm,i,j,k,nb)+ct*dudt(icrm,i,j,k,nc))
            dvdt(icrm,i,j,k,nc) = v(icrm,i,j,k) + dt3(na) *(at*dvdt(icrm,i,j,k,na)+bt*dvdt(icrm,i,j,k,nb)+ct*dvdt(icrm,i,j,k,nc))
            dwdt(icrm,i,j,k,nc) = w(icrm,i,j,k) + dt3(na) *(at*dwdt(icrm,i,j,k,na)+bt*dwdt(icrm,i,j,k,nb)+ct*dwdt(icrm,i,j,k,nc))
            u(icrm,i,j,k) = 0.5*(u(icrm,i,j,k)+dudt(icrm,i,j,k,nc)) * rhox
            v(icrm,i,j,k) = 0.5*(v(icrm,i,j,k)+dvdt(icrm,i,j,k,nc)) * rhoy
            misc(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc))
            w(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc)) * rhoz
          end do
        end do
      end do
    end do
  end subroutine adams
end module adams_mod
```

This was ported to C++, and this is the resulting code:

```C++
#include "adams.h"

extern "C" void adams() {
  // Adams-Bashforth scheme
  real dtdx = dtn/dx;
  real dtdy = dtn/dy;

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  Kokkos::parallel_for( "adams" , nzm*ny*nx*ncrms , KOKKOS_LAMBDA (int iGlob) {
    int k, j, i, icrm;
    yakl::unpackIndices( iGlob , nzm , ny , nx , ncrms , k , j , i , icrm );
    real dtdz = dtn/dz(icrm);
    real rhox = rho (k,icrm)*dtdx;
    real rhoy = rho (k,icrm)*dtdy;
    real rhoz = rhow(k,icrm)*dtdz;
    real utend = ( at*dudt(na-1,k,j,i,icrm) + bt*dudt(nb-1,k,j,i,icrm) + ct*dudt(nc-1,k,j,i,icrm) );
    real vtend = ( at*dvdt(na-1,k,j,i,icrm) + bt*dvdt(nb-1,k,j,i,icrm) + ct*dvdt(nc-1,k,j,i,icrm) );
    real wtend = ( at*dwdt(na-1,k,j,i,icrm) + bt*dwdt(nb-1,k,j,i,icrm) + ct*dwdt(nc-1,k,j,i,icrm) );
    dudt(nc-1,k,j,i,icrm) = u(k,j+offy_u,i+offx_u,icrm) + dt3(na-1) * utend;
    dvdt(nc-1,k,j,i,icrm) = v(k,j+offy_v,i+offx_v,icrm) + dt3(na-1) * vtend;
    dwdt(nc-1,k,j,i,icrm) = w(k,j+offy_w,i+offx_w,icrm) + dt3(na-1) * wtend;
    u   (k,j+offy_u,i+offx_u,icrm) = 0.5 * ( u(k,j+offy_u,i+offx_u,icrm) + dudt(nc-1,k,j,i,icrm) ) * rhox;
    v   (k,j+offy_v,i+offx_v,icrm) = 0.5 * ( v(k,j+offy_v,i+offx_v,icrm) + dvdt(nc-1,k,j,i,icrm) ) * rhoy;
    w   (k,j+offy_w,i+offx_w,icrm) = 0.5 * ( w(k,j+offy_w,i+offx_w,icrm) + dwdt(nc-1,k,j,i,icrm) ) * rhoz;
    misc(k,j       ,i       ,icrm) = 0.5 * ( w(k,j+offy_w,i+offx_w,icrm) + dwdt(nc-1,k,j,i,icrm) );
  }); 
}
```

# Turn on Kokkos Debugging

You'll find Kokkos debugging to be extremely helpful during development. It causes the code to run slower, but it tells you the exact line number where an illegal index occurs. To turn it on, set the environment variable `KOKKOS_DEBUG="yes"`.

# Mapping Fortran and C++ types

Types map in the following manner between Fortran and C++

| C++            | Fortran Interface               | Fortran Description                             |
|----------------|---------------------------------|-------------------------------------------------|
| `int    var`   | `integer(c_int) , value`        | `intent(in)` scalars must be passed by value    |
| `double var`   | `real(c_double) , value`        | ''                                              |
| `float  var`   | `real(c_float)  , value`        | ''                                              |
| `bool   var`   | `logical(c_bool), value`        | ''                                              |
| `int    &var`  | `integer(c_int)`                | `intent(out)` scalars whose values will change  |
| `float  &var`  | `real(c_float)`                 | ''                                              |
| `double &var`  | `real(c_double)`                | ''                                              |
| `bool   &var`  | `integer(c_bool)`               | ''                                              |
| `int    *var`  | `integer(c_int) , dimension(*)` | Valid for passing arrays of **any dimension**   |
| `float  *var`  | `real(c_float)  , dimension(*)` | ''                                              |
| `double *var`  | `real(c_double) , dimension(*)` | ''                                              |
| `bool   *var`  | `integer(c_bool), dimension(*)` | ''                                              |

# Running the Regression Test

Before running the tests, you'll need to load the git submodules. To do this, go to the root directory of the repo, and type `git submodule update --init`

To facilitate a fast C++ port, we chose to drive a standalone framework to test the answers in the model. Thus far, everything has been bit-for-bit the same while we're porting the work serially. This will change once we run the code in parallel, though. To run the regression test, you'll first need to set some environment variables.

## Test Environment Variables

**Mandatory Environment Variables**:

```bash
# If you put all netcdf in the same place, then specify the same path for both
# We need libnetcdf.a, libnetcdff.a, all include files, nf-config, and nc-config
NCHOME=/path/to/netcdf-c
NFHOME=/path/to/netcdf-fortran
```

**Optional Environment Variables**:

```bash
NCRMS=30     # Sets the number of CRMs to drive the test with
             # The lower this number, the faster the tests will run
FFLAGS=-O3   # Gives optimization Fortran flags
CXXFLAGS=-O3 # Gives optimization C++ flags
KOKKOS_DEBUG="yes" # Setting this will turn on Kokkos debugging
```

## Test Workflow

To run the regression tests, after you set your environment variables, the following is the workflow:

```bash
cd dynamics/sam++/test/build
./download_data.sh  #Only needed once for all time
./cmakescript.sh  2dinputfile.nc  3dinputfile.nc
make -j
./runtests.sh
```

When you add new source files (don't forget to `git add` them!), you'll need to reconfigure with `./cmakescript.sh`.

The role of `cmakescript.sh` is to grab information from the 2-D and 3-D CRM NetCDF files that drive the tests and to ensure the models are compiled with settings that are compatible with the driving NetCDF files. After your run make, you will have four executables created: `fortran2d`, `fortran3d`, `cpp2d`, and `cpp3d`, each of them in their own subdirectory. 

fortran - the original fortran code

cpp - progressively ported C++ version of above

test - testing framework
