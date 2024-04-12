(omega-design-shallow-water-omega0)=
# Omega v0: Shallow Water

## 1 Overview

This design document describes the first version of the Omega ocean model, Omega0. Overall, Omega is an unstructured-mesh ocean model based on TRiSK numerical methods ([Thuburn et al. 2009](https://www.sciencedirect.com/science/article/pii/S0021999109004434)) that is specifically designed for modern exascale computing architectures. The algorithms in Omega will be nearly identical to those in MPAS-Ocean, but it will be written in c++ rather than Fortran in order to take advantage of libraries to run on GPUs, such as Kokkos (https://github.com/kokkos).

The planned versions of Omega are:

1. **Omega-0: Shallow water equations with identical vertical layers and inactive tracers.** There is no vertical transport or advection. The tracer equation is horizontal advection-diffusion, but tracers do not feed back to dynamics. Pressure gradient is simply gradient of sea surface height. Capability is similar to [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780)
1. **Omega-1: Minimal stand-alone primitive equation eddy-permitting model.** This adds active temperature, salinity, density, and pressure as a function of depth. There is a true pressure gradient. Vertical velocity is from the continuity equation. An equation of state and simple vertical mixing scheme are needed. Capability is same as [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760)
1. **Omega-2: Eddying ocean with advanced parameterizations.** This model will have sufficient capability to run realistic global simulations, similar to E3SM V1 [Petersen et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018MS001373).


We will produce separate design documents for the time-stepping scheme and the tracer equations in Omega-0. Some of the requirements stated here are repeated or implied by the Omega framework design documents, but are included for clarity.

## 2 Requirements

### 2.1 Requirement: Omega-0 will solve the nonlinear shallow water equations, plus inactive tracers

The governing equations for Omega-0 are conservation of momentum, volume, and tracers in a single layer:

$$
\frac{\partial \boldsymbol{u}}{\partial t} + q\left(h\boldsymbol{u}^{\perp}\right) = -g\nabla(h+b) - \nabla K + \nu_2 \nabla^2 \boldsymbol{u} - \nu_4 \nabla^4 \boldsymbol{u} + C_D \frac{\boldsymbol{u}\left|\boldsymbol{u}\right|}{h} - C_W \frac{(\boldsymbol{u}_W - \boldsymbol{u})\left|\boldsymbol{u}_W - \boldsymbol{u}\right|}{h}
\hspace{1cm}   (1)
$$

$$
\frac{\partial h}{\partial t} + \nabla \cdot \left(h \boldsymbol{u}\right) = 0,
\hspace{1cm}   (2)
$$

$$
\frac{\partial h \phi}{\partial t} + \nabla \cdot \left(h \boldsymbol{u} \phi\right) = \kappa_2 h \nabla^2 \phi - \kappa_4 h \nabla^4 \phi.
\hspace{1cm}   (3)
$$

The first two equations are from [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780), equations 2 and 7, with additional viscosity, bottom drag, and wind forcing (see equation 1 in [Lilly et al. 2023](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022MS003327)).
This equation set does not include any vertical advection or diffusion. Omega-0 will have a vertical index for performance testing and future expansion, but vertical layers will be simply redundant.
Additional information on governing equations may be found in chapter 8 of the MPAS User's Guide ([Petersen et al. 2018](https://zenodo.org/record/1246893)).

| symbol              | name                        | units    | location | notes                                                        |
|---------------------|-----------------------------|----------|----------|--------------------------------------------------------------|
| $C_D$               | bottom drag                 | 1/m      | edge     |                                                              |
| $C_W$               | wind stress coefficient     | 1/m      | edge     |                                                              |
| $f$                 | Coriolis parameter          | 1/s      | vertex   |                                                              |
| $g$                 | gravitational acceleration  | m/s^2    | constant |                                                              |
| $H$                 | total unperturbed depth     | m        | cell     |                                                              |
| $h$                 | thickness                   | m        | cell     |                                                              |
| ${\boldsymbol k}$   | vertical unit vector        | unitless | none     |                                                              |
| $K$                 | kinetic energy              | m^2/s^2  | cell     |  $K = \left\| {\boldsymbol u} \right\|^2 / 2$                |
| $q$                 | potential vorticity         | 1/m/s    | vertex   | $q = \eta/h$                                                 |
| $t$                 | time                        | s        | none     |                                                              |
| ${\boldsymbol u}$   | velocity                    | m/s      | edge     |                                                              |
| ${\boldsymbol u}_W$ | wind velocity               | m/s      | edge     |                                                              |
| $\eta$              | absolute vorticity          | 1/s      | vertex   | $\eta={\boldsymbol k} \cdot \nabla \times {\boldsymbol u}+f$ |
| $\kappa_2$          | tracer diffusion            | m^2/2    | constant |                                                              |
| $\kappa_4$          | biharmonic tracer diffusion | m^4/2    | constant |                                                              |
| $\nu_2$             | viscosity                   | m^2/2    | constant |                                                              |
| $\nu_4$             | biharmonic viscosity        | m^4/2    | constant |                                                              |
| $\phi$              | tracer                      | varies   | cell     |                                                              |

Note: Table created with [markdown table generator](https://www.tablesgenerator.com/markdown_tables) and original [google sheet](https://docs.google.com/spreadsheets/d/1rz-QXDiwfemq5NpSR1XsvomI7aSKQ1myTNweCY4afcE/edit#gid=0).

Boundary conditions will include both no-slip and free-slip. The original MPAS-Ocean implementation had no-slip boundaries, implemented by setting vorticity to zero at the boundaries. The free-slip implementation is described by Darren Engwirda in a [github discussion](https://github.com/E3SM-Ocean-Discussion/E3SM/pull/49).

### 2.2 Requirement: Numerical method will be the TRiSK formulation

The horizontal discretization will be taken from
[Thuburn et al. 2009](https://www.sciencedirect.com/science/article/pii/S0021999109004434) and [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780), as described in the algorithmic formulation in Section 3 below. This is the same base formulation as in MPAS-Ocean. In addition, we will consider small alterations from the original MPAS-Ocean horizontal discretization and include them as options if they are of minimal additional effort. This includes the recent AUST formulation in [Calandrini et al. 2021](https://www.sciencedirect.com/science/article/pii/S146350032100161X) and simple vorticity averaging considered by the Omega team this past year.

### 2.3 Requirement: Omega-0 will use MPAS format unstructured-mesh domains

We will continue to use the MPAS-format netcdf files for input and output, with the same mesh variable names and dimensions. This will facilitate ease of use and interoperability between MPAS-Ocean and Omega.

Omega mesh information will be stored in a separate mesh file, and not be included with initial condition, restart, or output files. This way there is never redundant mesh data stored in these files.

To do: In the design section below, add list of minimum number of variables that defines the mesh, and the additional variables that can be computed on start-up.

### 2.4 Requirement: Omega-0 will interface with polaris for preprocessing and postprocessing

A substantial investment has already been made in the polaris tools. Continuing the use of polaris for Omega will speed development and encourage the documentation and long-term reproducibility of test cases.

The test cases relevant to this design document are in Section 5 below.

### 2.5 Requirement: Omega-0 will run portably on various DOE architectures (CPU and GPU nodes)

Omega will be able to run on all the upcoming DOE architectures and make good use of GPU hardware. This should occur with minimal alterations in the high-level PDE solver code for different platforms.

Options include: writing kernels directly for GPUs in CUDA; adding OpenACC pragmas for the GPUs; or calling libraries such as Kokkos ([Trott et al. 2022](https://ieeexplore.ieee.org/document/9485033)), Kokkos (https://github.com/kokkos) or [HIP](https://github.com/ROCm-Developer-Tools/HIP) that execute code optimized for specialized architectures on the back-end, while providing a simpler front-end interface for the domain scientist.

### 2.6 Requirement: Omega-0 will run on multi-node with domain decomposition

This is with MPI and halo communication, as described in framework design documents. Results must be bit-for-bit identical across different number of partitions. This may be demonstrated with the ''QU240 partition test'' in Polaris.

### 2.7 Requirement: Correct convergence rates of operators and exact solutions

Each operator will be tested individually for the proper convergence rate: divergence, gradient, curl interpolated to cell centers, and tangential velocity are all second order; curl at vertices is first order. The details of the convergence tests are explained in section 5.1 below.

### 2.8 Requirement: Conservation of volume and tracer

The total volume of the domain should be conserved to machine precision in the absence of surface volume fluxes. Likewise, the total tracer amount should be conserved to machine precision without surface fluxes. A simple test is to initialize an inactive test tracer with a uniform value of 1.0, and it should remain 1.0 throughout the domain.

### 2.9 Requirement: Performance will be at least as good or better than MPAS-Ocean

Performance will be assessed with the following metrics:

1. Single CPU throughput
1. Parallel CPU scalability to high node counts
1. Single GPU throughput

The first two items will be compared to MPAS-Ocean with the same test configuration. For Omega-0, these may be tested first with the layered shallow water equations (momentum and thickness only) and be compared directly to the results in [Bishnu et al. 2023a](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-57), using the inertia-gravity wave test case described in Section 6.2 below. After tracer transport capabilities are added, performance may be compared against MPAS-Ocean with test cases using active tracers in addition to momentum and thickness.

For these comparisons, MPAS-Ocean will be set up with the identical terms and functionality as Omega-0. This means that vertical advection, vertical mixing, and all parameterizations will be disabled. Comparisons with these terms will be made with Omega-1 and higher and will be described in future design documents.

### 2.10 Requirement: Full-node GPU throughput will be comparable or better than full-node CPU throughput

For GPU throughput, comparisons should be made between full-node CPU throughput and full-node GPU throughput. For example, [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/architecture/) has nodes with 64 and 128 CPU cores (AMD EPYC 7763) and 4 GPUs (NVIDIA A100). We expect the full-node GPU throughput to be at least as good as the full-node CPU throughput, and potentially a factor of four higher. These numbers depend on the performance specifications of the particular hardware.

Like the previous requirement, tests will first be conducted with layered shallow water equations and later with additional tracer advection. For reference, [Bishnu et al. 2023a](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-57) was able to obtain nearly identical throughput between 64 CPU cores and a single GPU on Perlmutter using a Julia code and a layered shallow water test case.

## 3 Algorithmic Formulation

**To do in next pull request:** Copy over algorithmic formulations from [Thuburn et al. 2009](https://www.sciencedirect.com/science/article/pii/S0021999109004434) and [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780)

## 4 Design

To do in later PR: Add all design details.

1. Index ordering will be: (time, tracer, horizontal index, vertical index) Note this is c-style indexing, with last index being fastest, i.e. contiguous in memory.
2. Time stepping will be multi-level, with an arbitrary number of levels for multistage methods.

You can include code blocks like this:

```c++
int var = value;
```

### 4.1 Data types and parameters

#### 4.1.1 Parameters

to do: List and define any configuration parameters or public constants.

#### 4.1.2 Class/structs/data types

to do: Describe any public data types and/or the class definition

### 4.2 Methods

to do: List and describe all public methods and their interfaces (actual code for
interface that would be in header file). Describe typical use cases.

## 5 Verification and Testing

### 5.1 Convergence of individual terms

The following terms can be tested with sine waves on periodic domains on a cartesian regular-hexagon mesh, as described in [Bishnu et al. 2023b](https://www.authorea.com/doi/full/10.22541/essoar.167100170.03833124/v1) (note: update links later when JAMES articles are published) and [Bishnu et al. 2021](https://doi.org/10.5281/zenodo.7439539).

1. Divergence operator at cell centers (2nd order)
1. Gradient operator normal to edges (2nd order)
1. Curl operator at vertices, i.e. vorticity from a vector field (1st order)
1. Curl operator interpolated to cell centers (2nd order)
1. Tangential velocity at edges, computed from normal velocity (2nd order)

Requirements for tests are:

- compares order of convergence against an expected threshold

These tests can be conducted at the earliest stages of dycore development, in tandem with the implementation of each operator. In the longer term, these tests will be subsumed by the inertia gravity wave, which tests the order of convergence of all these operators together.

Operator convergence tests may additionally be conducted on the sphere using spherical harmonics for the analytic solution. This was conducted on MPAS-Ocean by Hyun Kang in 2023.

### 5.2 Inertia Gravity Wave: linearized shallow water, no tracers

The inertia gravity wave test provides an exact solution in time for the linearized shallow water equations (momentum and thickness). It tests the time stepping scheme along with the linearized advection terms, the SSH gradient. It does not test diffusive terms or bottom drag. It is conducted on a doubly periodic cartesian mesh, so does not test boundary conditions. The numerical solution should converge to the analytic solution at 2nd order, as shown in [Bishnu et al. 2023b](https://www.authorea.com/doi/full/10.22541/essoar.167100170.03833124/v1) (add figure number and link when published). Also see the inertia gravity test case in Polaris.

Requirements for tests are:

- expected order of convergence
- conservation of total volume
- automation and reproducibility in polaris

This test should be conducted as soon as momentum and thickness equations and time-stepping is in place. It does not require any tracer infrastructure.

### 5.3 Manufactured Solution: full nonlinear shallow water, no tracers

The manufactured solution test provides an exact solution in time for the full nonlinear shallow water equations (momentum and thickness). It tests the time stepping scheme along with the the full advection terms and the SSH gradient.  It is conducted on a doubly periodic cartesian mesh, so does not test boundary conditions. Error should converge at 2nd order, as shown in [Bishnu et al. 2023b](https://www.authorea.com/doi/full/10.22541/essoar.167100170.03833124/v1) (add figure number and link when published). Also see the manufactured solution test case in Polaris.

Requirements for tests are:

- expected order of convergence
- conservation of total volume
- automation and reproducibility in polaris

This test should be conducted as soon as momentum and thickness equations and time-stepping is in place. It does not require any tracer infrastructure.

### 5.4 Tracer transport on a sphere

A test suite will be used to test horizontal transport schemes on the sphere comprised of test case 1 from [Williamson et al. (1992)](https://www.sciencedirect.com/science/article/pii/S0021999105800166) and several tests from [Lauritzen et al. (2012)](https://www.geosci-model-dev.net/5/887/2012/).
They test the tracer advection and tracer time stepping. Velocity and thickness are predefined fields and remain fixed, so this does not exercise those equations.

Requirements for tests are:

- compares order of convergence against an expected threshold
- conservation of total tracer amount
- produces visualization that allows the user to evaluate whether numerical mixing resembles real mixing (preserves functional relationships between tracers)
- tracer min and max values should remain bounded by initial bounds.
- automation and reproducibility in polaris

### 5.5 Performance testing

Tests can be conducted with inertia-gravity wave tests but with full non-linear terms and 100 identical layers. Domain will be cartesian resolutions from 64x64 up to 512x512. See performance test results using MPAS-Ocean on Perlmutter, using this setup, in [Bishnu et al. 2023a](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-57). Additionally, we could test spherical cases with a Williamson test case or similar.

Requirements for tests are:

- Single CPU performance should be as good or better than single CPU MPAS-Ocean.
- GPU tests on single node should be better than single node CPU tests.
- Scaling on CPUs multi-node, up to 4096 cores or higher, should be close to perfect scaling for 512x512 mesh, and as good or better than MPAS-Ocean.
- Scaling on GPUs for multi-node should be close to perfect scaling for 512x512 mesh.
- automation and reproducibility in polaris

### 5.6 Shallow water tests on spherical domains and with realistic coastlines

Global cases will facilitate the testing of shallow water dynamics on the sphere and of boundary conditions with realistic coastlines, which are not included in any of the previous tests. This includes testing conservation of volume and tracers in global domains.

Potential test cases include:

1. Surface gravity wave: The speed of the surface gravity wave can be compared to the theoretical expectation, as shown in [Pal et al. 2023](https://gmd.copernicus.org/articles/16/1297/2023) Appendix A. This could use an aquaplanet domain (flat bottom, no coastlines) to measure the gravity wave speed, and a realistic domain to test boundary conditions and variable-depth bathymetry.
1. Stommel double gyre:  This may be compared to an exact solution in the Cartesian case, as in [Pal et al. 2023](https://gmd.copernicus.org/articles/16/1297/2023) Appendix B, or qualitative comparisons for the spherical case using either idealized boundaries or an isolated Atlantic Basin domain.

Tests of realistic global circulation cannot be done with the shallow water equations of Omega-0, but will be part of Omega-1 development with the layered primitive equation model.

### 5.7 Further tests for Omega-0

Potential additional test cases include the following. These are useful to explore and validate model behavior, but are optional, depending on the time available.

1. [Williamson et al. 1992](https://www.sciencedirect.com/science/article/pii/S0021999105800166) as test case 2 or 3, global steady state solutions.
1. Global case with unstable jet [Galewsky et al. 2003](https://doi.org/10.3402/tellusa.v56i5.14436)
