(omega-design-my-class)=
# ShallowWaterOmega0

## 1 Overview

This design document describes the first version of the Omega ocean model, Omega0. Overall, Omega is an unstructured-mesh ocean model based on TRiSK numerical methods (ref) that is specifically designed for modern exascale computing architectures. The algorithms in Omega will be nearly identical to those in MPAS-Ocean, but it will be written in c++ rather than Fortran in order to take advantage of libraries to run on GPUs, such as YAKL (Yet Another Kernel Library, ref).

The planned versions of Omega are:

1. **Omega-0: Shallow water equations with identical vertical layers and inactive tracers.** There is no vertical transport or advection. The tracer equation is horizontal advection-diffusion, but tracers do not feed back to dynamics. Pressure gradient is simply gradient of sea surface height. Capability is similar to [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780)
1. **Omega-1: Minimal stand-alone primitive equation eddy-permitting model.** This adds active temperature, salinity, density, and pressure as a function of depth. There is a true pressure gradient. Vertical velocity is from the continuity equation. An equation of state and simple vertical mixing scheme are needed. Capability is same as [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760)
1. **Omega-2: Eddying ocean with advanced parameterizations.** This model will have sufficient capability to run realistic global simulations, similar to E3SM V1 [Petersen et al. 2019]().


We will produce separate design documents for the time-stepping scheme and the tracer equations in Omega-0. Some of the requirements stated here are repeated or implied by the Omega framework design documents, but are included for clarity.

## 2 Requirements

### 2.1 Requirement: Omega-0 will solve the nonlinear shallow water equations, plus inactive tracers

The governing equations for Omega-0 are conservation of momentum, volume, and tracers in a single layer:

$$
\frac{\partial \boldsymbol{u}}{\partial t} + q\left(h\boldsymbol{u}^{\perp}\right) = -g\nabla(h+b) - \nabla K + \nu_2 \nabla^2 \boldsymbol{u} - \nu_4 \nabla^4 \boldsymbol{u} - c\boldsymbol{u}\left|\boldsymbol{u}\right|
$$

$$
\frac{\partial h}{\partial t} + \nabla \cdot \left(h \boldsymbol{u}\right) = 0,
$$

$$
\frac{\partial h \phi}{\partial t} + \nabla \cdot \left(h \boldsymbol{u} \phi\right) = \kappa_2 h \nabla^2 \phi - \kappa_4 h \nabla^4 \phi
$$

where the first two equations are from [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780), equations 2 and 7, with additional viscosity and bottom drag terms.
This equation set does not include any vertical advection or diffusion. Omega-0 will have a vertical index for performance testing and future expansion, but vertical layers will be simply redundant.

To do: Complete table of symbol definitions. See [markdown table generator](https://www.tablesgenerator.com/markdown_tables)
| symbol | name      | units  |
|--------|-----------|--------|
| $u$    | velocity  | m/s    |
| $h$    | thickness | m      |
| $\phi$ | tracer    | varies |
| $t$    | time      | s      |

### 2.2 Requirement: Numerical method will be the TRiSK formulation

The horizontal discretization will be taken from 
[Thuburn et al. 2009](https://www.sciencedirect.com/science/article/pii/S0021999109004434) and [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780), as described in detail in the algorithmic formulation in Section 3 below. This is the same base formulation as in MPAS-Ocean. In addition, we will consider small alterations from the original MPAS-Ocean horizontal discretization and include them as options if they are of minimal additional effort. This includes the recent AUST formulation in [Calandrini et al. 2021](https://www.sciencedirect.com/science/article/pii/S146350032100161X) and simple vorticity averaging considered by the Omega team this past year.

### 2.3 Requirement: Omega-0 will use MPAS format unstructured-mesh domains

We will continue to use the MPAS-format netcdf files for input and output, with the same mesh variable names and dimensions. This will facilitate ease of use and interoperability between MPAS-Ocean and Omega.

Omega mesh information will be stored in a separate mesh file, and not be included with initial condition, restart, or output files. This way there is never redundant mesh data stored in these files.

To do: add list of minimum number of variables that defines the mesh, and the additional variables that can be computed on start-up.

### 2.4 Requirement: Omega-0 will interface with polaris for preprocessing and postprocessing

A substantial investment has already been made in the polaris tools. Continuing the use of polaris for Omega will speed development and encourage the documentation and long-term reproducibility of test cases.

The test cases relevant to this design document are in Section 5 below.

### 2.5 Requirement: Omega-0 will run portably on various DOE architectures (CPU and GPU nodes)

Omega will be able to run on all the upcoming DOE architectures and make good use of GPU hardware. This should occur with minimal alterations in the high-level PDE solver code for different platforms. 

The proposed solution is to use YAKL [(Norman et al. 2023)](https://link.springer.com/10.1007/s10766-022-00739-0) for within-node, shared memory computations. All large arrays will be YAKL objects, and computations that loop over array indices will use YAKL function calls.

### 2.6 Requirement: Omega-0 will run on multi-node with domain decomposition

This is with MPI and halo communication, as described in framework design documents.

### 2.7 Requirement: Correct convergence rates of operators and exact solutions

Each operator will be tested individually for the proper convergence rate: divergence, gradient, curl interpolated to cell centers, and tangential velocity are all second order; curl at vertices is first order. The details of the convergence tests are explained in section 5.1 below.

### 2.8 Requirement: Conservation of volume and tracer

The total volume of the domain should be conserved to machine precision in the absence of surface volume fluxes. Likewise, the total tracer amount should be conserved to machine precision without surface fluxes. A simple test is to initialize an inactive test tracer with a uniform value of 1.0, and it should remain 1.0 throughout the domain.

### 2.9 Requirement: Performance will be at least as good or better than MPAS-Ocean

Performance will assessed with the following metrics:

1. Single CPU throughput
1. Parallel CPU scalability to high node counts
1. Single GPU throughput

The first two items will be compared to MPAS-Ocean with the same test configuration.

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

The following terms can be tested with sine waves on periodic domains on a cartesian regular-hexagon mesh, as described in (ref Bishnu thesis) and (ref Julia paper).

1. Divergence operator at cell centers (2nd order)
1. Gradient operator normal to edges (2nd order)
1. Curl operator at vertices, i.e. vorticity from a vector field (1st order)
1. Curl operator interpolated to cell centers (2nd order)
1. Tangential velocity at edges, computed from normal velocity (2nd order)

Requirements for tests are:

- expected order of convergence

These tests can be conducted at the earliest stages of dycore development, in tandem with the implementation of each operator.

Note: consider applying this as a unit test.

Note: We may not need this group of tests long-term, because the inertia gravity wave tests all of them together.

Question: should we make tests for spherical domain here as well? Hyun has done those before with spherical harmonics.

### 5.2 Inertia Gravity Wave: linearized shallow water, no tracers

The inertia gravity wave test provides an exact solution in time for the linearized shallow water equations (momentum and thickness). It tests the time stepping scheme along with the linearized advection terms, the SSH gradient. It does not test diffusive terms or bottom drag. It is conducted on a doubly periodic cartesian mesh, so does not test boundary conditions. Error should converge at 1st order, as described by [Bishnu 2021](https://zenodo.org/record/7439539) and Bishnu et al. 2023. Also see the inertia gravity test case in Polaris.

Requirements for tests are:

- expected order of convergence
- conservation of total volume
- automation and reproducibility in polaris

This test should be conducted as soon as momentum and thickness equations and time-stepping is in place. It does not require any tracer infrastructure.

### 5.3 Manufactured Solution: full nonlinear shallow water, no tracers

The manufactured solution test provides an exact solution in time for the full nonlinear shallow water equations (momentum and thickness). It tests the time stepping scheme along with the the full advection terms and the SSH gradient.  It is conducted on a doubly periodic cartesian mesh, so does not test boundary conditions. Error should converge at 1st order, as described by [Bishnu 2021](https://zenodo.org/record/7439539) and Bishnu et al. 2023. Also see the manufactured solution test case in Polaris.

Requirements for tests are:

- expected order of convergence
- conservation of total volume
- automation and reproducibility in polaris

This test should be conducted as soon as momentum and thickness equations and time-stepping is in place. It does not require any tracer infrastructure.

### 5.4 Advection of cosine bell on a sphere

Advection of a tracer cosine bell around the sphere was described in [Williamson et al. 1992](https://www.sciencedirect.com/science/article/pii/S0021999105800166) as test case 1. It tests the tracer advection and tracer time stepping. Velocity and thickness are predefined fields and remain fixed, so this does not exercise those modules.

Requirements for tests are:

- expected order of convergence (check)
- conservation of total tracer amount
- tracer min and max values should remain bounded by initial bounds.
- automation and reproducibility in polaris

### 5.5 Performance testing

Tests can be conducted with inertia-gravity wave tests but with full non-linear terms and 100 identical layers. Domain will be cartesian resolutions from 64x64 up to 512x512. These can be compared with Bishnu 2023 (ref). Additionally, we could test spherical cases with a Williamson test case or similar.

Requirements for tests are:

- Single CPU performance should be as good or better than single CPU MPAS-Ocean.
- GPU tests on single node should be better than single node CPU tests.
- Scaling on CPUs multi-node, up to 4096 cores or higher, should be close to perfect scaling for 512x512 mesh, and as good or better than MPAS-Ocean.
- Scaling on GPUs for multi-node should be close to perfect scaling for 512x512 mesh.
- automation and reproducibility in polaris

### 5.6 Demonstration cases with realistic coastlines

Tests in realistic domains will facilitate the testing of boundary conditions on realistic coastlines, which are not included in any of the previous tests. They will also be useful for highlights and talks to show that Omega-0 can run in global realistic domains, albeit with shallow water equations only.

Potential test cases include:

1. Coastal inundation: Delaware Bay mesh, after ICOM tests, to demonstrate boundary conditions, resolution, and for 'splashy' images. Would tidal forcing be required for that?
1. Global realistic Tsunami to test boundary conditions and grid-scale noise. See Darren Engwirda's recent tests on Alaska's Aleutian Islands.
1. Test for boundary conditions - perhaps Stommel double gyre in [Pal et al. 2023](https://gmd.copernicus.org/articles/16/1297/2023). 

### 5.7 Further tests for Omega-0

Potential additional test cases:

1. [Williamson et al. 1992](https://www.sciencedirect.com/science/article/pii/S0021999105800166) as test case 2 or 3, global steady state solutions.
1. Global case with unstable jet [Galewsky et al. 2003](https://doi.org/10.3402/tellusa.v56i5.14436)
1. Sphere transport [Lauritzen et al. 2012](https://doi.org/10.5194/gmd-5-887-2012)
