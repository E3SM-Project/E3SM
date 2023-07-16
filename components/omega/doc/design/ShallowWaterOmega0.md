(omega-design-my-class)=
# ShallowWaterOmega0

## 1 Overview

This design document describes the first version of the Omega ocean model, Omega0. Overall, Omega is an unstructured-mesh ocean model based on TRiSK numerical methods (ref) that is specifically designed for modern exascale computing architectures. The algorithms in Omega will be nearly identical to those in MPAS-Ocean, but it will be written in c++ rather than Fortran in order to take advantage of libraries to run on GPUs, such as YAKL (Yet Another Kernal Library, ref). 

The planned versions of Omega are:
1. **Omega-0: Shallow water equations with identical vertical layers and inactive tracers.** There is no vertical transport or advection. The tracer equation is horizontal advection-diffusion, but tracers do not feed back to dynamics. Pressure gradient is simply gradient of sea surface height. Capability is similar to Ringler et al. 2010.
1. **Omega-1: Minimum stand-alone primitive equation eddy-permitting model.** This adds active temperature, salinity, density, and pressure as a function of depth. There is a true pressure gradient. Vertical velocity is from the continuity equation. An equation of state and simple vertical mixing scheme are needed. Capability is same as Ringler et al. 2013.
1. **Omega-2: Eddying ocean with advanced parameterizations.** This model will have sufficient capability to run realistic global simulations, similar to E3SM V1 (Petersen et al. 2019).


We will produce separate design documents for the time-stepping scheme and the tracer equations in Omega-0. Some of the requirements stated here are repeated or implied by the Omega framework design documents, but are included for clarity.

## 2 Requirements

### 2.1 Requirement: Omega-0 will solve the nonlinear shallow water equations, plus inactive tracers.

The governing equations for Omega-0 are conservation of momentum, volume, and tracers in a single layer:

$$
\frac{\partial \boldsymbol{u}}{\partial t} + q\left(h\boldsymbol{u}^{\perp}\right) = -g\nabla(h+b) - \nabla K + \nu \nabla^2 \boldsymbol{u},
$$

$$
\frac{\partial h}{\partial t} + \nabla \cdot \left(h \boldsymbol{u}\right) = 0,
$$

$$
\frac{\partial h \phi}{\partial t} + \nabla \cdot \left(h \boldsymbol{u} \phi\right) = \kappa h \nabla^2 \phi
$$

where the first two equations are from Ringler et al. 2010, equations 2 and 7. (Add symbol definitions)

This equation set does not include any vertical advection or diffusion. Omega-0 will have a vertical index for performance testing and future expansion, but vertical layers will be simply redundant.

### 2.2 Requirement: Numerical method will be the TRiSK formulation

Repeat numerical formulation here, from 2009, 2010 papers. Also add small options for extra recent formulations by Darren and Sara (AUST)

### 2.3 Requirement: Omega-0 will use MPAS format unstructured-mesh domains

We will continue to use the MPAS-format netcdf files for input and output, with the same mesh variable names and dimensions. This will facilitate ease of use and interoperability between MPAS-Ocean and Omega.

### 2.4 Requirement: Omega-0 will interface with compass/polaris for preprocessing and post-processing.

A substantial investment has already been made in the compass/polaris tools. Continuing the use of compass/polaris for Omega will speed development and encourage the documentation and long-term reproducability of test cases.

### 2.5 Requirement: Omega-0 will run on CPUs and GPUs with YAKL

All loops over horizontal and vertical indices will use YAKL to facilitate the use of both CPUs and GPUs.

### 2.6 Requirement: Omega-0 will run on multi-node with domain decomposition

This is with MPI and halo communication, as described in framework design documents.

## 3 Algorithmic Formulation

add text

## 4 Design

You can include code blocks like this:

```c++
int var = value;
```

### 4.1 Data types and parameters

#### 4.1.1 Parameters 

List and define any configuration parameters or public constants.

#### 4.1.2 Class/structs/data types

Describe any public data types and/or the class definition

### 4.2 Methods

List and describe all public methods and their interfaces (actual code for 
interface that would be in header file). Describe typical use cases.

## 5 Verification and Testing

### 5.1 Test xxx

Describe test including conditions for pass/fail
List which requirements it tests: 
  - tests requirement xxx

### 5.2 Test yyy

Describe test including conditions for pass/fail
List which requirements it tests: 
  - tests requirement zzz
