(omega-design-horz-mesh)=
# Mesh 

## 1 Overview

The mesh class will contain the YAKL arrays which describe the horizontal mesh and are used in the computation of the tendency terms in the discrete governing equations. OMEGA will separate the horizontal mesh variables from the vertical mesh information.

## 2 Requirements

### 2.1 Requirement: OMEGA will use the previously established MPAS Mesh Spec

The OMEGA mesh information should be compatible with the [MPAS Mesh Specification](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf).

### 2.2 Requirement: Functionality is needed to read the mesh on the host and transfer relevant data to the device for computation

Not all mesh information is required in computing the tendency terms on the device, e.g. lonCell, latCell, etc. 
However, other arrays will need to be allocated and copied to the device for use in tendency computation.
The mesh class will explicitly include host and device YAKL arrays for each variable.
A class method will be included to copy the relevant mesh information to the device.

### 2.3 Requirement: Zero-based cell, edge, and vertex numbering

Although the existing MPAS Mesh spec uses a one-based mesh numbering, zero-based mesh numbering is required to be compatible with the zero-based indexing used for YAKL arrays.

### 2.4 Requirement: Work with Decomp class to decompose mesh

The mesh class will reference the partitioned connectivity arrays created by the Decomp class.

### 2.5 Requirement: Mesh variables will be associated with metadata to describe data 

Following the Metadata and IO designs, the YAKL arrays for the mesh variables will be associated with information about the represented values.

### 2.6 Requirement: I/O to obtain mesh data

The Mesh class will have a method to read in the mesh information not obtained by the Decomp class. 

### 2.7 Desired: Ability to support multiple independent mesh objects

This flexibility is required to support future implementations of spatially split barotropic and baroclinic modes that are computed on different resolution meshes.
Additionally, this flexibility can be used to support separate domain decompositions for the barotropic and baroclinic meshes, which can help optimize the communication frequency for the barotropic subcycling via wide barotropic halos.

### 2.8 Desired: OMEGA can read in a reduced number of mesh variables and compute the remaining array information online.

Many of the mesh variables are not independent, e.g.  areaCell, weightsOnEdge, etc., and can be computed from a reduced set of mesh variables.
This functionality could be used to reduce mesh/restart file size for high resolution meshes.  
Building in this flexibility would allow for all mesh-related calculations to be available within the code base instead of spreading them amongst various utility programs in MPAS-Tools.
This will improve the ability to maintain and unit test the mesh calculations in MPAS-Tools.
The functions used to compute the dependent mesh information can also be used to verify any mesh information that is provided in the mesh input stream.

As standard practice, all necessary internally computed mesh information will be output in a single file for post-processing purposes.
A checksumming strategy will be implemented to avoid situations where simulation data and mesh information are mismatched during post processing. 

Where appropriate, some additional derived quantities (e.g. reciprocals) will also be included to improve the performance of device-side calculations.

## 3 Algorithmic Formulation


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
