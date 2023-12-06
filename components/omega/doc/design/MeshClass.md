(omega-design-mesh)=
# Mesh 

## 1 Overview

The mesh class will contain the YAKL arrays which describe the mesh and are used in the computation of the tendency terms in the discrete governing equations.

## 2 Requirements

### 2.1 Requirement: OMEGA will use the previously estabilished MPAS Mesh Spec

The OMEGA mesh information should be compatiable with the [MPAS Mesh Specification](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf).

### 2.2 Requirement: Functionality is needed to read the mesh on the host and transfer relevant data to the device for comutation

Not all mesh information is required in computing the tendency terms on the device, e.g. lonCell, latCell, etc. 
However, other arrays will need to allocated and copied to the device for use in tendency computation.
One option is to template the mesh class so a mesh object can be created for both the host and device, with a method to perform the copy between the two object types.

### 2.3 Requirement: Zero-based cell, edge, and vertex numbering

Although the existing MPAS Mesh spec uses a one-based mesh numbering, zero-based mesh numbering is required to be compatiable with the zero-based indexing used for YAKL arrays.

### 2.4 Desired: Ability to support multiple independent mesh objects

This flexibility is required to support future implmentations of spatially split barotropic and baroclinic modes that are computed on different resolution 

### 2.5 Desired: OMEGA can read in a reduced number of mesh variables and compute the remaning array information online.

Many of the mesh variables are not independent, e.g.  areaCell, weightsOnEdge, etc., and can be computed from a reduced set of mesh variables.
This functionality could be used to reduce mesh/restart file size for high resoultion meshes.  
Builiding in this flexibility would allow for all mesh-related calculations to be availiable within the code base instead of spreading them amonst various utility programs in MPAS-Tools.
The functions used to compute the dependent mesh information can also be used to verify mesh information provided in the mesh input stream.

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
