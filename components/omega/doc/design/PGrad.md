(omega-design-pressure-grad)=
# Pressure Gradient


## 1 Overview
The pressure gradient will be responsible for computing the horizontal gradients of both the pressure and geopotential terms for the non-Boussinesq primative equations implemented in Omega.
In the non-Boussinesq model, the vertical coordinate will be pressure as opposed to height.
In a pure pressure coordinate the pressure gradient term disappears (since the pressure does not vary along lines of constant pressure), just as how the geopotential term disappears in a pure z coordinate model.
However, similar to MPAS-Ocean's support for tilted height coordinates, Omega will allow for tilted pressure coordinates.
This means that Omega will need to compute both the pressure and geopotential gradients.

## 2 Requirements

### 2.1 Requirement: Support for tilted pressure coordinates for the non-Boussinesq primitive equations

The pressure gradient will compute the horizontal gradients of both the pressure and geopotential to support tilted pressure coordinates in the non-Boussinesq model.
This will allow for the use of a $p^\star$ coordinate, which functions similarly to the $z^\star$ in the Boussinesq MPAS-Ocean model.

### 2.2 Requirement: Initial support for a simple centered pressure gradient

For initial global cases without ice shelf cavities, the pressure and geopotential gradients will be computed with a simple centered difference approximation. In later versions of Omega, one ore more  high-order pressure gradients will be implemented and will replace the centered approach in production runs.
However, the centered pressure gradient will remain an option for use in idealized testing.

### 2.3 Requirement: Flexibility to support a high-order pressure gradient
The centered pressure gradient will be insufficient for future versions of Omega that include ice shelf cavities and high resolution shelf breaks.
The pressure gradient framework should be flexible enough to support a high-order pressure gradient in the future.

### 2.4 Requirement: Flexibility to support tidal forcing and sea level change
In later versions of Omega, the pressure gradient will need to be able to include tidal forcing in the geopotential term.
These tidal forcings include both the tidal potential and the self attraction and loading terms.
Additionally, other changes to the geoid

### 2.5 Requirement: Pressure gradient for barotropic mode

For split barotropic-baroclinic timestepping, the pressure gradient should provide the bottom pressure gradient tendency in the barotropic mode.

### Desired:

## 3 Algorithmic Formulation
The non-Boussinesq momentum equation is
$$ \frac{D \mathbf{u}_h}{D t } + f\boldsymbol{k}\times \mathbf{u}_h + \left(v\nabla_A p + \nabla_A \phi \right) = \boldsymbol{\mathcal{F}}. $$

where $\mathbf{u}_h$ is the horizontal velocity, $f$ is the Coriolis parameter, $v = \frac{1}{\rho}$ is the specific volume, $\rho = \rho(T,S,p)$ is the density, $p$ is the hydrostatic pressure, $\phi$ is the geopotential, and $\boldsymbol{\mathcal{F}}$ are the dissipative terms.
The operator $\nabla_A$ is the gradient along a constant surface, $A$, and the total derivative is

$$ \frac{D \mathbf{u}_h}{D t } = \left( \frac{\partial}{\partial t} \right)_A  + \mathbf{u}_h\cdot \nabla_A + \omega\frac{\partial}{\partial A}, $$

where $\omega$ is the cross coordinate flow.
In the layered non-Boussinesq equations, the prognostic variable is the pressure thickness $h_k$, so that the geometric thickness (in meters) is a diagnostic variable defined as:

$$ \Delta z_k = v_k h_k. $$

The pressure at vertical cell interfaces is the found by summing the pressure thicknesses:

$$ p_{K+1/2} = p_{surf} + g\sum_{k=1}^K h_k. $$

The geopotential at vertical cell integrates is found by summing the pressure thicknesses multiplied by the specific volume:

$$ \phi_{K+1/2} = g\left(z_b + \sum_{k=K}^{N} v_k h_k\right), $$

where $z_b$ is the (positive-up) bottom depth.

The discrete gradient operator at an edge is:

$$ \nabla {(\cdot)} = \frac{1}{d_e} \sum_{i\in CE(e)} -n_{e,i} (\cdot)_i $$

where $d_e$ is the distance between cell centers, $CE(e)$ are the cells on edge $e$, and $n_{e,i}$ is the sign of the edge normal with respect to cell $i$.
Therefore the centered pressure gradient will be calculated as:

$$ T^p_k = \frac{1}{d_e} \left( \widehat{v}_{k,e} \sum_{i \in CE(e)} n_{e,i} \overline{p}_{k,i} + \sum_{i\in CE(e)} n_{e,i} \overline{\phi}_{k,i}\right), $$

$$ = \frac{1}{d_e} \left(  \sum_{i \in CE(e)} n_{e,i} \left( \widehat{v}_{k,e}\overline{p}_{k,i} + \overline{\phi}_{k,i} \right) \right), $$

with the vertical averaging operator defined as:

$$ \overline{(\cdot)} = \frac{1}{2} \left((\cdot)_{k+1/2} + (\cdot)_{k-1/2} \right)$$

and the horizontal averaging operator:

$$\widehat{(\cdot)} = \frac{1}{2} \left( (\cdot)_{i=1} + (\cdot)_{i=2}\right)$$

## 4 Design

### 4.1 Data types and parameters

The `PressureGradient` class will be used to perform the horizontal gradients of the pressure and geopotential
```c++
class PressureGrad{
    public:
    private:
        PressureGradCentered CenteredPGrad;
        PressureGradHighOrder HighOrderPGrad;
        PressureGradType PressureGradChoice;
        I4 NVertLevels;
        I4 NChuncks;
        Array2DI4 CellsOnEdge;
        Array1DReal DvEdge;
        Array1DReal EdgeSignOnCell;
}
```
The user will select a pressure gradient option at runtime in the input yaml file under the pressure gradient section:
```yaml
    PressureGrad:
       PressureGradType: 'centered'
```
An `enum class` will be used to specify options for the pressure gradient used for an Omega simulation:
```c++
enum class PressureGradType{
   Centered,
   HighOrder
}
```
The functions to compute the centered and high order pressure gradient terms will be implemented as functors and the pressure gradient class will have private instances of these classes.
### 4.2 Methods

```c++
class PressureGrad{
    public:
        static PressureGrad *create();
        void computePressureGrad();
        void computePressureGradBtr();
    private:

}
```
#### 4.2.1 Creation

The constructor will be responsible for storing any static mesh information as private variables and handling the selection of the user-specified pressure gradient option.
```c++
PressureGrad::PressureGrad(const HorzMesh *Mesh, int NVertLevels, Config *Options);
```

The create method will take the same arguments as the constructor plus a name.
It calls the constructor to create a new pressure gradient instance, and put it in the static map of all pressure gradients.
It will return a pointer to the newly created object.
```c++
PressureGrad *PressureGrad::create(const std::string &Name, const HorzMesh *Mesh, int NVertLevels, Config *Options);
```

#### 4.2.2 Initialization

The init method will create the default pressure gradient and return an error code:
```c++
int PressureGrad::init();
```

#### 4.2.3 Retrieval

There will be methods for getting the default and non-default pressure gradient instances:
```c++
PressureGrad *PressureGrad::getDefault();
PressureGrad *PressureGrad::get(const std::string &Name);
```

#### 4.2.4 Computation

The public `computePressureGrad` method will rely on private methods for each specific pressure gradient option (centered and high order).
```c++
void PressureGrad::computePressureGrad(const Array2DReal &Tend,
                                       const Array2DReal &Pressure,
                                       const Array2DReal &Geopotential,
                                       const Array2DReal &SpecVol) {
OMEGA_SCOPE(LocCenteredPGrad, CenteredPGrad)
OMEGA_SCOPE(LocHighOrderPGrad, HighOrderPGrad)
    if (PressureGradChoice == PressureGradType::Centered){
       parallelFor("pgrad-centered", {NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocCenteredPGrad(Tend, Pressure, Geopotential, SpecVol);
       });
    }
    else if (PressureGradChoice == PressureGradType::HighOrder){
       parallelFor("pgrad-highorder", {NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocHighOrderPGrad(Tend, Pressure, Geopotential, SpecVol);
        });
    }
```

#### 4.2.5 Destruction and removal

No operations are needed in the destructor.
The erase method will remove a named pressure gradient instance, whereas the clear method will remove all of
them.
Both will call the destructor in the process.
```c++
void PressureGrad::erase(const std::string &Name);
void PressureGrad::clear();
```

## Verification and Testing

### Test: Spatial convergence to exact solution
For a given analytical $v$,  $h$, and $b$, the spatial convergence of the pressure gradient can be assessed by computing errors on progressively finer meshes.

### Test: Baroclinic gyre
The baroclinic gyre test case will test the pressure gradient term in the full non-Boussinesq equations.
