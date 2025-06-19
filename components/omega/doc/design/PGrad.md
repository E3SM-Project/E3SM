(omega-design-pressure-grad)=
# Pressure Gradient


## 1 Overview
The pressure gradient will be responsible for computing the horizontal gradients of both the pressure and geopotential terms for the non-Boussinesq primitive equations implemented in Omega.
In the non-Boussinesq model, the conserved quantity is mass rather than volume.
In Omega the prognostic variable $\tilde{h}$ is a pseudo thickness, rather than geometric thickness in m as in a Boussinesq model.
 Some non-Boussinesq models are written in pressure coordinates (e.g. [de Szoeke and Samelson 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml).
 However, Omega is written in general vertical coordinates and can reference either pressure $p$ or distance $z$ in the vertical.
In a pure pressure coordinate the pressure gradient term disappears (since the pressure does not vary along lines of constant pressure), just as how the geopotential term disappears in a pure z coordinate model.
However, similar to MPAS-Ocean's support for tilted height coordinates, Omega will allow for tilted pressure coordinates.
This means that Omega will need to compute both the pressure and geopotential gradients.

## 2 Requirements

### 2.1 Requirement: Support for tilted pressure coordinates for the non-Boussinesq primitive equations

The pressure gradient will compute the horizontal gradients of both the pressure and geopotential to support tilted pressure coordinates in the non-Boussinesq model.
This will allow for the use of a $p^\star$ coordinate, which functions similarly to the $z^\star$ in the Boussinesq MPAS-Ocean model.

### 2.2 Requirement: Initial support for a simple centered pressure gradient

For initial global cases without ice shelf cavities, the pressure and geopotential gradients will be computed with a simple centered difference approximation. In later versions of Omega, one or more  high-order pressure gradients will be implemented and will replace the centered approach in production runs.
However, the centered pressure gradient will remain an option for use in idealized testing.

### 2.3 Requirement: Flexibility to support a high-order pressure gradient
The centered pressure gradient will be insufficient for future versions of Omega that include ice shelf cavities and high resolution shelf breaks.
The pressure gradient framework should be flexible enough to support a high-order pressure gradient in the future.
The high order pressure gradient will be similar to [Adcroft et al. 2008](https://doi.org/10.1016/j.ocemod.2008.02.001).

### 2.4 Requirement: Flexibility to support tidal forcing and sea level change
In later versions of Omega, the pressure gradient will need to be able to include tidal forcing in the geopotential term.
These tidal forcings include both the tidal potential and the self attraction and loading terms.
Additionally, other long-term changes to the geoid can be included in the geopotential.

### 2.5 Disired: Pressure gradient for barotropic mode

For split barotropic-baroclinic timestepping, the pressure gradient should provide the bottom pressure gradient tendency in the barotropic mode.
This will be added in a future version when split time stepping is implemented.

### Desired:

## 3 Algorithmic Formulation
### 3.1 Centered Pressure Gradient
In the layered non-Boussinesq [momentum equation](OmegaV1GoverningEqns.md#discrete-momentum) solved in Omega, the pressure gradient tendency term for edge $e$ and level $k$, $T^p_{e,k}$, includes the gradient of the pressure and the gradient of the geopotential,

$$
T^p_{e,k} = -\left[ \alpha_{i,k} \right]_e \nabla p_{i,k} - \nabla \Phi_{i,k},
$$

where the second term is necessary to account for tilted layers that occur when using a general vertical coordinate.
In this equation, $\alpha_{i,k}$ is the specific volume for cell $i$ at the mid-point of level $k$, $p_{i,k}$ is the pressure, and $\Phi_{i,k}$ is the geopotential.
The discrete gradient operator at an edge is:

$$
 \nabla {(\cdot)} = \frac{1}{d_e} \sum_{i\in CE(e)} -n_{e,i} (\cdot)_i
$$

where $d_e$ is the distance between cell centers, $CE(e)$ are the cells on edge $e$, and $n_{e,i}$ is the sign of the edge normal with respect to cell $i$.
The horizontal averaging operator is:

$$
 [\cdot]_e = \frac{1}{2}\sum_{i\in CE(e)} (\cdot)_i
$$

Therefore, the centered pressure gradient will be calculated as:

$$
 T^p_{e,k} = \frac{1}{d_e} \left( [\alpha_{i,k}]_e \sum_{i \in CE(e)} n_{e,i} p_{k,i} + \sum_{i\in CE(e)} n_{e,i} \Phi_{k,i}\right),
$$

$$
 = \frac{1}{d_e} \left(  \sum_{i \in CE(e)} n_{e,i} \left( [\alpha_{i,k}]_ep_{k,i} + \Phi_{k,i} \right) \right),
$$

### 3.2 Barotropic Pressure Gradient

When split baroclinic-barotropic time stepping is implemented in the future, the barotropic pressure gradient will be calculated by the pressure gradient class.
The barotropic pressure gradient is found by depth integrating the pressure gradient.
The pressure is

$$
p(z) = p_b - g \int^z_{-h} \rho dz^\prime,
$$

where $p_b$ is the bottom pressure.
The bottom pressure is the sum of the atmospheric surface pressure, $p_s$, and the pressure contribution of the water column:

$$
p(z) &=  p_s + g\int_{-h}^\eta \rho dz - g \int^z_{-h} \rho dz^\prime, \\
     &=  p_s + g\rho_0\widetilde{H} - g \int^z_{-h} \rho dz^\prime,
$$

where the total water column pseudo height is expressed by

$$
\widetilde{H} = \int_{-h}^\eta \frac{\rho}{\rho_0} dz.
$$

$\widetilde{H}$ is the prognositc variable in the barotropic continuity equation.
The vertical integral of the pressure gradient is

$$
\frac{1}{\rho_0\widetilde{H}}\int^\eta_{-h} \nabla p dz &= \frac{1}{\rho_0\widetilde{H}}\int^\eta_{-h} \nabla \left( p_s + g\rho_0 \widetilde{H} - g \int^z_{-h} \rho dz^\prime \right) dz, \\
                           &= \frac{H}{\rho_0\widetilde{H}}\nabla p_s + \frac{gH}{\widetilde{H}}\nabla \widetilde{H} - \frac{g}{\rho_0\widetilde{H}} \int_{-h}^\eta \left( \nabla \int_{-h}^z \rho dz^\prime\right) dz,
$$

where the height of the water column is represented by $H$.
The $1/\rho_0\widetilde{H}$ factor comes vertically integrating the material derivative and expressing the resulting barotropic momentum equation in non-conservative form.

Therefore, the barotorpic pressure gradient term is discretized as:

$$
\overline{T}_e^p = g\left[ \frac{H_i}{\widetilde{H}_i} \right]_e\sum_{i \in CE(e)} n_{e,i}\widetilde{H}_e
$$

## 4 Design

### 4.1 Data types and parameters

The `PressureGradient` class will be used to perform the horizontal gradients of the pressure and geopotential
```c++
class PressureGrad{
    public:
    private:
        std::unique_ptr<PressureGrad> OmegaPressureGrad;
        PressureGradCentered CenteredPGrad;
        PressureGradHighOrder HighOrderPGrad; // To be implemented later
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
        static PressureGrad *init();
        static PressureGrad *get();
        void computePressureGrad();
    private:

}
```
#### 4.2.1 Creation

The constructor will be responsible for storing any static mesh information as private variables and handling the selection of the user-specified pressure gradient option.
```c++
PressureGrad::PressureGrad(const HorzMesh *Mesh, int NVertLevels, Config *Options);
```

#### 4.2.2 Initialization

The init method will create the pressure gradient and return an error code:
```c++
int PressureGrad::init();
```
It calls the constructor to create a new pressure gradient instance, producing a static managed (unique) pointer to the single instance.

#### 4.2.3 Retrieval

There will be a method for getting a pointer to the pressure gradient instance:
```c++
PressureGrad *PressureGrad::get();
```

#### 4.2.4 Computation

The public `computePressureGrad` method will rely on private methods for each specific pressure gradient option (centered and high order).
Note that the functors called by `computePressureGrad` are responsible for computing the sum of the pressure gradient and geopotential gradient accumulated in the `Tend` output array.
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
The clear method will remove the pressure gradient instance, calling the destructor in the process.
```c++
void PressureGrad::clear();
```

## Verification and Testing

### Test: Spatial convergence to exact solution
For a given analytical $v$,  $h$, and $b$, the spatial convergence of the pressure gradient can be assessed by computing errors on progressively finer meshes.

### Test: Baroclinic gyre
The baroclinic gyre test case will test the pressure gradient term in the full non-Boussinesq equations.
