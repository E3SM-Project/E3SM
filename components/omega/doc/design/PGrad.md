(omega-design-pressure-grad)=
# Pressure Gradient


## 1 Overview
The pressure gradient will be responsible for computing the horizontal gradients of both the pressure and geopotential terms for the non-Boussinesq primitive equations implemented in Omega.
In the non-Boussinesq model, the conserved quantity is mass rather than volume.
In Omega the prognostic variable $\tilde{h}$ is a pseudo thickness, rather than geometric thickness as in a Boussinesq model.
 Some non-Boussinesq models are written in pressure coordinates (e.g. [de Szoeke and Samelson 2002](https://journals.ametsoc.org/view/journals/phoc/32/7/1520-0485_2002_032_2194_tdbtba_2.0.co_2.xml).
 However, Omega is written in general vertical coordinates and can reference either pressure $p$ or distance $z$ in the vertical.
In a pure pressure coordinate the pressure gradient term disappears (since the pressure does not vary along lines of constant pressure), just as how the geopotential term disappears in a pure z coordinate model.
However, similar to MPAS-Ocean's support for tilted height coordinates, Omega will allow for tilted pressure coordinates.
This means that Omega will need to compute both the pressure and geopotential gradients.

## 2 Requirements

### 2.1 Requirement: Support for tilted pressure coordinates for the non-Boussinesq primitive equations

The pressure gradient will compute the horizontal gradients of both the pressure and geopotential to support tilted pressure coordinates in the non-Boussinesq model.
This will allow for the use of a $p^\star$ coordinate, which functions similarly to the $z^\star$ in the Boussinesq MPAS-Ocean model.
Note that we use the name "$p^\star$" to refer to the vertical coordinate, in Omega it will be expressed in terms of the pseudo-height, $\tilde{z}$, as opposed to pressure directly.

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
This requirement is satisfied by tidal forcing terms being included in the geopotential calculation in the `VertCoord` class.

### 2.5 Desired: Pressure gradient for barotropic mode

For split barotropic-baroclinic timestepping, the pressure gradient should provide the bottom pressure gradient tendency in the barotropic mode.
The details of the barotropic pressure gradient will be added in a future design document for split time stepping.

## 3 Algorithmic Formulation
### 3.1 Centered Pressure Gradient
In the layered non-Boussinesq  {ref}`momentum equation <omega-v1-momentum-eq>` solved in Omega, the pressure gradient tendency term for edge $e$ and layer $k$, $T^p_{e,k}$, includes the gradient of the geopotential, the gradient of a term involving pressure, and two terms evaluated at the cell interface:

$$
T^p_{e,k} = - \left(\nabla \Phi \right)_{e,k} - \frac{1}{\left[\tilde{h}_k\right]_e} \nabla \left( \tilde{h}_k \alpha_k p_k \right) + \frac{1}{\left[\tilde{h}_k\right]_e} \left\{ \left[ \alpha p \nabla \tilde{z}\right]_{e,k}^\text{top} -  \left[ \alpha p \nabla \tilde{z}\right]_{e,k+1}^\text{top} \right\}.
$$

The geopotential and interface terms are necessary to account for tilted layers that occur when using a general vertical coordinate, where the gradient operator is taken along layers.
In this equation, $\alpha_{i,k}$ specific volume, $p_{i,k}$ is the pressure, and $\Phi_{i,k}$ is the geopotential.
These three quantities are evaluated at the mid-point of layer $k$ of cell $i$ in the first two terms and at the cell interfaces in the third term along with the interface psudo-height, \tilde{z}.
The discrete gradient operator at an edge is:

$$
 \nabla {(\cdot)} = \frac{1}{d_e} \sum_{i\in CE(e)} -n_{e,i} (\cdot)_i
$$

where $d_e$ is the distance between cell centers, $CE(e)$ are the cells on edge $e$, and $n_{e,i}$ is the sign of the edge normal with respect to cell $i$.
The (cell-to-edge) horizontal averaging operator is:

$$
 [\cdot]_e = \frac{1}{2}\sum_{i\in CE(e)} (\cdot)_i
$$

Therefore, the centered pressure gradient will be calculated as:

$$
 T^p_{e,k} = \frac{1}{d_e}\sum_{i\in CE(e)} n_{e,i} \Phi_{i,k} + \frac{2}{d_e\sum_{i \in CE(e)}\tilde{h}_{i,k}}\left(\sum_{i \in CE(e)} n_{e,i} \tilde{h}_{i,k}\alpha_{i,k}p_{i,k}\right. \\ \left. - \frac{1}{2} \sum_{i\in CE(e)} \alpha_{i,k-1/2} p_{i,k-1/2}\sum_{i\in CE(e)} n_{e,i}\tilde{z}_{i,k-1/2} \right.\\ \left.+ \frac{1}{2} \sum_{i\in CE(e)} \alpha_{i,k+1/2} p_{i,k+1/2}\sum_{i\in CE(e)} n_{e,i}\tilde{z}_{i,k+1/2}\right),
$$
where $k-1/2$ and $k+1/2$ refer to the top and bottom of layer $k$, respectively.


### 3.2 High-order Pressure Gradient
The high order pressure gradient will be based on the {ref}`full volume integral form <omega-v1-vh-momentum-reynolds2>` of the geopotential and pressure terms:
$$
T^p &= - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, \left( \nabla \left<\Phi\right> \right) \, d\tilde{z} \, dA \\
& - \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \left(\left< \alpha \right> \left<p \right> + \left<\alpha^\prime p^\prime\right> \right) \, d\tilde{z} \right) dl \\
& - \int_A \rho_0 \left[ \left< \alpha \right> \left<p \nabla \tilde{z}_k^{\text{top}} \right> + \left<\alpha^\prime \left(p \nabla \tilde{z}_k^{\text{top}}\right)^\prime\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& + \int_A \rho_0 \left[ \left< \alpha \right> \left<p \nabla \tilde{z}_k^{\text{bot}} \right> + \left<\alpha^\prime \left(p \nabla \tilde{z}_k^{\text{bot}}\right)^\prime\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA.
$$
To obtain the expression that will be used, we neglect the turblent correlations and drop the Reynold's average notation for single variables:
$$
T^p &= - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, \left( \nabla \Phi \right) \, d\tilde{z} \, dA \\
& - \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \left(\alpha p \right) \, d\tilde{z} \right) dl \\
& - \int_A \rho_0 \left[ \alpha \left<p \nabla \tilde{z}_k^{\text{top}} \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& + \int_A \rho_0 \left[ \alpha \left<p \nabla \tilde{z}_k^{\text{bot}} \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA.
$$
These volume and area integrals will be computed using quadrature to account for the variablity of $\alpha$ with the recontructed values of temperature, salinity, and pressure atthe quadrature points.
The complete details for the high-order pressure gradient will be the subject of a future design document.

%### 3.3 Barotropic Pressure Gradient
%
%When split baroclinic-barotropic time stepping is implemented in the future, the barotropic pressure gradient will be calculated by the pressure gradient class.
%The barotropic pressure gradient is found by depth integrating the pressure gradient.
%The pressure is
%
%$$
%p(z) = p_b - g \int^z_{-h} \rho dz^\prime,
%$$
%
%where $p_b$ is the bottom pressure.
%The bottom pressure is the sum of the atmospheric surface pressure, $p_s$, and the pressure contribution of the water column:
%
%$$
%p(z) &=  p_s + g\int_{-h}^\eta \rho dz - g \int^z_{-h} \rho dz^\prime, \\
%     &=  p_s + g\rho_0\widetilde{H} - g \int^z_{-h} \rho dz^\prime,
%$$
%
%where the total water column pseudo height is expressed by
%
%$$
%\widetilde{H} = \int_{-h}^\eta \frac{\rho}{\rho_0} dz.
%$$
%
%$\widetilde{H}$ is the prognositc variable in the barotropic continuity equation.
%The vertical integral of the pressure gradient is
%
%$$
%\frac{1}{\rho_0\widetilde{H}}\int^\eta_{-h} \nabla p dz &= \frac{1}{\rho_0\widetilde{H}}\int^\eta_{-h} \nabla \left( p_s + g\rho_0 \widetilde{H} - g \int^z_{-h} \rho dz^\prime \right) dz, \\
%                           &= \frac{H}{\rho_0\widetilde{H}}\nabla p_s + \frac{gH}{\widetilde{H}}\nabla \widetilde{H} - \frac{g}{\rho_0\widetilde{H}} \int_{-h}^\eta \left( \nabla \int_{-h}^z \rho dz^\prime\right) dz,
%$$
%
%where the height of the water column is represented by $H$.
%The $1/\rho_0\widetilde{H}$ factor comes vertically integrating the material derivative and expressing the resulting barotropic momentum equation in non-conservative form.
%
%Therefore, the barotorpic pressure gradient term is discretized as:
%
%$$
%\overline{T}_e^p = g\left[ \frac{H_i}{\widetilde{H}_i} \right]_e\sum_{i \in CE(e)} n_{e,i}\widetilde{H}_e
%$$

## 4 Design

### 4.1 Data types and parameters

The `PressureGradient` class will be used to perform the horizontal gradients of the pressure and geopotential
```c++
class PressureGrad{
    public:
    private:
        std::unique_ptr<PressureGrad> OmegaPressureGrad;
        PressureGradCentered CenteredPGrad;
        PressureGradHighOrder HighOrderPGrad1; // To be implemented later
        PressureGradHighOrder HighOrderPGrad2; // Multiple high order options are likely in the future
        PressureGradType PressureGradChoice;
        I4 NVertLevels;
        I4 NChunks;
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
   HighOrder1,
   HighOrder2
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
PressureGrad::PressureGrad(const HorzMesh *Mesh, const VertCoord *VCoord, Config *Options);
```

#### 4.2.2 Initialization

The init method will create the default pressure gradient and return an error code:
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
const Array1DI4 &MinLyrEdgeBot   = VCoord->MinLayerEdgeBot;
const Array1DI4 &MaxLyrEdgeTop   = VCoord->MaxLayerEdgeTop;
    if (PressureGradChoice == PressureGradType::Centered){
       parallelForOuter(
          {NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int NChunks = computeNChunks(MinLyrEdgeBottom, MaxLyrEdgeTop, IEdge);
             parallelForInner(Team, NChunks, [=](const int KChunk) {
                LocCenteredPGrad(Tend, IEdge, KChunk, Pressure, Geopotential, SpecVol);
             });
       });
    }
    else if (PressureGradChoice == PressureGradType::HighOrder){
       parallelForOuter(
          {NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int NChunks = computeNChunks(MinLyrEdgeBottom, MaxLyrEdgeTop, IEdge);
             parallelForInner(Team, NChunks, [=](const int KChunk) {
                LocHighOrderPGrad(Tend, IEdge, KChunk, Pressure, Geopotential, SpecVol);
             });
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
For given analytical functions of $\alpha$, $h$, and $z$, the spatial convergence of the pressure gradient can be assessed by computing errors on progressively finer meshes.

### Test: Baroclinic gyre
The baroclinic gyre test case will test the pressure gradient term in the full non-Boussinesq equations.
