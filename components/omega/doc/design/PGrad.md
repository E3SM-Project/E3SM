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

For initial global cases without ice shelf cavities, the pressure and geopotential gradients will be computed with a simple centered difference approximation. In later versions of Omega, one or more high-order pressure gradients may be implemented to replace the centered approach in production runs.
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
T^p_{e,k} = - \left(\nabla \Phi \right)_{e,k} - \frac{1}{\left[\tilde{h}_k\right]_e} \nabla \left( \tilde{h}_k \alpha_k p_k \right) - \frac{1}{\left[\tilde{h}_k\right]_e} \left\{ \left[ \alpha p \nabla \tilde{z}\right]_{e,k}^\text{top} -  \left[ \alpha p \nabla \tilde{z}\right]_{e,k+1}^\text{top} \right\}.
$$

The geopotential and interface terms are necessary to account for tilted layers that occur when using a general vertical coordinate, where the gradient operator is taken along layers.
Even for layers at constant pseudo-height (rather than geometric height), the geoponential gradient would still be non-zero as will the second term (unless density is constant along layers).
In this equation, $\alpha_{i,k}$ is the specific volume, $p_{i,k}$ is the pressure, and $\Phi_{i,k}$ is the geopotential.
These three quantities are evaluated at the mid-point of layer $k$ of cell $i$ in the first two terms and at the cell interfaces in the third term along with the interface pseudo-height, $\tilde{z}$.
The discrete gradient operator at an edge is:

$$
 \nabla {(\cdot)} = \frac{1}{d_e} \sum_{i\in CE(e)} -n_{e,i} (\cdot)_i,
$$

where $d_e$ is the distance between cell centers, $CE(e)$ are the cells on edge $e$, and $n_{e,i}$ is the sign of the edge normal with respect to cell $i$.
The (cell-to-edge) horizontal averaging operator is:

$$
 [\cdot]_e = \frac{1}{2}\sum_{i\in CE(e)} (\cdot)_i.
$$

The $\gamma_{e,k-1/2} \equiv \left[ \alpha p \nabla \tilde{z}\right]_{e,k}^\text{top}/\left[\tilde{h}_k\right]_e$ term is computed using a four point average for all interior layer interfaces with a two-point averaged used for the top and bottom interfaces.
The interior interfaces are calculated as:
$$
\gamma_{e,k-1/2} = \frac{\sum_{i\in CE(e)} \alpha_{i,k} p_{i,k} h_{i,k} + \alpha_{i,k-1} p_{i,k-1} h_{i,k-1}}{d_e\sum_{i\in CE(e)} h_{i,k} + h_{i,k-1}} \sum_{i\in CE(e)} n_{e,i}\tilde{z}_{i,k-1/2} 
$$
with the top and bottom computed as:
$$
\gamma_{e,1/2} = \frac{\sum_{i\in CE(e)} \alpha_{i,k} p_{i,1} h_{i,1}}{d_e\sum_{i\in CE(e)} h_{i,1}} \sum_{i\in CE(e)} n_{e,i}\tilde{z}_{i,1/2}, \\
\gamma_{e,K+1/2} = \frac{\sum_{i\in CE(e)} \alpha_{i,K} p_{i,k} h_{i,K}}{d_e\sum_{i\in CE(e)} h_{i,K}} \sum_{i\in CE(e)} n_{e,i}\tilde{z}_{i,K+1/2}. 
$$

Therefore, the centered pressure gradient will be calculated as:

$$
 T^p_{e,k} = \frac{1}{d_e}\sum_{i\in CE(e)} n_{e,i} \Phi_{i,k} + \frac{2}{d_e\sum_{i \in CE(e)}\tilde{h}_{i,k}}\sum_{i \in CE(e)} n_{e,i} \tilde{h}_{i,k}\alpha_{i,k}p_{i,k} - \gamma_{e,k-1/2} + \gamma_{e,k+1/2}.
$$


### 3.2 High-order Pressure Gradient
The high order pressure gradient will be based on the {ref}`full volume integral form <omega-v1-vh-momentum-reynolds2>` of the geopotential and pressure terms:
$$
T^p &= - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, \left( \nabla \left<\Phi\right> \right) \, d\tilde{z} \, dA \\
& - \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \left(\left< \alpha \right> \left<p \right> + \left<\alpha^\prime p^\prime\right> \right) \, d\tilde{z} \right) dl \\
& - \int_A \rho_0 \left[ \left< \alpha \right> \left<p \nabla \tilde{z}_k^{\text{top}} \right> + \left<\alpha^\prime \left(p \nabla \tilde{z}_k^{\text{top}}\right)^\prime\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& + \int_A \rho_0 \left[ \left< \alpha \right> \left<p \nabla \tilde{z}_k^{\text{bot}} \right> + \left<\alpha^\prime \left(p \nabla \tilde{z}_k^{\text{bot}}\right)^\prime\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA.
$$
To obtain the expression that will be used, we neglect the turbulent correlations and drop the Reynold's average notation for single variables:
$$
T^p &= - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, \left( \nabla \Phi \right) \, d\tilde{z} \, dA \\
& - \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \left(\alpha p \right) \, d\tilde{z} \right) dl \\
& - \int_A \rho_0 \left[ \alpha \left<p \nabla \tilde{z}_k^{\text{top}} \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& + \int_A \rho_0 \left[ \alpha \left<p \nabla \tilde{z}_k^{\text{bot}} \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA.
$$
These volume and area integrals will be computed using quadrature to account for the variability of $\alpha$ with the reconstructed values of temperature, salinity, and pressure at the quadrature points.
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
        // Map of all pressure gradients
        static std::map<std::string, std::unique_ptr<PressureGrad>> AllPGrads;

        // Instances of functors
        PressureGradCentered CenteredPGrad;
        PressureGradHighOrder HighOrderPGrad1; // To be implemented later
        PressureGradHighOrder HighOrderPGrad2; // Multiple high order options are likely in the future

        // Pressure gradient choice from config
        PressureGradType PressureGradChoice;

        // Data required for computation (stored copies of HorzMesh/VCoord arrays)
        Array2DI4 CellsOnEdge;
        Array1DReal DvEdge;
        Array1DReal EdgeSignOnCell;
        Array1DI4 MinLayerEdgeBot;
        Array1DI4 MaxLaterEdgeTop;

        // Array for interface produce needed in centered pressure gradient
        Array2DReal InterfaceProduct;
}
```
The functions to compute the centered and high order pressure gradient terms will be implemented as functors and the pressure gradient class will have private instances of these classes.
The `PressureGrad` class will store pointers to the static mesh and vertical coordinate needed in the calculation of the pressure gradient, which will be initialized on construction.
The class will support multiple named instances to allow for the possibility of a barotropic pressure gradient calculation on a different mesh in the future.
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
### 4.2 Methods
The `PressureGrad` class will have methods for creating, retrieving, and removing instances similar to other Omega classes such as `Dcomp`, `HorzMesh`, etc.
Computation of the pressure gradient will be performed by `computePressureGrad`, which will select the user's chosen implementation from the available options.
It will also have a helper method to pre-compute the interface term, $\gamma_{e,k-1/2}$, needed for the centered pressure gradient to separate layer interface and midpoint calculations.
```c++
class PressureGrad{
    public:
        // Creation methods
        static PressureGrad *init();
        static PressureGrad *create(const std::string &Name, const HorzMesh *Mesh,
                                    const VertCoord *VCoord, Config *Options);

        // Retrival methods
        static PressureGrad *get(const std::string &Name);
        static PressureGrad *getDefault();

        // Removal methods
        static void clear();
        static void erase(const std::string &Name);
        ~PressureGrad();

        // Main compute method
        void computePressureGrad(Array2DReal &Tend, const OceanState *State,
                                 const VertCoord *VCoord, const Eos *EqState,
                                 const int TimeLevel);
    private:
        // Helper compute method
        void computeInterfaceProduct(const Array2DReal &PressureMid,
                                     const Array2DReal &SpecVol,
                                     const Array2DReal &LayerThick,
                                     const Array2DReal &ZInterface);

}
```
#### 4.2.1 Creation

The `create` method will be responsible for calling the constructor which will, store any static mesh information as private variables and handle the retrieval of the user-specified pressure gradient option.
```c++
PressureGrad *PressureGrad::create(const std::string &Name, 
                                   const HorzMesh *Mesh,
                                   const VertCoord *VCoord,
                                   Config *Options);
```

#### 4.2.2 Initialization

The init method will create the default pressure gradient and return an error code:
```c++
static PressureGrad::init();
```
It calls the constructor to create a new pressure gradient instance, producing a static managed (unique) pointer to the single instance.

#### 4.2.3 Retrieval

There will be a method for getting a pointer to the default pressure gradient instance:
```c++
PressureGrad *PressureGrad::getDefault();
```
or a named instance:
```c++
PressureGrad *PressureGrad::get(const std::string &Name);
```

#### 4.2.4 Computation

The public `computePressureGrad` method will rely on private methods for each specific pressure gradient option (centered and high order).
Note that the functors called by `computePressureGrad` are responsible for computing the sum of the pressure gradient and geopotential gradient accumulated in the `Tend` output array.
```c++
void PressureGrad::computePressureGrad(const Array2DReal &Tend,
                                       const OceanState *State,
                                       const VertCoord *VCoord,
                                       const Eos *EqState,
                                       const I4 TimeLevel) {

   OMEGA_SCOPE(LocCenteredPGrad, CenteredPGrad);
   OMEGA_SCOPE(LocHighOrderPGrad, HighOrderPGrad);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);

   const Array2DReal &PressureMid       = VCoord->PressureMid;
   const Array2DReal &PressureInterface = VCoord->PressureInterface;
   const Array2DReal &Geopotential      = VCoord->GeopotentialMid;
   const Array2DReal &SpecVol           = EqState->SpecVol;
   const Array2DReal &ZInterface        = VCoord->ZInterface;
   Array2DReal LayerThick;
   State->getLayerThickness(LayerThick, TimeLevel);

   if (PressureGradChoice == PressureGradType::Centered) {

      // computes alpha*p*grad(z) term at edge interfaces
      computeInterfaceProduct(PressureMid, SpecVol,
                              LayerThick, ZInterface);

      // computes centered geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-centered", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin = LocMinLayerEdgeBot(IEdge);
             const int KMax = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRange(KMin, KMax);

             parallelForInner(
                 Team, KRange,
                 INNER_LAMBDA(int KChunk) {
                    LocCenteredPGrad(Tend, IEdge, KChunk, PressureMid,
                                     Geopotential, LayerThick,
                                     LocInterfaceProduct, SpecVol);
                 });
          });

   } else {
      
      // computes high-order geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-highorder", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin = MinLayerEdgeBot(IEdge);
             const int KMax = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRange(KMin, KMax);

             parallelForInner(
                 Team, KRange,
                 INNER_LAMBDA(int KChunk) {
                    LocHighOrderPGrad(Tend, IEdge, KChunk, PressureMid,
                                      Geopotential, SpecVol);
                 });
          });

   }
} // end compute pressure gradient
```

#### 4.2.5 Destruction and removal

No operations are needed in the destructor.
The there will be a  method to remove all pressure gradient instances,
```c++
static void PressureGrad::clear();
```
or a named instance
```c++
static void erase(const std::string &Name);
```


## Verification and Testing

### Test: Hydrostatic test with tilted layers
The pressure gradient term will be computed on an edge with two adjacent cells that have the same sea surface height, but internal layers that are tilted in terms of geometric height.
The same linear temperature and salinity profiles will be used for both cells, thus the pressure gradient should be zero on the edge. 
The pseudo-height (and specific volume) of the layers will be computed via iteration, where the pressure is updated in the EOS evaluation until convergence is reached.
The pressure gradient term will be compared to the zero exact solution at progressively finer resolutions to evaluate convergence.

### Test: Spatial convergence to exact solution
For given analytical functions of $\alpha$, $h$, and $z$, the spatial convergence of the pressure gradient can be assessed by computing errors on progressively finer meshes.
This will be compared to an exact solution computed using a high-order quadrature.

### Test: Seamount test
The seamount test in Polaris will be used to verify the pressure gradient's ability to preserve the resting state of fluid in a case with tilted layers.

### Test: Baroclinic gyre
The baroclinic gyre test case will test the pressure gradient term in the full non-Boussinesq equations.
