(omega-design-eos)=

# Equation of State (EOS)

## 1 Overview
The equation of state relates density to the prognostic state variables temperature and salinity. In the case of a non-Boussinesq model, it is also dependent on prognostic pressure. The prognostic temperature, $\theta$, is either potential temperature or Conservative Temperature, depending on the chosen equation of state, and  $S$ is the absolute salinity (linear eos could use practical salinity). Given that the Omega governing equations are  non-Boussinesq, the equation of state class will provide specific volume from the state variables. It will provide methods for computing the specific volume, and its first derivatives.

## 2 Requirements

### 2.1 Requirement: calculate specific volume and density from prognostic temperature and salinity
Omega will use TEOS-10 as the primary equation of state. The reference implementation of TEOS-10 is contained in the Gibbs Sea Water (GSW) toolbox function. The Omega implementation will be equivalent to the computation in the GSW toolbox function: `gsw_specvol(sa,ct,p)`

### 2.2 Requirement: allow a choice of equations of state
Omega-1 should include the option for two EOS options. The first is a linear eos, for benchmarking against MPAS-O and to use in simplified cases. The other option should be a well-established version of TEOS-10 (e.g. Roquet 75-term polynomial).

### 2.3 Requirement: computational efficiency
As one of the most expensive kernels, the EOS must be as efficient as possible. Developers should be discouraged from calling the EOS too frequently. The fraction of compute time used be the EOS within Omega1 will be compared to fraction used in MPAS-Ocean.

### 2.4 Requirement: abide by license for the GSW Oceanographic Toolbox
The GSW toolbox can only be used without modification, so a Omega-specific port of certain functions will be required to support Kokkos/GPU enabled computation. The toolbox without modifications can still be used for verification testing against the Omega implementation.

### 2.5 Requirement: verification
The Omega implementation should be compared with the test value published in Roquet et al 2015, and the official GSW Oceanographic Toolbox (supported by TEOS-10) to verify correctness.

### 2.6 Requirement: check value range
The EOS code should include a check on the range of ocean properties to assess use of TEOS-10 (equivalent to the ``ocean funnel'' check), and a truncation of state variables to use the polynomial only in the valid range (similar to the truncation in MPAS-Ocean).

### 2.7 Requirement: provide adiabatically-displaced density
The EOS code should be able to be used to calculate density adiabatically displaced to i) the layer below, ii) the surface. This is needed to calculate the Brunt-Vaisala and for mixed-layer parameterizations. This operation should be done as efficiently as possible (e.g. we will explore reusing parts of the polynomial or leveraging Kernel Fusion).

### 2.8 Desired: freezing temperature calculation
Later versions should include calculating the freezing temperature of seawater.

### 2.9 Desired: output in situ temperature
Later versions should include calculating the in situ  temperature needed for coupling (at surface including at non-zero pressure).

### 2.10 Desired: output first derivatives
The EOS will need to be able to produce alpha and beta (drho/dT and drho/dS) for linear expansions in the higher-order pressure gradient and for some mixing parameterizations.

### 2.11 Desired: variable conversion
Later versions should include functions to convert between conservative and potential temperatures, and absolute and practical salinity. These will be used to convert the initial conditions or compare to ocean states from other sources (MPAS-Ocean, reanalysis, other models etc.) These functions may be used offline in pre- or post-processing but need to be consistent with the EOS implementation.

### 2.12 Desired: enthalpy conversion
Later versions should include functions to calculate the layer potential enthalpy from enthalpy fluxes to convert coupling fluxes to changes in conservative temperature in a rigorous manner.

## 3 Algorithmic Formulation
The formulation of the specific volume calculation will be the 75-term polynomial as documented in [Roquet et al. 2015](https://www.sciencedirect.com/science/article/pii/S1463500315000566). The implementation of the displaced density calculation may be altered to improve performance by reusing some polynomial coefficients, but should produce the same results as the full 75-term polynomial.

## 4 Design
The GSW-C toolbox will be incorporated into Omega as a submodule. In order to abide by the license of the toolbox, the toolbox will only be used in testing to check our implementation. The coefficients for the Roquet et al. 2015 expansion will be based on the published values (Appendix). A Kokkos implementation of the function to compute specific volume will be ported to Omega. The GSW-C toolbox submodule will serve as a baseline reference for our ports in unit tests.

To optimize performance, calls to the EOS function should be able to support whole array computation. The EOS class will include a public method that calculate specific volume given the state variables, and a public variable defining the chosen eos option (changeable by the user).

### 4.1 Data types and parameters
The `Eos` class will be used to perform the necessary eos operations such as calculating specific volume:
```c++
class Eos{
    public:
       Array SpecVol;
       Array SpecVolDisplaced;
       void computeSpecVol();
       static Eos *create();
    private:
       EOSType eosChoice;
       I4 NCellsAll;
       I4 NChunks = NVertLayers / VecLength;
       void computeSpecVolTEOS10Poly75t();
       void computeSpecVolLinear();
       void truncateTempSal();
       static Eos *DefaultEos;
       static std::map<std::string, std::unique_ptr<Eos>> AllEos;

```
The user will select the EOS type at runtime in the input yaml file under an eos section:
```yaml
    eos:
       eosType: 'teos10'
```

An `enum class` will be used to specify options for the EOS used for an Omega simulation:
```c++
enum class EOSType{
   Linear,  /// Linear equation of state
   TEOS10Poly75t  /// Roquet et al. 2015 75 term expansion
}
```

We assume that the call to the eos is in a timestepping routine where input tracers for this timestep are defined, eg.
```c++
   Array2DReal ConservativeTemperature;
   Array2DReal AbsoluteSalinity;
   Tracers->getByName(ConservativeTemperature, VelTimeLevel, "temperature");
   Tracers->getByName(AbsoluteSalinity, VelTimeLevel, "salinity");
```


### 4.2 Methods
There will be a constructor and destructor for the class, as well as several public and private
methods. Modeled on the `Tendencies` class, the constructor will be private. A static `create` method is used to ensure every eos instance is properly stored in the static map of eos instances. This is to help ensure that the eos initialized in the `init` step of the model set-up is then retrievable during the `forward` step.


#### 4.2.1 Creation

The constructor will be responsible for:
  * determining the selected eos option
  * allocating arrays

```c++
Eos::Eos(const HorzMesh *Mesh, int NVertLayers, Config *Options);
```

The create method will take the same arguments as the constructor plus a name. It calls the constructor to
create a new eos instance, and put it in the static map of all eos.
It will return a pointer to the newly created object.
```c++
Eos *Eos::create(const std::string &Name, const HorzMesh *Mesh, int NVertLayers, Config *Options);
```

#### 4.2.2 Initialization

The init method will create the default eos and return an error code:
```c++
int Eos::init();
```

#### 4.2.3 Retrieval

There will be methods for getting the default and non-default eos instances:
```c++
Eos *Eos::getDefault();
Eos *Eos::get(const std::string &Name);
```

#### 4.2.4 Computation

The public `computeSpecVol` method will rely on private methods for each specific EOS option (linear and TEOS-10).
```c++
void Eos::computeSpecVol(const Array2DReal &SpecVol,
                         const Array2DReal &ConservativeTemperature,
                         const Array2DReal &AbsoluteSalinity,
                         const Array2DReal &Pressure) {
OMEGA_SCOPE(...)
    if (eosChoice == EOSType::Linear){
       parallelFor("eos-linear", {NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocComputeSpecVolLinear(...);
       });
    }
    else if (eosChoice == EOSType::TEOS10Poly75t){
       parallelFor("eos-teos10", {NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocComputeSpecVolTEOS10Poly75t(...);
        });
    }...
```



#### 4.2.5 Destruction and removal
No operations are needed in the destructor as the Kokkos arrays are removed
when they are no longer in scope. The erase method
will remove a named eos instance, whereas the clear method will remove all of
them. Both will call the destructor in the process.
```c++
void Eos::erase(const std::string &Name);
void Eos::clear();
```
## 5 Verification and Testing

### 5.1 Test: Verification of GSW-C submodule
A unit test will verify that the result of the GSW-C toolbox (used as a submodule, without modifications) matches the expected value of specific volume published in Roquet et al. 2015 within machine precision. The publication only includes one data point (Sa, Ct, P) to check.

### 5.2 Test: Verification of our TEOS-10 75-term polynomial implementation
A unit test will verify that the result of the Omega implementation of the Roquet et al. 2015 75-term polynomial calculation of the specific volume matches the expected value of specific volume published in Roquet et al. 2015 within machine precision. The publication only includes one data point (Sa, Ct, P) to check.

### 5.3 Test: Verification of TEOS-10 with GSW-C toolbox
A unit test will verify that the result of the Omega and GSW-C toolbox implementations for the Roquet et al. 2015 expansion are within machine precision over a range of conservative temperature, absolute salinity, and pressure ranges.

### 5.4 Test: Verification of linear EOS
The linear EOS should be close, but not equal to the TEOS-10 result. The calculated density should be exactly linear in T,S. It will be verified against known values from the MPAS-Ocean model for a range of temperature and salinity values.

### 5.5 Test: conservation of potential enthalpy
The total potential enthalpy, or equivalently, the conservative temperature should be conserved in the absence of external energy fluxes. This could be based on the tracer conservation test (e.g. merry-go-round), applied to active tracers.
