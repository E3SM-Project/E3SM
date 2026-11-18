# Vertical Mixing Coefficients

## 1 Overview

Representation of unresolved vertical fluxes of momentum, heat, salt, and biogeochemical tracers in ocean models is essential to simulations fidelity. Models of turbulent fluxes spans a wide range of complexity, but the models generally fall into a few categories: simply polynomial relationships, equilibrium turbulence models, and prognostic turbulence models.

## 2 Requirements

### 2.1 Requirement: Vertical mixing coefficient interface should be modular

To prepare for additional mixing closures and later improvements, the interface to vertical mixing must allow for easy connection of closures. It is assumed that the strength of vertical mixing from each closure is additive, such that the final vertical diffusion coefficient is the sum of the coefficients and potential non local terms from each closure.

### 2.2 Requirement: Vertical mixing models should be have local and/or gradient free terms

For optimal performance, initial models of vertical turbulent fluxes for Omega should be cast as an implicit, down gradient mixing and an explicit, gradient free, component, e.g.,

$$
\overline{w' \phi'}\approx \kappa(\overline{\rho},\overline{u},\overline{v}) \left(\frac{\partial \overline{\phi}}{\partial z} + \gamma(\overline{\rho},\overline{u},\overline{v}) \right)\,.
$$

Here, $\phi$ is a generic tracer, $\overline{w'\phi'}$ is the vertical turbulent flux of that tracer, $\kappa$ is the vertical diffusivity, and $\gamma$ is the gradient free portion of the flux (often referred to as 'non-local'). The vertical diffusivity $\kappa$ can be as simple as a constant value, to complex functions of quantities like shear and stratification. A similar equation can be written for turbulent momentum fluxes and vertical viscosity ($\nu$).

Invocation of the gradient diffusion hypothesis is a simplifying assumption for turbulence closures that casts turbulence as a purely diffusive process and allows the mixing problem to be cast implicitly in a tridagonal solve.  For more physical mixing closures, other methods can be explored, e.g., direct and iterative implicit solvers or subcycling.

The initial design and implementation of these vertical mixing coefficients will include only the local, down gradient portion of the flux (first right-hand-side term in Eq. (1)). Implementation of the non-local, gradient free portion of the vertical flux will come later. However, it is important to note here that $\gamma$ is only non-zero within the boundary layer, and when non-local flux formulations (such as [KPP (Large et al., 1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)) are being used, $\kappa$ value contributions such as the convective mixing flux detailed in Section 3.2 are only computed outside of the boundary layer and remain zero within the boundary layer.

### 2.3 Requirement: Vertical mixing models cannot ingest single columns at a time

Many standard vertical mixing libraries (e.g. CVMix) receive and operate on one column at a time. This results in poor performance, especially on accelerators. Any utilized model of vertical mixing must ingest and operate on multiple columns concurrently.

## 3 Algorithmic Formulation

For Omega-1 a few simple vertical mixing algorithms will be used. When computing $\kappa$ and $\nu$, linear superposition when combining the various parameterizations is assumed.

### 3.1 Richardson number dependent mixing

One of the simplest class of models of vertical turbulence are polynomial functions of the gradient Richardson number. As in MPAS-Ocean, we choose to define the Richardson number as

$$
Ri = \frac{N^2}{\left|\frac{\partial \mathbf{U}}{\partial z}\right|^2}\,,
$$

where $N^2$ is the Brunt Vaisala Frequency, which for non Boussinesq flows is defined as

$$
N^2 = \frac{g}{\rho_0}\frac{\partial \rho}{\partial z}\,.
$$

If we assume this term is defined at the top of the cell, it is discretized as

$$
N^2(k) = g \rho_0 \frac{\rho_{DD}(k) - \rho(k)}{z_m(k-1) - z_m(k)}\,.
$$

Here $\rho_{DD}$ is the density of the fluid of layer k-1 displaced to layer k adiabatically.

The shear in the denominator, also defined at the top of the cell, is discretized as

$$
\left| \frac{\partial \mathbf{U}}{\partial z}\right|^2 = \left(\frac{U_n(k-1)-U_n(k)}{z_m(k-1) - z_m(k)}\right)^2 + \left(\frac{U_t(k-1)-U_t(k)}{z_m(k-1) - z_m(k)}\right)^2\,.
$$

Where the subscripts *n* and *t* are the normal and tangential velocities respectively.

With the definition of the gradient richardson number, the viscosity and diffusivity are defined according to [Pacanowski and Philander (1981)](https://journals.ametsoc.org/view/journals/phoc/11/11/1520-0485_1981_011_1443_povmin_2_0_co_2.xml?tab_body=pdf) as

$$
\nu = \frac{\nu_o}{(1+\alpha Ri)^n} + \nu_b\,,
$$

$$
\kappa = \frac{\nu}{(1+\alpha Ri)} + \kappa_b\,.
$$

$\nu$ and $\kappa$ are defined at the top of the cell, similarly to $Ri$ and $N^2$. In these formulae, $\alpha$ and *n* are tunable parameters, most often chosen as 5 and 2 respectively. $\nu_b$ and $\kappa_b$ are the background viscosity and diffusivity respectively. MPAS-Ocean has an additional $\kappa_{b,passive}$ for applying to just the passive tracers, leaving the active tracers with $\kappa_b$.

### 3.2 Convective Instability Mixing

Commonly (and how this is currently handled by CVMix in MPAS-O), mixing due to convective instability is treated as a step function for the diffusivity and viscosity. This is often represented as

$$
\kappa =
\begin{cases}
\kappa_{b} + \kappa_{conv} \quad \text{ if } N^2 \leq N^2_{crit}\\
\kappa_{b} \quad \text{ if } N^2 > N^2_{crit}
\end{cases}
$$

A similar expression is utilized for viscosity. The effect of this formula is to homogenize T, S, passive tracers, and normal velocity, for unstable stratification, but also in neutral stratification, as $N^2_{crit}$ is typically chosen to be 0. The behavior in unstable stratification is likely correct, but for neutral stratification, homogenization of momentum is questionable. In Omega, we will slightly modify the algorithm above to,

$$
\kappa =
\begin{cases}
\kappa_{b} + \kappa_{conv} \quad \text{ if } N^2 < N^2_{crit}\\
\kappa_{b} \quad \text{ if } N^2 \geq N^2_{crit}
\end{cases}
$$

So that homogenization only occurs for unstable stratification (assuming $N^2_{crit} = 0$). Values for $\kappa_{conv}$ are typically large ($\sim 1 \, \mathrm{m^2/s}$). Again, when non-local, gradient free flux formulations are included in Eq. (1) in future developments, this convective instability contribution to the mixing will not be calculated within the boundary layer where $\gamma$ is active.

## 4 Design

Vertical chunking, as is, does not work well for vertical derivatives, so a first design of the vertical mixing coefficients computation does not do vertical chunking and just computes with the full-depth column.

Additionally, to start, only the down gradient contribution to vertical mixing will be added, with contributions to the vertical viscosity and diffusivity coming from the shear (Richardson number mixing), convective, and background mixing models detailed in the prior section. However, in future developments, the K Profile Parameterization [(KPP; Large et al., 1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872) and other non-local and/or higher-order mixing models will be added. Some of these parameterizations require vertical derivatives, full-depth integrals, and/or formulate $\kappa$ and $\gamma$ in Eq. (1) Section 2.1 as functions of surface forcing, which require either full-depth column or surface forcing information at depth, thus we design the mixing coefficients routine with this in mind. A solution that allows chunking with vertical operations such as derivatives and integrals and/or surface forcing values would be advantageous, but is outside of the scope of this design document.

### 4.1 Data types and parameters

Several parameters will be needed for the three different contributions to the local, down gradient vertical mixing coefficients. These parameters should be run-time configurable through the YAML config. A new section in the YAML config should be made for these parameters, as more run-time configurable vertical mixing parameters will be added to this section when the non-local formulations are added.

#### 4.1.1 Parameters

The following config options should be included for the Richardson number dependent -- Pacanowski and Philander (1981) based -- vertical mixing:

1. `EnableShearMix`. If true, shear-based mixing is computed based upon Pacanowski and Philander (1981) and applied to velocity components and tracer quantities. \[.true./.false.\]
2. `ShearNuZero`. Numerator of Pacanowski and Philander (1981) Eq (1). \[0.005 $m^2 s^{-1}$\]
3. `ShearAlpha`. Alpha value used in Pacanowski and Philander (1981) Eqs (1) and (2). \[5\]
4. `ShearExponent`. Exponent used in denominator of Pacanowski and Philander (1981) Eqs (1). \[2\]

The following config options should be included for convective vertical mixing:

1. `EnableConvectiveMix`. If true, convective mixing is computed and applied to velocity components and tracer quantities. \[.true./.false.\]
2. `ConvectiveDiffusivity`. Convective vertical viscosity and diffusivity applied to velocity components and tracer quantities, respectively. \[1.0 $m^2 s^{-1}$\]
3. `ConvectiveTriggerBVF`. Value of Brunt Viasala frequency squared below which convective mixing is triggered. \[0.0\]

The following config options should be included for background vertical mixing:

1. `BackgroundViscosity`. Background vertical viscosity applied to velocity components. \[1.0e-4 $m^2 s^{-1}$\]
2. `BackgroundDiffusivity`. Background vertical diffusivity applied to all tracer quantities. \[1.0e-5 $m^2 s^{-1}$\]

### 4.2 Methods

<!--- List and describe all public methods and their interfaces (actual code for
interface that would be in header file). Describe typical use cases. -->

It is assumed that viscosity and diffusivity will be stored at cell centers. Additionally, since they will be used elsewhere, it is assumed that zMid and N2 have been computed elsewhere (N2 in the density routine) and stored. Quantities such as squared shear and Richardson number are not needed outside of the vertical mixing routines, so are computed locally. Note: how topography will be dealt with has not been finalized yet, but will be adopted into the following once finalized.

```c++
parallelFor(
    {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int ICell, int K) {

        // Add background contribution to viscosity and diffusivity
        VertViscosity(ICell, K) = BackgroundViscosity;
        VertDiffusivity(ICell, K) = BackgroundDiffusivity;

        // If shear mixing true, add shear contribution to viscosity and diffusivity
        if (EnableShearMix) {
            const Real InvAreaCell = 1._Real / AreaCell(ICell);
            // Calculate the square of the shear
            for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                const Real Factor = 0.5_Real * DcEdge(JEdge) * DvEdge(JEdge) * InvAreaCell;
                const Real DelU2 = pow(NormalVelocity(JEdge,K-1) - NormalVelocity(JEdge,K), 2._Real) + pow(TangentialVelocity(JEdge,K-1) - TangentialVelocity(JEdge,K), 2._Real);
                ShearSquared(ICell, K) = ShearSquared(ICell, K) + Factor * DelU2;
            }
            ShearSquared(ICell, K) = ShearSquared(ICell, K) / pow(zMid(ICell, K-1) - zMid(ICell, K), 2._Real);

            // Calculate Richardson number
            RichardsonNum(ICell, K) = BruntVaisalaFreq(ICell, K) / max(ShearSquared(ICell, K), 1.0e-12_Real);

            // Add in shear contribution to vertical mixing
            VertViscosity(ICell, K) = VertViscosity(ICell, K) + ShearNuZero / pow(1._Real + ShearAlpha * RichardsonNum(ICell, K), ShearExponent);
            VertDiffusivity(ICell, K) = VertDiffusivity(ICell, K) + VertViscosity(ICell, K) / (1 + ShearAlpha * RichardsonNum(ICell, K));
        }

        // If conv mixing true, add conv contribution to viscosity and diffusivity
        if (EnableConvectiveMix) {
            if (BruntVaisalaFreq(ICell, K) < ConvectiveTriggerBVF) {
                VertViscosity(ICell, K) = VertViscosity(ICell, K) + ConvectiveDiffusivity;
                VertDiffusivity(ICell, K) = VertDiffusivity(ICell, K) + ConvectiveDiffusivity;
            }
        }

    });
```

A no-flux condition should be applied to both the vertical mixing of momentum and tracers through setting the vertical viscosity and diffusivity to zero at the boundaries. From here, the vertical viscosity and vertical diffusivity will enter into the [tridiagonal solver](./TridiagonalSolver.md) to compute the vertical flux of momentum and tracers.

## 5 Verification and Testing

Unit tests can be initialized with linear-with-depth initial conditions (velocity and density) and use a linear equation of state. Expected values of the Richardson number, vertical viscosity, and vertical diffusion coefficients can be computed and compared to. Tests for both the shear, convective, and combined mixing contributions should be made to test requirement (2.1).

This assumes that the displaced density has already been tested. Vertical fluxes of momentum and tracers will be tested separately with the [tridiagonal solver](./TridiagonalSolver.md).
