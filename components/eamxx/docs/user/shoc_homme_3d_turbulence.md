# SHOC-HOMME 3D Turbulence

This page describes how to enable the 3D-turbulence coupling between SHOC and
HOMME, along with the main user-facing parameters associated with it.

## Enabling 3D Turbulence Coupling Between SHOC and HOMME

EAMxx exposes the 3D-turbulence option through the HOMME control-namelist
parameter `ctl_nl::do_3d_turbulence`. This is the setting that users should
change.

To enable the feature, run

``` {.shell .copy}
./atmchange ctl_nl::do_3d_turbulence=true
```

After changing the setting, rebuild the generated namelists before the next
run in the usual way, for example via `case.submit` or by re-running the case
setup/build workflow used for your case.

### What This Switch Does

When `ctl_nl::do_3d_turbulence=true`:

1. HOMME enables its 3D-turbulence path.
2. `eamxx_buildnml.py` mirrors that value into the locked EAMxx parameter
   `homme::do_3d_turbulence_homme`.
3. During atmosphere-driver initialization, that HOMME-facing flag is copied
   into SHOC's internal runtime option `do_3d_turbulence_shoc`.
4. SHOC computes horizontal eddy diffusivities for heat and momentum and passes
   them back to HOMME.
5. HOMME computes horizontal shear components and passes them to SHOC, which
   uses them to form the 3D shear-production term in the TKE equation.

Users should NOT edit `homme::do_3d_turbulence_homme` directly. That parameter
is intentionally locked and is maintained automatically from
`ctl_nl::do_3d_turbulence`.

### Related SHOC Parameters

The following SHOC parameters remain user-configurable and are relevant when
3D turbulence is enabled:

- `shoc::coeff_kh_horiz`: horizontal eddy-diffusivity coefficient for heat.
- `shoc::coeff_km_horiz`: horizontal eddy-diffusivity coefficient for momentum.
- `shoc::coeff_kh`: vertical eddy-diffusivity coefficient for heat.
- `shoc::coeff_km`: vertical eddy-diffusivity coefficient for momentum.

For example:

``` {.shell .copy}
./atmquery shoc::coeff_kh_horiz
./atmquery shoc::coeff_km_horiz
./atmchange shoc::coeff_kh_horiz=0.1
./atmchange shoc::coeff_km_horiz=0.1
```

The default values of `coeff_kh_horiz` and `coeff_km_horiz` are both `0.1`.
These parameters are tunable, but early testing suggests that `0.1` provides a
good balance for preserving isotropic and anisotropic turbulence behavior as
resolution changes.

A more exhaustive parameter study is still needed. At present, values much
larger than `1.0` appear to reduce the model's effective resolution
substantially and should generally be avoided.
