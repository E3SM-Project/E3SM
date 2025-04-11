# Full-model Testing

Full model system testing of EAMxx is done via CIME test cases
(much like the rest of E3SM).

We offer a number of test suites, including:

- `e3sm_eamxx_v1_lowres`
      - A handful of reasonably quickly-running tests run at
        ultra-low-resolution (`ne4` grid).
- `e3sm_eamxx_mam4xx_v1_lowres`
      - Tests designed for the MAM4xx aerosol library and are separate from
        the above (mostly `ne4`, though also includes `ne30`).
- `e3sm_eamxx_v1_dp-eamxx`[^dpxx]
      - A quickly-running test that employs the doubly-periodic (dp)
        configuration of EAMxx at low-resolution (`ne30)`.
- `e3sm_eamxx_v1_medres`
      - Tests intended to be a middle-ground (`ne30`) between the
        `_lowres` and `_hires` in that it is a more comprehensive test of
        capabilities than low-resolution but runs more quickly than high-resolution.
- `e3sm_eamxx_v1_hires`
      - A small number of larger, longer-running (high-resolution) tests
        to measure performance.

***Note:*** See the collapsed section "Lists of Test Cases and Suites..." below
for more information on these test suite configurations.

```{ .shell .copy title="Running a Test Suite"}
cd $repo/cime/scripts
./create_test e3sm_scream_v1_at --wait
```

```{ .shell .copy title="Running a Single Test Case"}
cd $repo/cime/scripts
./create_test SMS.ne4_ne4.F2010-SCREAMv1 --wait
```

- There are many behavioral tweaks you can make to a test case, like
changing the run length, test type, etc.
      - Most of this is not specific to EAMxx and works for any CIME case.
- This general information and much more can be found in the [CIME docs](http://esmci.github.io/cime/versions/master/html/users_guide/testing.html).

## EAMxx-specific Model Configuration

The main model-level configuration options for EAMxx are:

1. ***grids***
1. ***compsets***
1. ***testmods***

### Common EAMxx Grids[^grid_nomenclature]

- `ne4_ne4`:[^ne4-testing-only] ultra-low-resolution
- `ne4pg2_ne4pg2`:[^ne4-testing-only] ultra-low-resolution with phys grid
- `ne30_ne30`: low-resolution
- `ne30pg2_ne30pg2`: low-resolution with phys grid
- `ne1024pg2_ne1024pg2`: ultra-high-resolution with phys grid

More information about grids may be found on this [Confluence Page](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/933986549/ATM+Grid+Resolution+Summary).

### Common EAMxx Compsets

- `F2010-SCREAMv1`: V1 standard compset with EAMxx V1 atmosphere
- `FIOP-SCREAMv1-DP`: V1 with `dpxx`[^dpxx]
- `F2010-SCREAMv1-noAero`: V1 without aerosol forcing

Full info on supported compsets may be found by taking a look at the EAMxx
[`config_compsets.xml`](https://github.com/E3SM-Project/E3SM/blob/master/components/eamxx/cime_config/config_compsets.xml)
source file.

### Common EAMxx Testmods

- `small_kernels`
      - Enable smaller-granularity kernels, which can improve performance
      on some systems.
- `scream-output-preset-{i}, i = 1,..., 6`
      - 6 output presets for EAMxx
      - These turn some combination of our three output streams
      (`phys_dyn`, `phys`, and `diags`), various remaps, etc.
- `bfbhash`
      - Turns on bit-for-bit hash output.
      - See [this Confluence page](https://acme-climate.atlassian.net/wiki/spaces/NGDNA/pages/3831923056/EAMxx+BFB+hashing)
      for additional information.

!!! Abstract "Further Reading"
    Additional information about running EAMxx may be found at
    [this webpage](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3386015745/How+To+Run+EAMxx+SCREAMv1).

<!-- false-positive spaces in code fence marker -->
<!-- markdownlint-disable MD038 -->
??? Info "List of Selected Test Cases and Suites from `cime_config/tests.py`"
    ```{.python .copy}
    "e3sm_eamxx_v1_lowres" : {
            "time"  : "01:00:00",
            "inherit" : ("e3sm_eamxx_mam4xx_v1_lowres"),
            "tests" : (
                "ERP_D_Lh4.ne4_ne4.F2010-SCREAMv1.eamxx-output-preset-1",
                "ERS_Ln9.ne4_ne4.F2000-SCREAMv1-AQP1.eamxx-output-preset-2",
                "SMS_D_Ln9.ne4_ne4.F2010-SCREAMv1-noAero.eamxx-output-preset-3",
                "ERP_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.eamxx-output-preset-4",
                "ERS_D_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.eamxx-rad_frequency_2--eamxx-output-preset-5",
                "ERS_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.eamxx-small_kernels--eamxx-output-preset-5",
                "ERS_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.eamxx-small_kernels_p3--eamxx-output-preset-5",
                "ERS_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.eamxx-small_kernels_shoc--eamxx-output-preset-5",
                "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.eamxx-mam4xx-all_mam4xx_procs",
                )
        },
    "e3sm_eamxx_mam4xx_v1_lowres" : {
            "time"  : "01:00:00",
            "tests" : (
                "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.eamxx-mam4xx-optics",
                "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.eamxx-mam4xx-aci",
                "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.eamxx-mam4xx-wetscav",
                "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.eamxx-mam4xx-drydep",
                "SMS_D_Ln5.ne30pg2_oECv3.F2010-SCREAMv1-MPASSI.eamxx-mam4xx-remap_emiss_ne4_ne30"
            )
        },
    "e3sm_eamxx_v1_dp-eamxx" : {
        "time"  : "01:00:00",
        # each test runs with 225 dynamics and 100 physics columns,
        # roughly size of ne2
        "tests" : (
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FIOP-SCREAMv1-DP.eamxx-dpxx-dycomsrf01",
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FIOP-SCREAMv1-DP.eamxx-dpxx-arm97",
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FIOP-SCREAMv1-DP.eamxx-dpxx-comble",
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FRCE-SCREAMv1-DP",
            )
    },
    "e3sm_eamxx_v1_medres" : {
        "time"  : "02:00:00",
        "tests" : (
            "ERS_Ln22.ne30_ne30.F2010-SCREAMv1.eamxx-internal_diagnostics_level--eamxx-output-preset-3",
            "PEM_Ln90.ne30pg2_ne30pg2.F2010-SCREAMv1.eamxx-spa_remap--eamxx-output-preset-4",
            "ERS_Ln90.ne30pg2_ne30pg2.F2010-SCREAMv1.eamxx-small_kernels--eamxx-output-preset-5",
            "ERP_Ln22.conusx4v1pg2_r05_oECv3.F2010-SCREAMv1-noAero.eamxx-bfbhash--eamxx-output-preset-6",
            "ERS_Ln22.ne30pg2_ne30pg2.F2010-SCREAMv1.eamxx-L128--eamxx-output-preset-4",
            "REP_Ld5.ne30pg2_ne30pg2.F2010-SCREAMv1.eamxx-L128--eamxx-output-preset-6",
            "SMS.ne30pg2_EC30to60E2r2.WCYCLXX2010",
            "ERS_Ln90.ne30pg2_ne30pg2.F2010-SCREAMv1.eamxx-L128--eamxx-sl_nsubstep2",
            )
    },
    # Used to track performance
    "e3sm_eamxx_v1_hires" : {
        "time"  : "01:00:00",
        "tests" : (
            "SMS_Ln300.ne30pg2_ne30pg2.F2010-SCREAMv1.eamxx-perf_test--eamxx-output-preset-1"
            )
    }
    ```
<!-- markdownlint-enable MD038 -->

[^dpxx]: An EAMxx compset configuration employing doubly-periodic lateral boundary conditions.
[^grid_nomenclature]:
    See below for some additional details about ***Grid Nomenclature***
    ??? Note "Notes on Grid Nomenclature"

        <!-- markdownlint-disable-line MD046 -->
        <!-- false-positive inconsistency in fenced/indented code blocks -->
        The convention of E3SM and EAMxx is to characterize grids in terms of:
    
        - The **N**umber of spectral **E**lements (`ne<X>`) along each dimension
        (e.g., length and width) of a face of the "cubed sphere" model used to
        represent the earth.
            - That is, an `ne4_ne4` grid has 16 elements on each face of the
            cubed sphere.
            - As such, the number of elements tiling the cubed-sphere earth is
            <br>
            $N_{\text{elements}} = N_{\text{sides}}\times [\text{ne}]^2 = 6 [\text{ne}]$.
            - Note that the grid resolution may be written denoting an
            `np<X>` term
            that indicates the number of Gauss-Legendre-Lobatto nodes used to
            discretize each dimension of the spectral element.
            - For `ne<X>` grids, `np` is essentially fixed at 4 for modern E3SM,
            meaning each element is discretized on an 4$\times$4 grid.
            - Thus, the total number of GLL nodes discretizing the surface of the
            cubed sphere is
            <br>
            $N_{\text{points}} = N_{\text{elements}}$
            $\times \left([\text{np}] - 1\right)^2 + 2$

[^ne4-testing-only]: ***Note:*** The `ne4` (ultra-low resolution) grid is
intended to only be used for unit testing or debugging.

<!-- markdownlint-disable MD033 -->
<!-- html error--I'm not aware of another way to make a figure and feel like
this is helpful to have -->
<figure markdown="span">
  ![Cubed-sphere-earth](https://acme-climate.atlassian.net/wiki/download/thumbnails/34113147/NE4-SEgrid_p0000.png?version=1&modificationDate=1442089084122&cacheVersion=1&api=v2&width=564&height=564){width=250}
  <figcaption>Cubed-sphere representation of the earth depicting<br> spectral
    elements (blue boxes)<br> and Gauss-Legendre-Lobatto interpolation points
    (green dots).</figcaption>
</figure>
<!-- markdownlint-enable -->
