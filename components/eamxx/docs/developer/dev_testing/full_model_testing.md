# Full-model Testing

Full model system testing of EAMxx is done via CIME test cases
(much like the rest of E3SM).

We offer a number of test suites, including:

- `e3sm_scream_v0`: Test the full set of V0 (pre-C++) tests
- `e3sm_scream_v1`: Test the full set of V1 (C++) tests
- `e3sm_scream_v1_at`: A smaller and quicker set of tests for autotesting (AT)
- `e3sm_scream_hires`: A small number of bigger, longer-running
(high-resolution) tests to measure performance

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

- `ne4_ne4`: low resolution
- `ne4pg2_ne4pg2`: low resolution with phys grid
- `ne30_ne30`: medium resolution
- `ne30pg2_ne30pg2`: medium resolution with phys grid
- `ne1024pg2_ne1024pg2`: ultra-high with phys grid

More information about grids may be found on this [Confluence Page](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/933986549/ATM+Grid+Resolution+Summary).

### Common EAMxx Compsets

- `F2010-SCREAM-LR`: V0 low res compset with EAMxx V0 atmosphere
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

??? Info "Lists of Test Cases and Suites from `cime_config/tests.py`"
    ```{.python .copy}
    "e3sm_atm_developer" : {
        "inherit" : ("eam_theta_pg2"),
        "tests"   : (
            "ERP_Ld3.ne4pg2_oQU480.F2010",
            "SMS_Ln9.ne4pg2_oQU480.F2010.eam-outfrq9s",
            "SMS.ne4pg2_oQU480.F2010.eam-cosplite",
            "SMS_R_Ld5.ne4_ne4.FSCM-ARM97.eam-scm",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010",
            "SMS_Ln5.ne4pg2_oQU480.F2010",
            "ERS_D.ne4pg2_oQU480.F2010.eam-hommexx",
            "SMS_Ln9_P24x1.ne4_ne4.FDPSCREAM-ARM97",
        )
    },
    "e3sm_atm_stealth" : {
        "tests"   : (
            "ERP_Ln18.ne4_oQU240.F2010.eam-cflx_cpl_2",
            "SMS_D_Ln5.ne4_oQU240.F2010.eam-cflx_cpl_2",
            "ERS.ne4pg2_oQU480.F2010.eam-p3",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-p3",
        )
    },
    "e3sm_atm_integration" : {
        "inherit" : ("eam_preqx", "eam_theta"),
        "tests" : (
            "ERP_Ln9.ne4pg2_ne4pg2.FAQP",
            "SMS_Ld1.ne4pg2_ne4pg2.FAQP.eam-clubb_only",
            "ERP_Ln9.ne4pg2_ne4pg2.FRCE",
            "PET_Ln5.ne4pg2_oQU480.F2010.allactive-mach-pet",
            "PEM_Ln5.ne4pg2_oQU480.F2010",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-cosplite_nhtfrq5",
            "SMS_Ln1.ne4pg2_oQU480.F2010.eam-chem_pp",
            "SMS_Ln5.ne30pg2_r05_IcoswISC30E3r5.BGCEXP_LNDATM_CNPRDCTC_20TR",
            "SMS_Ln5.ne30pg2_r05_IcoswISC30E3r5.BGCEXP_LNDATM_CNPRDCTC_1850",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-clubb_sp",
            "ERS_Ld5.ne4pg2_oQU480.F2010.eam-rrtmgp",
            "ERS_Ld5.ne4pg2_oQU480.F2010.eam-rrtmgpxx",
            "REP_Ln5.ne4pg2_oQU480.F2010",
            "SMS_Ld3.ne4pg2_oQU480.F2010.eam-thetahy_sl_pg2_mass",
            "ERP_Ld3.ne4pg2_ne4pg2.FIDEAL.allactive-pioroot1",
            "ERS_Ld5.ne4pg2_oQU480.F2010.eam-sathist_F2010",
        )
    },
    #atmopheric tests for extra coverage
    "e3sm_atm_extra_coverage" : {
        "tests" : (
            "SMS_Lm1.ne4pg2_oQU480.F2010",
            "ERS_Ld31.ne4pg2_oQU480.F2010",
            "ERP_Lm3.ne4pg2_oQU480.F2010",
            "SMS_D_Ln5.ne30pg2_r05_IcoswISC30E3r5.F2010",
            "ERP_Ld3.ne30pg2_r05_IcoswISC30E3r5.F2010.allactive-pioroot1",
            "SMS_Ly1.ne4pg2_oQU480.F2010",
        "SMS_D_Ln5.ne45pg2_ne45pg2.FAQP",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-implicit_stress",
            "ERS_Ld5.ne30pg2_r05_IcoswISC30E3r5.F2010.eam-implicit_stress",
            "ERP_Ld3.ne4pg2_oQU480.F2010.eam-condidiag_dcape",
            "ERP_Ld3.ne4pg2_oQU480.F2010.eam-condidiag_rhi",
        )
    },
    #atmopheric tests for hi-res
    "e3sm_atm_hi_res" : {
        "time" : "01:30:00",
        "tests" : "SMS.ne120pg2_r025_RRSwISC6to18E3r5.F2010"
        },
    #atmopheric tests to mimic low res production runs
    "e3sm_atm_prod" : {
        "tests" : (
            "SMS_Ln5.ne30pg2_r05_IcoswISC30E3r5.F2010.eam-wcprod_F2010",
            "SMS_Ld1.ne30pg2_r05_IcoswISC30E3r5.F20TR.eam-wcprod_F20TR",
        )
    },
    #atmopheric nbfb tests
    "e3sm_atm_nbfb" : {
        "tests" : (
            "PGN_P1x1.ne4pg2_oQU480.F2010",
            "TSC_PS.ne4pg2_oQU480.F2010",
            "MVK_PS.ne4pg2_oQU480.F2010",
        )
    }
    <!-- false-positive spaces in code fence marker -->
    <!-- markdownlint-disable-next-line MD038 -->
    ```

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

[^dpxx]: An EAMxx compset configuration employing doubly-periodic lateral boundary conditions.

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
