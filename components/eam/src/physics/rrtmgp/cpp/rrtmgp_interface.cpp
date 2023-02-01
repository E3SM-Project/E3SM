#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include "mo_load_coefficients.h"
#include "mo_rte_sw.h"
#include "mo_rte_lw.h"
#include "mo_optical_props.h"
#include "rrtmgp_const.h"
#include "mo_fluxes_byband.h"

// Prototypes
extern "C" int get_nband_sw();
extern "C" int get_nband_lw();
extern "C" int get_ngpt_sw();
extern "C" int get_ngpt_lw();
extern "C" double get_min_temperature();
extern "C" double get_max_temperature();
extern "C" void get_gpoint_bands_sw(int *gpoint_bands);
extern "C" void get_gpoint_bands_lw(int *gpoint_bands);
extern "C" void rrtmgp_finalize();
extern "C" void rrtmgp_run_sw (
        int ngas, int ncol, int nlay,
        double *gas_vmr_p, double *pmid_p      , double *tmid_p      , double *pint_p,
        double *coszrs_p , double *albedo_dir_p, double *albedo_dif_p,
        double *cld_tau_gpt_p, double *cld_ssa_gpt_p, double *cld_asm_gpt_p,
        double *aer_tau_bnd_p, double *aer_ssa_bnd_p, double *aer_asm_bnd_p,
        double *allsky_flux_up_p    , double *allsky_flux_dn_p    , double *allsky_flux_net_p    , double *allsky_flux_dn_dir_p,
        double *allsky_bnd_flux_up_p, double *allsky_bnd_flux_dn_p, double *allsky_bnd_flux_net_p, double *allsky_bnd_flux_dn_dir_p,
        double *clrsky_flux_up_p    , double *clrsky_flux_dn_p    , double *clrsky_flux_net_p    , double *clrsky_flux_dn_dir_p,
        double *clrsky_bnd_flux_up_p, double *clrsky_bnd_flux_dn_p, double *clrsky_bnd_flux_net_p, double *clrsky_bnd_flux_dn_dir_p,
        double tsi_scaling
        );
extern "C" void rrtmgp_run_lw (
        int ngas, int ncol, int nlay,
        double *gas_vmr_p           , 
        double *pmid_p              , double *tmid_p              , double *pint_p               , double *tint_p,
        double *emis_sfc_p          ,
        double *cld_tau_gpt_p       , double *aer_tau_bnd_p       ,
        double *allsky_flux_up_p    , double *allsky_flux_dn_p    , double *allsky_flux_net_p    ,
        double *allsky_bnd_flux_up_p, double *allsky_bnd_flux_dn_p, double *allsky_bnd_flux_net_p,
        double *clrsky_flux_up_p    , double *clrsky_flux_dn_p    , double *clrsky_flux_net_p    ,
        double *clrsky_bnd_flux_up_p, double *clrsky_bnd_flux_dn_p, double *clrsky_bnd_flux_net_p
        );

// Objects live here in file scope because they need to be initialized just *once*
GasOpticsRRTMGP k_dist_sw;
GasOpticsRRTMGP k_dist_lw;

// Vector of strings to hold active gas names. 
string1d active_gases;

extern "C" void rrtmgp_initialize_cxx(int ngas, char *gas_names[], char const *coefficients_file_sw, char const *coefficients_file_lw) {
    // First, make sure yakl has been initialized
    if (!yakl::isInitialized()) {
        yakl::init();
    }

    // Read gas optics coefficients from file
    // Need to initialize available_gases here! The only field of the
    // available_gases type that is used in the kdist initialize is
    // available_gases%gas_name, which gives the name of each gas that would be
    // present in the ty_gas_concs object. So, we can just set this here, rather
    // than trying to fully populate the ty_gas_concs object here, which would be
    // impossible from this initialization routine because I do not think the
    // rad_cnst objects are setup yet.
    // the other tasks!
    active_gases = string1d("active_gases", ngas);
    for (int igas=0; igas<ngas; igas++) {
        active_gases(igas+1) = gas_names[igas];
    }
    GasConcs available_gases;
    available_gases.init(active_gases, 1, 1);
    load_and_init(k_dist_sw, coefficients_file_sw, available_gases);
    load_and_init(k_dist_lw, coefficients_file_lw, available_gases);
}

extern "C" void rrtmgp_finalize() {
    k_dist_sw.finalize();
    k_dist_lw.finalize();
    yakl::finalize();
}

extern "C" int get_nband_sw() {
    return k_dist_sw.get_nband();
}

extern "C" int get_nband_lw() {
    return k_dist_lw.get_nband();
}

extern "C" int get_ngpt_sw() {
    return k_dist_sw.get_ngpt();
}

extern "C" int get_ngpt_lw() {
    return k_dist_lw.get_ngpt();
}

extern "C" double get_min_temperature() {
    return std::min(k_dist_sw.temp_ref_min, k_dist_lw.temp_ref_min);
}

extern "C" double get_max_temperature() {
    return std::max(k_dist_sw.temp_ref_max, k_dist_lw.temp_ref_max);
}

extern "C" void get_gpoint_bands_sw(int *gpoint_bands_p) {
    auto gpoint_bands_sw = intHost1d("gpoint_bands", gpoint_bands_p, k_dist_sw.get_ngpt());
    auto tmp = k_dist_sw.get_gpoint_bands();
    tmp.deep_copy_to(gpoint_bands_sw);
    yakl::fence();
}

extern "C" void get_gpoint_bands_lw(int *gpoint_bands_p) {
    auto gpoint_bands_lw = intHost1d("gpoint_bands", gpoint_bands_p, k_dist_lw.get_ngpt());
    auto tmp = k_dist_lw.get_gpoint_bands();
    tmp.deep_copy_to(gpoint_bands_lw);
    yakl::fence();
}

extern "C" void rrtmgp_run_sw (
        int ngas, int ncol, int nlay,
        double *gas_vmr_p, double *pmid_p      , double *tmid_p      , double *pint_p,
        double *coszrs_p , double *albedo_dir_p, double *albedo_dif_p,
        double *cld_tau_gpt_p, double *cld_ssa_gpt_p, double *cld_asm_gpt_p,
        double *aer_tau_bnd_p, double *aer_ssa_bnd_p, double *aer_asm_bnd_p,
        double *allsky_flux_up_p    , double *allsky_flux_dn_p    , double *allsky_flux_net_p    , double *allsky_flux_dn_dir_p,
        double *allsky_bnd_flux_up_p, double *allsky_bnd_flux_dn_p, double *allsky_bnd_flux_net_p, double *allsky_bnd_flux_dn_dir_p,
        double *clrsky_flux_up_p    , double *clrsky_flux_dn_p    , double *clrsky_flux_net_p    , double *clrsky_flux_dn_dir_p,
        double *clrsky_bnd_flux_up_p, double *clrsky_bnd_flux_dn_p, double *clrsky_bnd_flux_net_p, double *clrsky_bnd_flux_dn_dir_p,
        double tsi_scaling
        ) {

    using yakl::fortran::parallel_for;
    using yakl::fortran::Bounds;

    // Wrap pointers in YAKL arrays
    int nswbands = k_dist_sw.get_nband();
    int nswgpts  = k_dist_sw.get_ngpt();
    auto gas_vmr_host                = realHost3d("gas_vmr", gas_vmr_p, ngas, ncol, nlay);
    auto pmid_host                   = realHost2d("pmid", pmid_p, ncol, nlay);
    auto tmid_host                   = realHost2d("tmid", tmid_p, ncol, nlay);
    auto pint_host                   = realHost2d("pint", pint_p, ncol, nlay+1);
    auto coszrs_host                 = realHost1d("coszrs", coszrs_p, ncol);
    auto albedo_dir_host             = realHost2d("albedo_dir", albedo_dir_p, nswbands, ncol);
    auto albedo_dif_host             = realHost2d("albedo_dif", albedo_dif_p, nswbands, ncol);
    auto cld_tau_gpt_host            = realHost3d("cld_tau_gpt", cld_tau_gpt_p, ncol, nlay, nswgpts);
    auto cld_ssa_gpt_host            = realHost3d("cld_ssa_gpt", cld_ssa_gpt_p, ncol, nlay, nswgpts);
    auto cld_asm_gpt_host            = realHost3d("cld_asm_gpt", cld_asm_gpt_p, ncol, nlay, nswgpts);
    auto aer_tau_bnd_host            = realHost3d("aer_tau_bnd", aer_tau_bnd_p, ncol, nlay, nswbands);
    auto aer_ssa_bnd_host            = realHost3d("aer_ssa_bnd", aer_ssa_bnd_p, ncol, nlay, nswbands);
    auto aer_asm_bnd_host            = realHost3d("aer_asm_bnd", aer_asm_bnd_p, ncol, nlay, nswbands);
    auto allsky_flux_up_host         = realHost2d("allsky_flux_up", allsky_flux_up_p, ncol, nlay+1);
    auto allsky_flux_dn_host         = realHost2d("allsky_flux_dn", allsky_flux_dn_p, ncol, nlay+1);
    auto allsky_flux_dn_dir_host     = realHost2d("allsky_flux_dn_dir", allsky_flux_dn_dir_p, ncol, nlay+1);
    auto allsky_flux_net_host        = realHost2d("allsky_flux_net", allsky_flux_net_p, ncol, nlay+1);
    auto clrsky_flux_up_host         = realHost2d("clrsky_flux_up", clrsky_flux_up_p, ncol, nlay+1);
    auto clrsky_flux_dn_host         = realHost2d("clrsky_flux_dn", clrsky_flux_dn_p, ncol, nlay+1);
    auto clrsky_flux_dn_dir_host     = realHost2d("clrsky_flux_dn_dir", clrsky_flux_dn_dir_p, ncol, nlay+1);
    auto clrsky_flux_net_host        = realHost2d("clrsky_flux_net", clrsky_flux_net_p, ncol, nlay+1);
    auto allsky_bnd_flux_up_host     = realHost3d("allsky_bnd_flux_up", allsky_bnd_flux_up_p, ncol, nlay+1, nswbands);
    auto allsky_bnd_flux_dn_host     = realHost3d("allsky_bnd_flux_dn", allsky_bnd_flux_dn_p, ncol, nlay+1, nswbands);
    auto allsky_bnd_flux_dn_dir_host = realHost3d("allsky_bnd_flux_dn_dir", allsky_bnd_flux_dn_dir_p, ncol, nlay+1, nswbands);
    auto allsky_bnd_flux_net_host    = realHost3d("allsky_bnd_flux_net", allsky_bnd_flux_net_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_up_host     = realHost3d("clrsky_bnd_flux_up", clrsky_bnd_flux_up_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_dn_host     = realHost3d("clrsky_bnd_flux_dn", clrsky_bnd_flux_dn_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_dn_dir_host = realHost3d("clrsky_bnd_flux_dn_dir", clrsky_bnd_flux_dn_dir_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_net_host    = realHost3d("clrsky_bnd_flux_net", clrsky_bnd_flux_net_p, ncol, nlay+1, nswbands);

    real3d gas_vmr               ("gas_vmr", ngas, ncol, nlay);
    real2d pmid                  ("pmid", ncol, nlay);
    real2d tmid                  ("tmid", ncol, nlay);
    real2d pint                  ("pint", ncol, nlay+1);
    real1d coszrs                ("coszrs", ncol);
    real2d albedo_dir            ("albedo_dir", nswbands, ncol);
    real2d albedo_dif            ("albedo_dif", nswbands, ncol);
    real3d cld_tau_gpt           ("cld_tau_gpt", ncol, nlay, nswgpts);
    real3d cld_ssa_gpt           ("cld_ssa_gpt", ncol, nlay, nswgpts);
    real3d cld_asm_gpt           ("cld_asm_gpt", ncol, nlay, nswgpts);
    real3d aer_tau_bnd           ("aer_tau_bnd", ncol, nlay, nswbands);
    real3d aer_ssa_bnd           ("aer_ssa_bnd", ncol, nlay, nswbands);
    real3d aer_asm_bnd           ("aer_asm_bnd", ncol, nlay, nswbands);
    real2d allsky_flux_up        ("allsky_flux_up", ncol, nlay+1);
    real2d allsky_flux_dn        ("allsky_flux_dn", ncol, nlay+1);
    real2d allsky_flux_dn_dir    ("allsky_flux_dn_dir", ncol, nlay+1);
    real2d allsky_flux_net       ("allsky_flux_net", ncol, nlay+1);
    real2d clrsky_flux_up        ("clrsky_flux_up", ncol, nlay+1);
    real2d clrsky_flux_dn        ("clrsky_flux_dn", ncol, nlay+1);
    real2d clrsky_flux_dn_dir    ("clrsky_flux_dn_dir", ncol, nlay+1);
    real2d clrsky_flux_net       ("clrsky_flux_net", ncol, nlay+1);
    real3d allsky_bnd_flux_up    ("allsky_bnd_flux_up", ncol, nlay+1, nswbands);
    real3d allsky_bnd_flux_dn    ("allsky_bnd_flux_dn", ncol, nlay+1, nswbands);
    real3d allsky_bnd_flux_dn_dir("allsky_bnd_flux_dn_dir", ncol, nlay+1, nswbands);
    real3d allsky_bnd_flux_net   ("allsky_bnd_flux_net", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_up    ("clrsky_bnd_flux_up", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_dn    ("clrsky_bnd_flux_dn", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_dn_dir("clrsky_bnd_flux_dn_dir", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_net   ("clrsky_bnd_flux_net", ncol, nlay+1, nswbands);

    // TODO: Only copy in the inputs
    gas_vmr_host               .deep_copy_to(gas_vmr               );
    pmid_host                  .deep_copy_to(pmid                  );
    tmid_host                  .deep_copy_to(tmid                  );
    pint_host                  .deep_copy_to(pint                  );
    coszrs_host                .deep_copy_to(coszrs                );
    albedo_dir_host            .deep_copy_to(albedo_dir            );
    albedo_dif_host            .deep_copy_to(albedo_dif            );
    cld_tau_gpt_host           .deep_copy_to(cld_tau_gpt           );
    cld_ssa_gpt_host           .deep_copy_to(cld_ssa_gpt           );
    cld_asm_gpt_host           .deep_copy_to(cld_asm_gpt           );
    aer_tau_bnd_host           .deep_copy_to(aer_tau_bnd           );
    aer_ssa_bnd_host           .deep_copy_to(aer_ssa_bnd           );
    aer_asm_bnd_host           .deep_copy_to(aer_asm_bnd           );
    //allsky_flux_up_host        .deep_copy_to(allsky_flux_up        );
    //allsky_flux_dn_host        .deep_copy_to(allsky_flux_dn        );
    //allsky_flux_dn_dir_host    .deep_copy_to(allsky_flux_dn_dir    );
    //allsky_flux_net_host       .deep_copy_to(allsky_flux_net       );
    //clrsky_flux_up_host        .deep_copy_to(clrsky_flux_up        );
    //clrsky_flux_dn_host        .deep_copy_to(clrsky_flux_dn        );
    //clrsky_flux_dn_dir_host    .deep_copy_to(clrsky_flux_dn_dir    );
    //clrsky_flux_net_host       .deep_copy_to(clrsky_flux_net       );
    //allsky_bnd_flux_up_host    .deep_copy_to(allsky_bnd_flux_up    );
    //allsky_bnd_flux_dn_host    .deep_copy_to(allsky_bnd_flux_dn    );
    //allsky_bnd_flux_dn_dir_host.deep_copy_to(allsky_bnd_flux_dn_dir);
    //allsky_bnd_flux_net_host   .deep_copy_to(allsky_bnd_flux_net   );
    //clrsky_bnd_flux_up_host    .deep_copy_to(clrsky_bnd_flux_up    );
    //clrsky_bnd_flux_dn_host    .deep_copy_to(clrsky_bnd_flux_dn    );
    //clrsky_bnd_flux_dn_dir_host.deep_copy_to(clrsky_bnd_flux_dn_dir);
    //clrsky_bnd_flux_net_host   .deep_copy_to(clrsky_bnd_flux_net   );


    // Populate gas concentrations object
    GasConcs gas_concs;
    gas_concs.init(active_gases, ncol, nlay);
    real2d tmp2d;
    tmp2d = real2d("tmp", ncol, nlay);
    for (int igas = 1; igas <= ngas; igas++) {
        parallel_for(Bounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
            tmp2d(icol,ilay) = gas_vmr(igas,icol,ilay);
        });
        gas_concs.set_vmr(active_gases(igas), tmp2d);
    }

    // Do gas optics
    // TODO: should we avoid allocating here?
    OpticalProps2str combined_optics;
    combined_optics.alloc_2str(ncol, nlay, k_dist_sw);
    bool top_at_1 = pmid_host(1, 1) < pmid_host (1, 2);
    real2d toa_flux("toa_flux", ncol, nswgpts);
    k_dist_sw.gas_optics(ncol, nlay, top_at_1, pmid, pint, tmid, gas_concs, combined_optics, toa_flux);

    // Apply TOA flux scaling
    parallel_for(Bounds<2>(nswgpts,ncol), YAKL_LAMBDA (int igpt, int icol) {
        toa_flux(icol, igpt) = tsi_scaling * toa_flux(icol, igpt);
    });

    // Add in aerosol
    // TODO: should we avoid allocating here?
    OpticalProps2str aerosol_optics;
    auto &aerosol_optics_tau = aerosol_optics.tau;
    auto &aerosol_optics_ssa = aerosol_optics.ssa;
    auto &aerosol_optics_g   = aerosol_optics.g  ;
    if (true) {
        aerosol_optics.alloc_2str(ncol, nlay, k_dist_sw);
        auto gpt_bnd = aerosol_optics.get_gpoint_bands();
        parallel_for(Bounds<3>(nswgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
            aerosol_optics_tau(icol,ilay,igpt) = aer_tau_bnd(icol,ilay,gpt_bnd(igpt));
            aerosol_optics_ssa(icol,ilay,igpt) = aer_ssa_bnd(icol,ilay,gpt_bnd(igpt));
            aerosol_optics_g  (icol,ilay,igpt) = aer_asm_bnd(icol,ilay,gpt_bnd(igpt));
        });
    } else {
        aerosol_optics.alloc_2str(ncol, nlay, k_dist_sw.get_band_lims_wavenumber());
        parallel_for(Bounds<3>(nswbands,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol) {
            aerosol_optics_tau(icol,ilay,ibnd) = aer_tau_bnd(icol,ilay,ibnd);
            aerosol_optics_ssa(icol,ilay,ibnd) = aer_ssa_bnd(icol,ilay,ibnd);
            aerosol_optics_g  (icol,ilay,ibnd) = aer_asm_bnd(icol,ilay,ibnd);
        });
    }
    aerosol_optics.delta_scale();
    aerosol_optics.increment(combined_optics);

    // Do the clearsky calculation before adding in clouds
    FluxesByband fluxes_clrsky;
    fluxes_clrsky.flux_up = real2d("clrsky_flux_up", ncol, nlay+1); // clrsky_flux_up;
    fluxes_clrsky.flux_dn = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_dn;
    fluxes_clrsky.flux_dn_dir = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_dn_dir;
    fluxes_clrsky.flux_net = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_net;
    fluxes_clrsky.bnd_flux_up = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_up;
    fluxes_clrsky.bnd_flux_dn = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_dn;
    fluxes_clrsky.bnd_flux_dn_dir = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_dn_dir;
    fluxes_clrsky.bnd_flux_net = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_net;
    rte_sw(combined_optics, top_at_1, coszrs, toa_flux, albedo_dir, albedo_dif, fluxes_clrsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_clrsky.flux_up.deep_copy_to(clrsky_flux_up);
    fluxes_clrsky.flux_dn.deep_copy_to(clrsky_flux_dn);
    fluxes_clrsky.flux_dn_dir.deep_copy_to(clrsky_flux_dn_dir);
    fluxes_clrsky.flux_net.deep_copy_to(clrsky_flux_net);
    fluxes_clrsky.bnd_flux_up.deep_copy_to(clrsky_bnd_flux_up);
    fluxes_clrsky.bnd_flux_dn.deep_copy_to(clrsky_bnd_flux_dn);
    fluxes_clrsky.bnd_flux_dn_dir.deep_copy_to(clrsky_bnd_flux_dn_dir);
    fluxes_clrsky.bnd_flux_net.deep_copy_to(clrsky_bnd_flux_net);

    // Add in clouds
    OpticalProps2str cloud_optics;
    cloud_optics.alloc_2str(ncol, nlay, k_dist_sw);
    auto &cloud_optics_tau = cloud_optics.tau;
    auto &cloud_optics_ssa = cloud_optics.ssa;
    auto &cloud_optics_g   = cloud_optics.g  ;
    parallel_for(Bounds<3>(nswgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
        cloud_optics_tau(icol,ilay,igpt) = cld_tau_gpt(icol,ilay,igpt);
        cloud_optics_ssa(icol,ilay,igpt) = cld_ssa_gpt(icol,ilay,igpt);
        cloud_optics_g  (icol,ilay,igpt) = cld_asm_gpt(icol,ilay,igpt);
    });
    cloud_optics.delta_scale();
    cloud_optics.increment(combined_optics);

    // Call SW flux driver
    FluxesByband fluxes_allsky;
    fluxes_allsky.flux_up = real2d("allsky_flux_up", ncol, nlay+1); // allsky_flux_up;
    fluxes_allsky.flux_dn = real2d("allsky_flux_up", ncol, nlay+1); //allsky_flux_dn;
    fluxes_allsky.flux_dn_dir = real2d("allsky_flux_up", ncol, nlay+1); //allsky_flux_dn_dir;
    fluxes_allsky.flux_net = real2d("allsky_flux_up", ncol, nlay+1); //allsky_flux_net;
    fluxes_allsky.bnd_flux_up = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_up;
    fluxes_allsky.bnd_flux_dn = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_dn;
    fluxes_allsky.bnd_flux_dn_dir = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_dn_dir;
    fluxes_allsky.bnd_flux_net = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_net;
    rte_sw(combined_optics, top_at_1, coszrs, toa_flux, albedo_dir, albedo_dif, fluxes_allsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_allsky.flux_up.deep_copy_to(allsky_flux_up);
    fluxes_allsky.flux_dn.deep_copy_to(allsky_flux_dn);
    fluxes_allsky.flux_dn_dir.deep_copy_to(allsky_flux_dn_dir);
    fluxes_allsky.flux_net.deep_copy_to(allsky_flux_net);
    fluxes_allsky.bnd_flux_up.deep_copy_to(allsky_bnd_flux_up);
    fluxes_allsky.bnd_flux_dn.deep_copy_to(allsky_bnd_flux_dn);
    fluxes_allsky.bnd_flux_dn_dir.deep_copy_to(allsky_bnd_flux_dn_dir);
    fluxes_allsky.bnd_flux_net.deep_copy_to(allsky_bnd_flux_net);

    // TODO: Only copy out the outputs
    //gas_vmr               .deep_copy_to(gas_vmr_host               );
    //pmid                  .deep_copy_to(pmid_host                  );
    //tmid                  .deep_copy_to(tmid_host                  );
    //pint                  .deep_copy_to(pint_host                  );
    //coszrs                .deep_copy_to(coszrs_host                );
    //albedo_dir            .deep_copy_to(albedo_dir_host            );
    //albedo_dif            .deep_copy_to(albedo_dif_host            );
    //cld_tau_gpt           .deep_copy_to(cld_tau_gpt_host           );
    //cld_ssa_gpt           .deep_copy_to(cld_ssa_gpt_host           );
    //cld_asm_gpt           .deep_copy_to(cld_asm_gpt_host           );
    //aer_tau_bnd           .deep_copy_to(aer_tau_bnd_host           );
    //aer_ssa_bnd           .deep_copy_to(aer_ssa_bnd_host           );
    //aer_asm_bnd           .deep_copy_to(aer_asm_bnd_host           );
    allsky_flux_up        .deep_copy_to(allsky_flux_up_host        );
    allsky_flux_dn        .deep_copy_to(allsky_flux_dn_host        );
    allsky_flux_dn_dir    .deep_copy_to(allsky_flux_dn_dir_host    );
    allsky_flux_net       .deep_copy_to(allsky_flux_net_host       );
    clrsky_flux_up        .deep_copy_to(clrsky_flux_up_host        );
    clrsky_flux_dn        .deep_copy_to(clrsky_flux_dn_host        );
    clrsky_flux_dn_dir    .deep_copy_to(clrsky_flux_dn_dir_host    );
    clrsky_flux_net       .deep_copy_to(clrsky_flux_net_host       );
    allsky_bnd_flux_up    .deep_copy_to(allsky_bnd_flux_up_host    );
    allsky_bnd_flux_dn    .deep_copy_to(allsky_bnd_flux_dn_host    );
    allsky_bnd_flux_dn_dir.deep_copy_to(allsky_bnd_flux_dn_dir_host);
    allsky_bnd_flux_net   .deep_copy_to(allsky_bnd_flux_net_host   );
    clrsky_bnd_flux_up    .deep_copy_to(clrsky_bnd_flux_up_host    );
    clrsky_bnd_flux_dn    .deep_copy_to(clrsky_bnd_flux_dn_host    );
    clrsky_bnd_flux_dn_dir.deep_copy_to(clrsky_bnd_flux_dn_dir_host);
    clrsky_bnd_flux_net   .deep_copy_to(clrsky_bnd_flux_net_host   );
    yakl::fence();

} 

extern "C" void rrtmgp_run_lw (
        int ngas, int ncol, int nlay,
        double *gas_vmr_p           , 
        double *pmid_p              , double *tmid_p              , double *pint_p               , double *tint_p,
        double *emis_sfc_p          ,
        double *cld_tau_gpt_p       , double *aer_tau_bnd_p       ,
        double *allsky_flux_up_p    , double *allsky_flux_dn_p    , double *allsky_flux_net_p    ,
        double *allsky_bnd_flux_up_p, double *allsky_bnd_flux_dn_p, double *allsky_bnd_flux_net_p,
        double *clrsky_flux_up_p    , double *clrsky_flux_dn_p    , double *clrsky_flux_net_p    ,
        double *clrsky_bnd_flux_up_p, double *clrsky_bnd_flux_dn_p, double *clrsky_bnd_flux_net_p
        ) {

    using yakl::fortran::parallel_for;
    using yakl::fortran::Bounds;

    // Wrap pointers in YAKL arrays
    int nlwbands = k_dist_lw.get_nband();
    int nlwgpts  = k_dist_lw.get_ngpt();
    auto gas_vmr_host             = realHost3d("gas_vmr", gas_vmr_p, ngas, ncol, nlay);
    auto pmid_host                = realHost2d("pmid", pmid_p, ncol, nlay);
    auto tmid_host                = realHost2d("tmid", tmid_p, ncol, nlay);
    auto pint_host                = realHost2d("pint", pint_p, ncol, nlay+1);
    auto tint_host                = realHost2d("tint", tint_p, ncol, nlay+1);
    auto emis_sfc_host            = realHost2d("emis_sfc", emis_sfc_p, nlwbands, ncol);
    auto cld_tau_gpt_host         = realHost3d("cld_tau_gpt", cld_tau_gpt_p, ncol, nlay, nlwgpts);
    auto aer_tau_bnd_host         = realHost3d("aer_tau_bnd", aer_tau_bnd_p, ncol, nlay, nlwbands);
    auto allsky_flux_up_host      = realHost2d("allsky_flux_up", allsky_flux_up_p, ncol, nlay+1);
    auto allsky_flux_dn_host      = realHost2d("allsky_flux_dn", allsky_flux_dn_p, ncol, nlay+1);
    auto allsky_flux_net_host     = realHost2d("allsky_flux_net", allsky_flux_net_p, ncol, nlay+1);
    auto clrsky_flux_up_host      = realHost2d("clrsky_flux_up", clrsky_flux_up_p, ncol, nlay+1);
    auto clrsky_flux_dn_host      = realHost2d("clrsky_flux_dn", clrsky_flux_dn_p, ncol, nlay+1);
    auto clrsky_flux_net_host     = realHost2d("clrsky_flux_net", clrsky_flux_net_p, ncol, nlay+1);
    auto allsky_bnd_flux_up_host  = realHost3d("allsky_bnd_flux_up", allsky_bnd_flux_up_p, ncol, nlay+1, nlwbands);
    auto allsky_bnd_flux_dn_host  = realHost3d("allsky_bnd_flux_dn", allsky_bnd_flux_dn_p, ncol, nlay+1, nlwbands);
    auto allsky_bnd_flux_net_host = realHost3d("allsky_bnd_flux_net", allsky_bnd_flux_net_p, ncol, nlay+1, nlwbands);
    auto clrsky_bnd_flux_up_host  = realHost3d("clrsky_bnd_flux_up", clrsky_bnd_flux_up_p, ncol, nlay+1, nlwbands);
    auto clrsky_bnd_flux_dn_host  = realHost3d("clrsky_bnd_flux_dn", clrsky_bnd_flux_dn_p, ncol, nlay+1, nlwbands);
    auto clrsky_bnd_flux_net_host = realHost3d("clrsky_bnd_flux_net", clrsky_bnd_flux_net_p, ncol, nlay+1, nlwbands);

    real3d gas_vmr            ("gas_vmr", ngas, ncol, nlay);
    real2d pmid               ("pmid", ncol, nlay);
    real2d tmid               ("tmid", ncol, nlay);
    real2d pint               ("pint", ncol, nlay+1);
    real2d tint               ("tint", ncol, nlay+1);
    real2d emis_sfc           ("emis_sfc", nlwbands, ncol);
    real3d cld_tau_gpt        ("cld_tau_gpt", ncol, nlay, nlwgpts);
    real3d aer_tau_bnd        ("aer_tau_bnd", ncol, nlay, nlwbands);
    real2d allsky_flux_up     ("allsky_flux_up", ncol, nlay+1);
    real2d allsky_flux_dn     ("allsky_flux_dn", ncol, nlay+1);
    real2d allsky_flux_net    ("allsky_flux_net", ncol, nlay+1);
    real2d clrsky_flux_up     ("clrsky_flux_up", ncol, nlay+1);
    real2d clrsky_flux_dn     ("clrsky_flux_dn", ncol, nlay+1);
    real2d clrsky_flux_net    ("clrsky_flux_net", ncol, nlay+1);
    real3d allsky_bnd_flux_up ("allsky_bnd_flux_up", ncol, nlay+1, nlwbands);
    real3d allsky_bnd_flux_dn ("allsky_bnd_flux_dn", ncol, nlay+1, nlwbands);
    real3d allsky_bnd_flux_net("allsky_bnd_flux_net", ncol, nlay+1, nlwbands);
    real3d clrsky_bnd_flux_up ("clrsky_bnd_flux_up", ncol, nlay+1, nlwbands);
    real3d clrsky_bnd_flux_dn ("clrsky_bnd_flux_dn", ncol, nlay+1, nlwbands);
    real3d clrsky_bnd_flux_net("clrsky_bnd_flux_net", ncol, nlay+1, nlwbands);

    // TODO: Only copy in the inputs
    gas_vmr_host            .deep_copy_to(gas_vmr            );
    pmid_host               .deep_copy_to(pmid               );
    tmid_host               .deep_copy_to(tmid               );
    pint_host               .deep_copy_to(pint               );
    tint_host               .deep_copy_to(tint               );
    emis_sfc_host           .deep_copy_to(emis_sfc           );
    cld_tau_gpt_host        .deep_copy_to(cld_tau_gpt        );
    aer_tau_bnd_host        .deep_copy_to(aer_tau_bnd        );
    allsky_flux_up_host     .deep_copy_to(allsky_flux_up     );
    allsky_flux_dn_host     .deep_copy_to(allsky_flux_dn     );
    allsky_flux_net_host    .deep_copy_to(allsky_flux_net    );
    clrsky_flux_up_host     .deep_copy_to(clrsky_flux_up     );
    clrsky_flux_dn_host     .deep_copy_to(clrsky_flux_dn     );
    clrsky_flux_net_host    .deep_copy_to(clrsky_flux_net    );
    allsky_bnd_flux_up_host .deep_copy_to(allsky_bnd_flux_up );
    allsky_bnd_flux_dn_host .deep_copy_to(allsky_bnd_flux_dn );
    allsky_bnd_flux_net_host.deep_copy_to(allsky_bnd_flux_net);
    clrsky_bnd_flux_up_host .deep_copy_to(clrsky_bnd_flux_up );
    clrsky_bnd_flux_dn_host .deep_copy_to(clrsky_bnd_flux_dn );
    clrsky_bnd_flux_net_host.deep_copy_to(clrsky_bnd_flux_net);

    // Populate gas concentrations
    GasConcs gas_concs;
    gas_concs.init(active_gases, ncol, nlay);
    real2d tmp2d;
    tmp2d = real2d("tmp", ncol, nlay);
    for (int igas = 1; igas <= ngas; igas++) {
        parallel_for(Bounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
            tmp2d(icol,ilay) = gas_vmr(igas,icol,ilay);
        });
        gas_concs.set_vmr(active_gases(igas), tmp2d);
    }

    //  Boundary conditions
    SourceFuncLW lw_sources;
    lw_sources.alloc(ncol, nlay, k_dist_lw);

    // Weights and angle secants for first order (k=1) Gaussian quadrature.
    //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
    //   after Abramowitz & Stegun 1972, page 921
    int constexpr max_gauss_pts = 4;
    realHost2d gauss_Ds_host ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    gauss_Ds_host(1,1) = 1.66      ; gauss_Ds_host(2,1) =         0.; gauss_Ds_host(3,1) =         0.; gauss_Ds_host(4,1) =         0.;
    gauss_Ds_host(1,2) = 1.18350343; gauss_Ds_host(2,2) = 2.81649655; gauss_Ds_host(3,2) =         0.; gauss_Ds_host(4,2) =         0.;
    gauss_Ds_host(1,3) = 1.09719858; gauss_Ds_host(2,3) = 1.69338507; gauss_Ds_host(3,3) = 4.70941630; gauss_Ds_host(4,3) =         0.;
    gauss_Ds_host(1,4) = 1.06056257; gauss_Ds_host(2,4) = 1.38282560; gauss_Ds_host(3,4) = 2.40148179; gauss_Ds_host(4,4) = 7.15513024;

    realHost2d gauss_wts_host("gauss_wts",max_gauss_pts,max_gauss_pts);
    gauss_wts_host(1,1) = 0.5         ; gauss_wts_host(2,1) = 0.          ; gauss_wts_host(3,1) = 0.          ; gauss_wts_host(4,1) = 0.          ;
    gauss_wts_host(1,2) = 0.3180413817; gauss_wts_host(2,2) = 0.1819586183; gauss_wts_host(3,2) = 0.          ; gauss_wts_host(4,2) = 0.          ;
    gauss_wts_host(1,3) = 0.2009319137; gauss_wts_host(2,3) = 0.2292411064; gauss_wts_host(3,3) = 0.0698269799; gauss_wts_host(4,3) = 0.          ;
    gauss_wts_host(1,4) = 0.1355069134; gauss_wts_host(2,4) = 0.2034645680; gauss_wts_host(3,4) = 0.1298475476; gauss_wts_host(4,4) = 0.0311809710;

    real2d gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    real2d gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
    gauss_Ds_host .deep_copy_to(gauss_Ds );
    gauss_wts_host.deep_copy_to(gauss_wts);

    // Populate optical property objects
    OpticalProps1scl combined_optics;
    combined_optics.alloc_1scl(ncol, nlay, k_dist_lw);
    bool top_at_1 = pmid_host(1, 1) < pmid_host (1, 2);
    real1d t_sfc("t_sfc", ncol);
    parallel_for(Bounds<1>(ncol), YAKL_LAMBDA (int icol) {
        t_sfc(icol) = tint(icol,nlay+1);
    });
    k_dist_lw.gas_optics(ncol, nlay, top_at_1, pmid, pint, tmid, t_sfc, gas_concs, combined_optics, lw_sources, real2d(), tint);

    // Add in aerosol; we can define this by bands or gpoints. If we define by
    // bands, then internally when increment() is called it will map these to
    // gpoints. Not sure if there is a beneift one way or another.
    OpticalProps1scl aerosol_optics;
    auto &aerosol_optics_tau = aerosol_optics.tau;
    if (false) {
        aerosol_optics.alloc_1scl(ncol, nlay, k_dist_lw);
        auto gpt_bnd = aerosol_optics.get_gpoint_bands();
        parallel_for(Bounds<3>(nlwgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
            aerosol_optics_tau(icol,ilay,igpt) = aer_tau_bnd(icol,ilay,gpt_bnd(igpt));
        });
    } else {
        aerosol_optics.alloc_1scl(ncol, nlay, k_dist_lw.get_band_lims_wavenumber());
        parallel_for(Bounds<3>(nlwbands,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol) {
            aerosol_optics_tau(icol,ilay,ibnd) = aer_tau_bnd(icol,ilay,ibnd);
        });
    }
    aerosol_optics.increment(combined_optics);

    // Do the clearsky calculation before adding in clouds
    FluxesByband fluxes_clrsky;
    fluxes_clrsky.flux_up = real2d("clrsky_flux_up", ncol, nlay+1); // clrsky_flux_up;
    fluxes_clrsky.flux_dn = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_dn;
    fluxes_clrsky.flux_net = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_net;
    fluxes_clrsky.bnd_flux_up = real3d("clrsky_flux_up", ncol, nlay+1, nlwbands); //clrsky_bnd_flux_up;
    fluxes_clrsky.bnd_flux_dn = real3d("clrsky_flux_up", ncol, nlay+1, nlwbands); //clrsky_bnd_flux_dn;
    fluxes_clrsky.bnd_flux_net = real3d("clrsky_flux_up", ncol, nlay+1, nlwbands); //clrsky_bnd_flux_net;
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, combined_optics, top_at_1, lw_sources, emis_sfc, fluxes_clrsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_clrsky.flux_up.deep_copy_to(clrsky_flux_up);
    fluxes_clrsky.flux_dn.deep_copy_to(clrsky_flux_dn);
    fluxes_clrsky.flux_net.deep_copy_to(clrsky_flux_net);
    fluxes_clrsky.bnd_flux_up.deep_copy_to(clrsky_bnd_flux_up);
    fluxes_clrsky.bnd_flux_dn.deep_copy_to(clrsky_bnd_flux_dn);
    fluxes_clrsky.bnd_flux_net.deep_copy_to(clrsky_bnd_flux_net);

    // Add in clouds
    OpticalProps1scl cloud_optics;
    cloud_optics.alloc_1scl(ncol, nlay, k_dist_lw);
    auto &cloud_optics_tau = cloud_optics.tau;
    parallel_for(Bounds<3>(nlwgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
        cloud_optics_tau(icol,ilay,igpt) = cld_tau_gpt(icol,ilay,igpt);
    });
    cloud_optics.increment(combined_optics);

    // Call LW flux driver
    FluxesByband fluxes_allsky;
    fluxes_allsky.flux_up = real2d("flux_up", ncol, nlay+1); //allsky_flux_up;
    fluxes_allsky.flux_dn = real2d("flux_dn", ncol, nlay+1); //allsky_flux_dn;
    fluxes_allsky.flux_net = real2d("flux_net", ncol, nlay+1); //allsky_flux_net;
    fluxes_allsky.bnd_flux_up = real3d("flux_up", ncol, nlay+1, nlwbands); //allsky_bnd_flux_up;
    fluxes_allsky.bnd_flux_dn = real3d("flux_dn", ncol, nlay+1, nlwbands); //allsky_bnd_flux_dn;
    fluxes_allsky.bnd_flux_net = real3d("flux_net", ncol, nlay+1, nlwbands); //allsky_bnd_flux_net;
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, combined_optics, top_at_1, lw_sources, emis_sfc, fluxes_allsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_allsky.flux_up.deep_copy_to(allsky_flux_up);
    fluxes_allsky.flux_dn.deep_copy_to(allsky_flux_dn);
    fluxes_allsky.flux_net.deep_copy_to(allsky_flux_net);
    fluxes_allsky.bnd_flux_up.deep_copy_to(allsky_bnd_flux_up);
    fluxes_allsky.bnd_flux_dn.deep_copy_to(allsky_bnd_flux_dn);
    fluxes_allsky.bnd_flux_net.deep_copy_to(allsky_bnd_flux_net);

    // TODO: Only copy out the outputs
    gas_vmr            .deep_copy_to(gas_vmr_host            );
    pmid               .deep_copy_to(pmid_host               );
    tmid               .deep_copy_to(tmid_host               );
    pint               .deep_copy_to(pint_host               );
    tint               .deep_copy_to(tint_host               );
    emis_sfc           .deep_copy_to(emis_sfc_host           );
    cld_tau_gpt        .deep_copy_to(cld_tau_gpt_host        );
    aer_tau_bnd        .deep_copy_to(aer_tau_bnd_host        );
    allsky_flux_up     .deep_copy_to(allsky_flux_up_host     );
    allsky_flux_dn     .deep_copy_to(allsky_flux_dn_host     );
    allsky_flux_net    .deep_copy_to(allsky_flux_net_host    );
    clrsky_flux_up     .deep_copy_to(clrsky_flux_up_host     );
    clrsky_flux_dn     .deep_copy_to(clrsky_flux_dn_host     );
    clrsky_flux_net    .deep_copy_to(clrsky_flux_net_host    );
    allsky_bnd_flux_up .deep_copy_to(allsky_bnd_flux_up_host );
    allsky_bnd_flux_dn .deep_copy_to(allsky_bnd_flux_dn_host );
    allsky_bnd_flux_net.deep_copy_to(allsky_bnd_flux_net_host);
    clrsky_bnd_flux_up .deep_copy_to(clrsky_bnd_flux_up_host );
    clrsky_bnd_flux_dn .deep_copy_to(clrsky_bnd_flux_dn_host );
    clrsky_bnd_flux_net.deep_copy_to(clrsky_bnd_flux_net_host);
    yakl::fence();

}
