#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include "mo_load_coefficients.h"
#include "mo_rte_sw.h"
#include "mo_rte_lw.h"
#include "mo_optical_props.h"
#include "const.h"

using yakl::intrinsics::minval;
using yakl::intrinsics::maxval;

// Prototypes
extern "C" void add_gas_name(char const *gas_name);
extern "C" void convert_gas_names(string1d &gas_names);
extern "C" void vect_to_string1d(std::vector<std::string> vect, string1d strarr);
extern "C" int get_nband_sw();
extern "C" int get_nband_lw();
extern "C" int get_ngpt_sw();
extern "C" int get_ngpt_lw();
extern "C" double get_min_temperature();
extern "C" double get_max_temperature();
extern "C" void get_gpoint_bands_sw(int *gpoint_bands);
extern "C" void get_gpoint_bands_lw(int *gpoint_bands);
extern "C" void rrtmgpxx_finalize();
extern "C" void rrtmgpxx_run_sw_cpp (
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
extern "C" void rrtmgpxx_run_lw (
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

GasOpticsRRTMGP k_dist_sw;
GasOpticsRRTMGP k_dist_lw;

// Vector of strings to hold active gas names. 
// These need to be added at runtime, one by one,
// via the add_gas_name function.
std::vector<std::string> gas_names_vect;

extern "C" void rrtmgpxx_initialize_cpp(char const *coefficients_file_sw, char const *coefficients_file_lw) {
    // First, make sure yakl has been initialized
    if (!yakl::isInitialized()) {
        yakl::init();
    }

    // Read gas optics coefficients from file
    // Need to initialize available_gases here! The only field of the
    // available_gases type that is used int he kdist initialize is
    // available_gases%gas_name, which gives the name of each gas that would be
    // present in the ty_gas_concs object. So, we can just set this here, rather
    // than trying to fully populate the ty_gas_concs object here, which would be
    // impossible from this initialization routine because I do not thing the
    // rad_cnst objects are setup yet.
    // the other tasks!
    // TODO: This needs to be fixed to ONLY read in the data if masterproc, and then broadcast
    //gas_names = string1d("gas_names", 8);
    //
    // Let us cheat for a moment and hard-code the gases.
    // TODO: fix this!
    string1d gas_names("gas_names", gas_names_vect.size());
    convert_gas_names(gas_names);
    GasConcs available_gases;
    available_gases.init(gas_names, 1, 1);
    load_and_init(k_dist_sw, coefficients_file_sw, available_gases);
    load_and_init(k_dist_lw, coefficients_file_lw, available_gases);
}

extern "C" void rrtmgpxx_finalize() {
    k_dist_sw.finalize();
    k_dist_lw.finalize();
    yakl::finalize();
}

extern "C" void add_gas_name(char const *gas_name) {
    gas_names_vect.push_back(std::string(gas_name));
}

extern "C" void convert_gas_names(string1d &gas_names) {
    int ngas = gas_names_vect.size();
    if (ngas == 0) {
        throw "No active gases; are you sure you initialized gas_names_vect?";
    }
    for (int i = 1; i <= ngas; i++) {
        gas_names(i) = gas_names_vect[i-1];
    }
}

extern "C" void vect_to_string1d(std::vector<std::string> vect, string1d strarr) {
    int n = vect.size();
    for (int i = 0; i < n; i++) {
        strarr(i+1) = vect[i];
    }
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
    return min(k_dist_sw.temp_ref_min, k_dist_lw.temp_ref_min);
}

extern "C" double get_max_temperature() {
    return max(k_dist_sw.temp_ref_max, k_dist_lw.temp_ref_max);
}

extern "C" void get_gpoint_bands_sw(int *gpoint_bands_p) {
    auto gpoint_bands_sw = int1d("gpoint_bands", gpoint_bands_p, k_dist_sw.get_ngpt());
    auto tmp = k_dist_sw.get_gpoint_bands();
    tmp.deep_copy_to(gpoint_bands_sw);
}

extern "C" void get_gpoint_bands_lw(int *gpoint_bands_p) {
    auto gpoint_bands_lw = int1d("gpoint_bands", gpoint_bands_p, k_dist_lw.get_ngpt());
    auto tmp = k_dist_lw.get_gpoint_bands();
    tmp.deep_copy_to(gpoint_bands_lw);
}

extern "C" void rrtmgpxx_run_sw_cpp (
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
    // Wrap pointers in YAKL arrays
    int nswbands = k_dist_sw.get_nband();
    int nswgpts = k_dist_sw.get_ngpt();
    auto gas_vmr = real3d("gas_vmr", gas_vmr_p, ngas, ncol, nlay);
    auto pmid = real2d("pmid", pmid_p, ncol, nlay);
    auto tmid = real2d("tmid", tmid_p, ncol, nlay);
    auto pint = real2d("pint", pint_p, ncol, nlay+1);
    auto coszrs = real1d("coszrs", coszrs_p, ncol);
    auto albedo_dir = real2d("albedo_dir", albedo_dir_p, nswbands, ncol);
    auto albedo_dif = real2d("albedo_dif", albedo_dif_p, nswbands, ncol);
    auto cld_tau_gpt = real3d("cld_tau_gpt", cld_tau_gpt_p, ncol, nlay, nswgpts);
    auto cld_ssa_gpt = real3d("cld_ssa_gpt", cld_ssa_gpt_p, ncol, nlay, nswgpts);
    auto cld_asm_gpt = real3d("cld_asm_gpt", cld_asm_gpt_p, ncol, nlay, nswgpts);
    auto aer_tau_bnd = real3d("aer_tau_bnd", aer_tau_bnd_p, ncol, nlay, nswbands);
    auto aer_ssa_bnd = real3d("aer_ssa_bnd", aer_ssa_bnd_p, ncol, nlay, nswbands);
    auto aer_asm_bnd = real3d("aer_asm_bnd", aer_asm_bnd_p, ncol, nlay, nswbands);
    auto allsky_flux_up = real2d("allsky_flux_up", allsky_flux_up_p, ncol, nlay+1);
    auto allsky_flux_dn = real2d("allsky_flux_dn", allsky_flux_dn_p, ncol, nlay+1);
    auto allsky_flux_dn_dir = real2d("allsky_flux_dn_dir", allsky_flux_dn_dir_p, ncol, nlay+1);
    auto allsky_flux_net = real2d("allsky_flux_net", allsky_flux_net_p, ncol, nlay+1);
    auto clrsky_flux_up = real2d("clrsky_flux_up", clrsky_flux_up_p, ncol, nlay+1);
    auto clrsky_flux_dn = real2d("clrsky_flux_dn", clrsky_flux_dn_p, ncol, nlay+1);
    auto clrsky_flux_dn_dir = real2d("clrsky_flux_dn_dir", clrsky_flux_dn_dir_p, ncol, nlay+1);
    auto clrsky_flux_net = real2d("clrsky_flux_net", clrsky_flux_net_p, ncol, nlay+1);

    // Populate gas concentrations object
    string1d gas_names("gas_names", gas_names_vect.size());
    convert_gas_names(gas_names);
    GasConcs gas_concs;
    gas_concs.init(gas_names, ncol, nlay);
    real2d tmp2d;
    tmp2d = real2d("tmp", ncol, nlay);
    for (int igas = 1; igas <= ngas; igas++) {
        for (int icol = 1; icol <= ncol; icol++) {
            for (int ilay = 1; ilay <= nlay; ilay++) {
                tmp2d(icol,ilay) = gas_vmr(igas,icol,ilay);
            }
        }
        gas_concs.set_vmr(gas_names(igas), tmp2d);
    }

    // Do gas optics
    OpticalProps2str combined_optics;
    combined_optics.alloc_2str(ncol, nlay, k_dist_sw);
    auto pmid_host = pmid.createHostCopy();
    bool top_at_1 = pmid_host(1, 1) < pmid_host (1, 2);
    real2d toa_flux("toa_flux", ncol, nswgpts);
    k_dist_sw.gas_optics(top_at_1, pmid, pint, tmid, gas_concs, combined_optics, toa_flux);

    // Apply TOA flux scaling
    parallel_for(Bounds<2>(nswgpts,ncol), YAKL_LAMBDA (int igpt, int icol) {
        toa_flux(icol, igpt) = tsi_scaling * toa_flux(icol, igpt);
    });

    // Add in aerosol
    OpticalProps2str aerosol_optics;
    if (true) {
        aerosol_optics.alloc_2str(ncol, nlay, k_dist_sw);
        auto gpt_bnd = aerosol_optics.get_gpoint_bands();
        parallel_for(Bounds<3>(nswgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
            aerosol_optics.tau(icol,ilay,igpt) = aer_tau_bnd(icol,ilay,gpt_bnd(igpt));
            aerosol_optics.ssa(icol,ilay,igpt) = aer_ssa_bnd(icol,ilay,gpt_bnd(igpt));
            aerosol_optics.g  (icol,ilay,igpt) = aer_asm_bnd(icol,ilay,gpt_bnd(igpt));
        });
    } else {
        aerosol_optics.alloc_2str(ncol, nlay, k_dist_sw.get_band_lims_wavenumber());
        parallel_for(Bounds<3>(nswbands,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol) {
            aerosol_optics.tau(icol,ilay,ibnd) = aer_tau_bnd(icol,ilay,ibnd);
            aerosol_optics.ssa(icol,ilay,ibnd) = aer_ssa_bnd(icol,ilay,ibnd);
            aerosol_optics.g  (icol,ilay,ibnd) = aer_asm_bnd(icol,ilay,ibnd);
        });
    }
    aerosol_optics.increment(combined_optics);

    // Do the clearsky calculation before adding in clouds
    FluxesBroadband fluxes_clrsky;
    fluxes_clrsky.flux_up = clrsky_flux_up;
    fluxes_clrsky.flux_dn = clrsky_flux_dn;
    fluxes_clrsky.flux_dn_dir = clrsky_flux_dn_dir;
    fluxes_clrsky.flux_net = clrsky_flux_net;
    rte_sw(combined_optics, top_at_1, coszrs, toa_flux, albedo_dir, albedo_dif, fluxes_clrsky);

    // Add in clouds
    OpticalProps2str cloud_optics;
    cloud_optics.alloc_2str(ncol, nlay, k_dist_sw);
    parallel_for(Bounds<3>(nswgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
        cloud_optics.tau(icol,ilay,igpt) = cld_tau_gpt(icol,ilay,igpt);
        cloud_optics.ssa(icol,ilay,igpt) = cld_ssa_gpt(icol,ilay,igpt);
        cloud_optics.g  (icol,ilay,igpt) = cld_asm_gpt(icol,ilay,igpt);
    });
    cloud_optics.increment(combined_optics);

    // Call SW flux driver
    FluxesBroadband fluxes_allsky;
    fluxes_allsky.flux_up = allsky_flux_up;
    fluxes_allsky.flux_dn = allsky_flux_dn;
    fluxes_allsky.flux_dn_dir = allsky_flux_dn_dir;
    fluxes_allsky.flux_net = allsky_flux_net;
    rte_sw(combined_optics, top_at_1, coszrs, toa_flux, albedo_dir, albedo_dif, fluxes_allsky);
} 

extern "C" void rrtmgpxx_run_lw (
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
    // Wrap pointers in YAKL arrays
    int nlwbands = k_dist_lw.get_nband();
    int nlwgpts = k_dist_lw.get_ngpt();
    auto gas_vmr = real3d("gas_vmr", gas_vmr_p, ngas, ncol, nlay);
    auto pmid = real2d("pmid", pmid_p, ncol, nlay);
    auto tmid = real2d("tmid", tmid_p, ncol, nlay);
    auto pint = real2d("pint", pint_p, ncol, nlay+1);
    auto tint = real2d("tint", tint_p, ncol, nlay+1);
    auto emis_sfc = real2d("emis_sfc", emis_sfc_p, nlwbands, ncol);
    auto cld_tau_gpt = real3d("cld_tau_gpt", cld_tau_gpt_p, ncol, nlay, nlwgpts);
    auto aer_tau_bnd = real3d("aer_tau_bnd", aer_tau_bnd_p, ncol, nlay, nlwbands);
    auto allsky_flux_up = real2d("allsky_flux_up", allsky_flux_up_p, ncol, nlay+1);
    auto allsky_flux_dn = real2d("allsky_flux_dn", allsky_flux_dn_p, ncol, nlay+1);
    auto allsky_flux_net = real2d("allsky_flux_net", allsky_flux_net_p, ncol, nlay+1);
    auto clrsky_flux_up = real2d("clrsky_flux_up", clrsky_flux_up_p, ncol, nlay+1);
    auto clrsky_flux_dn = real2d("clrsky_flux_dn", clrsky_flux_dn_p, ncol, nlay+1);
    auto clrsky_flux_net = real2d("clrsky_flux_net", clrsky_flux_net_p, ncol, nlay+1);

    // Populate gas concentrations
    string1d gas_names("gas_names", gas_names_vect.size());
    convert_gas_names(gas_names);
    GasConcs gas_concs;
    gas_concs.init(gas_names, ncol, nlay);
    real2d tmp2d;
    tmp2d = real2d("tmp", ncol, nlay);
    for (int igas = 1; igas <= ngas; igas++) {
        for (int icol = 1; icol <= ncol; icol++) {
            for (int ilay = 1; ilay <= nlay; ilay++) {
                tmp2d(icol,ilay) = gas_vmr(igas,icol,ilay);
            }
        }
        gas_concs.set_vmr(gas_names(igas), tmp2d);
    }

    //  Boundary conditions
    SourceFuncLW lw_sources;
    lw_sources.alloc(ncol, nlay, k_dist_lw);

    // Weights and angle secants for first order (k=1) Gaussian quadrature.
    //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
    //   after Abramowitz & Stegun 1972, page 921
    int constexpr max_gauss_pts = 4;
    realHost2d gauss_Ds_host ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    gauss_Ds_host(1,1) = 1.66_wp      ; gauss_Ds_host(2,1) =         0._wp; gauss_Ds_host(3,1) =         0._wp; gauss_Ds_host(4,1) =         0._wp;
    gauss_Ds_host(1,2) = 1.18350343_wp; gauss_Ds_host(2,2) = 2.81649655_wp; gauss_Ds_host(3,2) =         0._wp; gauss_Ds_host(4,2) =         0._wp;
    gauss_Ds_host(1,3) = 1.09719858_wp; gauss_Ds_host(2,3) = 1.69338507_wp; gauss_Ds_host(3,3) = 4.70941630_wp; gauss_Ds_host(4,3) =         0._wp;
    gauss_Ds_host(1,4) = 1.06056257_wp; gauss_Ds_host(2,4) = 1.38282560_wp; gauss_Ds_host(3,4) = 2.40148179_wp; gauss_Ds_host(4,4) = 7.15513024_wp;

    realHost2d gauss_wts_host("gauss_wts",max_gauss_pts,max_gauss_pts);
    gauss_wts_host(1,1) = 0.5_wp         ; gauss_wts_host(2,1) = 0._wp          ; gauss_wts_host(3,1) = 0._wp          ; gauss_wts_host(4,1) = 0._wp          ;
    gauss_wts_host(1,2) = 0.3180413817_wp; gauss_wts_host(2,2) = 0.1819586183_wp; gauss_wts_host(3,2) = 0._wp          ; gauss_wts_host(4,2) = 0._wp          ;
    gauss_wts_host(1,3) = 0.2009319137_wp; gauss_wts_host(2,3) = 0.2292411064_wp; gauss_wts_host(3,3) = 0.0698269799_wp; gauss_wts_host(4,3) = 0._wp          ;
    gauss_wts_host(1,4) = 0.1355069134_wp; gauss_wts_host(2,4) = 0.2034645680_wp; gauss_wts_host(3,4) = 0.1298475476_wp; gauss_wts_host(4,4) = 0.0311809710_wp;

    real2d gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    real2d gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
    gauss_Ds_host .deep_copy_to(gauss_Ds );
    gauss_wts_host.deep_copy_to(gauss_wts);

    // Populate optical property objects
    OpticalProps1scl combined_optics;
    combined_optics.alloc_1scl(ncol, nlay, k_dist_lw);
    auto pmid_host = pmid.createHostCopy();
    bool top_at_1 = pmid_host(1, 1) < pmid_host (1, 2);
    real1d t_sfc("t_sfc", ncol);
    for (int icol=1; icol<=ncol; icol++) {
        t_sfc(icol) = tint(icol,nlay+1);
    }
    //k_dist_lw.gas_optics(top_at_1, pmid, pint, tmid, t_sfc, gas_concs, combined_optics, lw_sources, real2d(), real2d());
    k_dist_lw.gas_optics(top_at_1, pmid, pint, tmid, t_sfc, gas_concs, combined_optics, lw_sources, real2d(), tint);

    // Add in aerosol; we can define this by bands or gpoints. If we define by
    // bands, then internally when increment() is called it will map these to
    // gpoints. Not sure if there is a beneift one way or another.
    OpticalProps1scl aerosol_optics;
    if (false) {
        aerosol_optics.alloc_1scl(ncol, nlay, k_dist_lw);
        auto gpt_bnd = aerosol_optics.get_gpoint_bands();
        parallel_for(Bounds<3>(nlwgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
            aerosol_optics.tau(icol,ilay,igpt) = aer_tau_bnd(icol,ilay,gpt_bnd(igpt));
        });
    } else {
        aerosol_optics.alloc_1scl(ncol, nlay, k_dist_lw.get_band_lims_wavenumber());
        parallel_for(Bounds<3>(nlwbands,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol) {
            aerosol_optics.tau(icol,ilay,ibnd) = aer_tau_bnd(icol,ilay,ibnd);
        });
    }
    aerosol_optics.increment(combined_optics);

    // Do the clearsky calculation before adding in clouds
    FluxesBroadband fluxes_clrsky;
    fluxes_clrsky.flux_up = clrsky_flux_up;
    fluxes_clrsky.flux_dn = clrsky_flux_dn;
    fluxes_clrsky.flux_net = clrsky_flux_net;
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, combined_optics, top_at_1, lw_sources, emis_sfc, fluxes_clrsky);

    // Add in clouds
    OpticalProps1scl cloud_optics;
    cloud_optics.alloc_1scl(ncol, nlay, k_dist_lw);
    parallel_for(Bounds<3>(nlwgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
        cloud_optics.tau(icol,ilay,igpt) = cld_tau_gpt(icol,ilay,igpt);
    });
    cloud_optics.increment(combined_optics);

    // Call LW flux driver
    FluxesBroadband fluxes_allsky;
    fluxes_allsky.flux_up = allsky_flux_up;
    fluxes_allsky.flux_dn = allsky_flux_dn;
    fluxes_allsky.flux_net = allsky_flux_net;
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, combined_optics, top_at_1, lw_sources, emis_sfc, fluxes_allsky);
}
