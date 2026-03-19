#include "physics/rrtmgp/eamxx_rrtmgp_interface.hpp"
#include "physics/rrtmgp/rrtmgp_test_utils.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/core/eamxx_session.hpp"

// From RRTMGP submodule
#include <cpp/rrtmgp/mo_gas_concentrations.h>
#include <mo_garand_atmos_io.h>

#include <ekat_logger.hpp>

#include <iostream>
#include <cmath>

using namespace scream;

// Names of input files we will need.
// NOTE: use full spectral resolution absorption data for consistency with reference problem
std::string coefficients_file_sw = SCREAM_DATA_DIR "/init/rrtmgp-data-sw-g224-2018-12-04.nc";
std::string coefficients_file_lw = SCREAM_DATA_DIR "/init/rrtmgp-data-lw-g256-2018-12-04.nc";
std::string cloud_optics_file_sw = SCREAM_DATA_DIR "/init/rrtmgp-cloud-optics-coeffs-sw.nc";
std::string cloud_optics_file_lw = SCREAM_DATA_DIR "/init/rrtmgp-cloud-optics-coeffs-lw.nc";

int main (int argc, char** argv) {
    MPI_Init(&argc,&argv);

    using namespace ekat::logger;
    using logger_t = Logger<LogNoFile,LogRootRank>;

    using layout_t    = Kokkos::LayoutLeft;
    using interface_t = scream::rrtmgp::rrtmgp_interface<Real, layout_t>;
    using utils_t     = rrtmgpTest::rrtmgp_test_utils<Real, layout_t>;
    using gas_concs_t = typename interface_t::gas_concs_t;
    using r1d         = typename interface_t::real1dk;
    using r2d         = typename interface_t::real2dk;
    using r3d         = typename interface_t::real3dk;
    using MDRP        = typename interface_t::MDRP;

    ekat::Comm comm(MPI_COMM_WORLD);
    auto logger = std::make_shared<logger_t>("",LogLevel::info,comm);

    // Get filenames from command line
    if (argc < 3) {
      std::string msg = "Missing required inputs. Usage:\n";
      msg += argv[0];
      msg += " inputfile baseline\n";
      logger->error(msg);
      return 1;
    }
    std::string inputfile(argv[1]);
    std::string baseline(argv[2]);

    // Initialize kernel launcher
    scream::init_kls();

    // Get reference fluxes from input file; do this here so we can get ncol dimension
    r2d sw_flux_up_ref;
    r2d sw_flux_dn_ref;
    r2d sw_flux_dn_dir_ref;
    r2d lw_flux_up_ref;
    r2d lw_flux_dn_ref;
    logger->info("read_fluxes...");
    utils_t::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

    // Get dimension sizes
    const int ncol = sw_flux_up_ref.extent(0);
    const int nlev = sw_flux_up_ref.extent(1);
    const int nlay = nlev - 1;

    // Read in dummy Garand atmosphere; if this were an actual model simulation,
    // these would be passed as inputs to the driver
    // NOTE: set ncol to size of col_flx dimension in the input file. This is so
    // that we can compare to the reference data provided in that file. Note that
    // this will copy the first column of the input data (the first profile) ncol
    // times. We will then fill some fraction of these columns with clouds for
    // the test problem.
    r2d p_lay ("p_lay", ncol, nlay);
    r2d t_lay ("t_lay", ncol, nlay);
    r2d p_lev ("p_lev", ncol, nlay+1);
    r2d t_lev ("t_lev", ncol, nlay+1);
    r2d col_dry;
    gas_concs_t gas_concs;
    read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

    // Initialize the RRTMGP interface; this will read in the k-distribution
    // data that contains information about absorption coefficients for gases
    logger->info("rrtmgp_initialize...");
    interface_t::rrtmgp_initialize(gas_concs, coefficients_file_sw, coefficients_file_lw, cloud_optics_file_sw, cloud_optics_file_lw, logger, 2.0);

    // Setup dummy all-sky problem
    r1d sfc_alb_dir_vis ("sfc_alb_dir_vis", ncol);
    r1d sfc_alb_dir_nir ("sfc_alb_dir_nir", ncol);
    r1d sfc_alb_dif_vis ("sfc_alb_dif_vis", ncol);
    r1d sfc_alb_dif_nir ("sfc_alb_dif_nir", ncol);
    r1d mu0 ("mu0", ncol);
    r2d lwp ("lwp", ncol, nlay);
    r2d iwp ("iwp", ncol, nlay);
    r2d rel ("rel", ncol, nlay);
    r2d rei ("rei", ncol, nlay);
    r2d cld ("cld", ncol, nlay);
    utils_t::dummy_atmos(
        inputfile, ncol, p_lay, t_lay,
        sfc_alb_dir_vis, sfc_alb_dir_nir,
        sfc_alb_dif_vis, sfc_alb_dif_nir,
        mu0,
        lwp, iwp, rel, rei, cld
    );

    // Setup flux outputs; In a real model run, the fluxes would be
    // input/outputs into the driver (persisting between calls), and
    // we would just have to setup the pointers to them in the
    // FluxesBroadband object
    const auto nswbands = interface_t::k_dist_sw_k->get_nband();
    const auto nlwbands = interface_t::k_dist_lw_k->get_nband();
    r2d sw_flux_up ("sw_flux_up" , ncol, nlay+1);
    r2d sw_flux_dn ("sw_flux_dn" , ncol, nlay+1);
    r2d sw_flux_dn_dir("sw_flux_dn_dir", ncol, nlay+1);
    r2d lw_flux_up ("lw_flux_up" , ncol, nlay+1);
    r2d lw_flux_dn ("lw_flux_dn" , ncol, nlay+1);
    r2d sw_clnclrsky_flux_up ("sw_clnclrsky_flux_up" , ncol, nlay+1);
    r2d sw_clnclrsky_flux_dn ("sw_clnclrsky_flux_dn" , ncol, nlay+1);
    r2d sw_clnclrsky_flux_dn_dir("sw_clnclrsky_flux_dn_dir", ncol, nlay+1);
    r2d sw_clrsky_flux_up ("sw_clrsky_flux_up" , ncol, nlay+1);
    r2d sw_clrsky_flux_dn ("sw_clrsky_flux_dn" , ncol, nlay+1);
    r2d sw_clrsky_flux_dn_dir("sw_clrsky_flux_dn_dir", ncol, nlay+1);
    r2d sw_clnsky_flux_up ("sw_clnsky_flux_up" , ncol, nlay+1);
    r2d sw_clnsky_flux_dn ("sw_clnsky_flux_dn" , ncol, nlay+1);
    r2d sw_clnsky_flux_dn_dir("sw_clnsky_flux_dn_dir", ncol, nlay+1);
    r2d lw_clnclrsky_flux_up ("lw_clnclrsky_flux_up" , ncol, nlay+1);
    r2d lw_clnclrsky_flux_dn ("lw_clnclrsky_flux_dn" , ncol, nlay+1);
    r2d lw_clrsky_flux_up ("lw_clrsky_flux_up" , ncol, nlay+1);
    r2d lw_clrsky_flux_dn ("lw_clrsky_flux_dn" , ncol, nlay+1);
    r2d lw_clnsky_flux_up ("lw_clnsky_flux_up" , ncol, nlay+1);
    r2d lw_clnsky_flux_dn ("lw_clnsky_flux_dn" , ncol, nlay+1);
    r3d sw_bnd_flux_up ("sw_bnd_flux_up" , ncol, nlay+1, nswbands);
    r3d sw_bnd_flux_dn ("sw_bnd_flux_dn" , ncol, nlay+1, nswbands);
    r3d sw_bnd_flux_dir("sw_bnd_flux_dir", ncol, nlay+1, nswbands);
    r3d lw_bnd_flux_up ("lw_bnd_flux_up" ,ncol, nlay+1, nlwbands);
    r3d lw_bnd_flux_dn ("lw_bnd_flux_dn" ,ncol, nlay+1, nlwbands);

    // Compute band-by-band surface_albedos.
    r2d sfc_alb_dir("sfc_alb_dir", ncol, nswbands);
    r2d sfc_alb_dif("sfc_alb_dif", ncol, nswbands);
    interface_t::compute_band_by_band_surface_albedos(
      ncol, nswbands,
      sfc_alb_dir_vis, sfc_alb_dir_nir,
      sfc_alb_dif_vis, sfc_alb_dif_nir,
      sfc_alb_dir, sfc_alb_dif);

    // Setup some dummy aerosol optical properties
    auto aer_tau_sw = r3d("aer_tau_sw", ncol, nlay, nswbands);
    auto aer_ssa_sw = r3d("aer_ssa_sw", ncol, nlay, nswbands);
    auto aer_asm_sw = r3d("aer_asm_sw", ncol, nlay, nswbands);
    auto aer_tau_lw = r3d("aer_tau_lw", ncol, nlay, nlwbands);
    Kokkos::parallel_for( MDRP::template get<3>({nswbands,nlay,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilay, int icol) {
        aer_tau_sw(icol,ilay,ibnd) = 0;
        aer_ssa_sw(icol,ilay,ibnd) = 0;
        aer_asm_sw(icol,ilay,ibnd) = 0;
    });
    Kokkos::parallel_for( MDRP::template get<3>({nlwbands,nlay,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilay, int icol) {
        aer_tau_lw(icol,ilay,ibnd) = 0;
    });

    // These are returned as outputs now from rrtmgp_main
    // TODO: provide as inputs consistent with how aerosol is treated?
    const auto nswgpts = interface_t::k_dist_sw_k->get_ngpt();
    const auto nlwgpts = interface_t::k_dist_lw_k->get_ngpt();
    auto cld_tau_sw_bnd = r3d("cld_tau_sw_bnd", ncol, nlay, nswbands);
    auto cld_tau_lw_bnd = r3d("cld_tau_lw_bnd", ncol, nlay, nlwbands);
    auto cld_tau_sw = r3d("cld_tau_sw", ncol, nlay, nswgpts);
    auto cld_tau_lw = r3d("cld_tau_lw", ncol, nlay, nlwgpts);

    // Run RRTMGP standalone codes and compare with AD run
    // Do something interesting here...
    // NOTE: these will get replaced with AD stuff that handles these
    logger->info("rrtmgp_main...");
    const Real tsi_scaling = 1;
    interface_t::rrtmgp_main(
        ncol, nlay,
        p_lay, t_lay, p_lev, t_lev, gas_concs,
        sfc_alb_dir, sfc_alb_dif, mu0,
        lwp, iwp, rel, rei, cld,
        aer_tau_sw, aer_ssa_sw, aer_asm_sw, aer_tau_lw,
        cld_tau_sw_bnd, cld_tau_lw_bnd,
        cld_tau_sw, cld_tau_lw,  // outputs
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn,
        sw_clnclrsky_flux_up, sw_clnclrsky_flux_dn, sw_clnclrsky_flux_dn_dir,
        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir,
        sw_clnsky_flux_up, sw_clnsky_flux_dn, sw_clnsky_flux_dn_dir,
        lw_clnclrsky_flux_up, lw_clnclrsky_flux_dn,
        lw_clrsky_flux_up, lw_clrsky_flux_dn,
        lw_clnsky_flux_up, lw_clnsky_flux_dn,
        sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir,
        lw_bnd_flux_up, lw_bnd_flux_dn, tsi_scaling,
        logger
    );

    // Write fluxes
    logger->info("write_fluxes...");
    utils_t::write_fluxes(
        baseline,
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn
    );

    logger->info("cleaning up...");
    // Clean up from test; this is probably not necessary, these things
    // should be deallocated when they fall out of scope, but we should be
    // good citizens and clean up our mess.

    gas_concs.reset();
    interface_t::rrtmgp_finalize();
    scream::finalize_kls();

    MPI_Finalize();

    return 0;
}
