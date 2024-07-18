#include <mam4xx/mam4.hpp>

namespace scream::impl {

using mam4::utils::min_max_bound;

using HostView1D    = haero::DeviceType::view_1d<Real>::HostMirror;
using HostViewInt1D = haero::DeviceType::view_1d<int>::HostMirror;

//-------------------------------------------------------------------------
//                    Reading the photolysis table
//-------------------------------------------------------------------------
// This logic is currently implemented using serial NetCDF calls for
// clarity of purpose. We should probably read the data for the photolysis
// table using SCREAM's SCORPIO interface instead, but I wanted to make
// clear what we're trying to do in terms of "elementary" operations first.

// ON HOST (MPI root rank only), reads the dimension of a NetCDF variable from
// the file with the given ID
int nc_dimension(const char *file, int nc_id, const char *dim_name) {
  int dim_id;
  int result = nc_inq_dimid(nc_id, dim_name, &dim_id);
  EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't fetch " << dim_name <<
    " dimension ID from NetCDF file '" << file << "'\n");
  size_t dim;
  result = nc_inq_dimlen(nc_id, dim_id, &dim);
  EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't fetch " << dim_name <<
    " dimension from NetCDF file '" << file << "'\n");
  return static_cast<int>(dim);
}

// ON HOST (MPI root rank only), reads data from the given NetCDF variable from
// the file with the given ID into the given Kokkos host View
template <typename V>
void read_nc_var(const char *file, int nc_id, const char *var_name, V host_view) {
  int var_id;
  int result = nc_inq_varid(nc_id, var_name, &var_id);
  EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't fetch ID for variable '" << var_name <<
    "' from NetCDF file '" << file << "'\n");
  result = nc_get_var(nc_id, var_id, host_view.data());
  EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't read data for variable '" << var_name <<
    "' from NetCDF file '" << file << "'\n");
}

// ON HOST (MPI root rank only), reads data from the NetCDF variable with the
// given ID, from the file with the given ID, into the given Kokkos host View
template <typename V>
void read_nc_var(const char *file, int nc_id, int var_id, V host_view) {
  int result = nc_get_var(nc_id, var_id, host_view.data());
  EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't read data for variable with ID " <<
    var_id << " from NetCDF file '" << file << "'\n");
}

// ON HOST (MPI root only), sets the lng_indexer and pht_alias_mult_1 host views
// according to parameters in our (hardwired) chemical mechanism
void set_lng_indexer_and_pht_alias_mult_1(const char *file, int nc_id,
                                          HostViewInt1D lng_indexer,
                                          HostView1D pht_alias_mult_1) {
  // NOTE: it seems that the chemical mechanism we're using
  // NOTE: 1. sets pht_alias_lst to a blank string [1]
  // NOTE: 2. sets pht_alias_mult_1 to 1.0 [1]
  // NOTE: 3. sets rxt_tag_lst to ['jh2o2', 'usr_HO2_HO2', 'usr_SO2_OH', 'usr_DMS_OH'] [2]
  // NOTE: References:
  // NOTE: [1] (https://github.com/eagles-project/e3sm_mam4_refactor/blob/refactor-maint-2.0/components/eam/src/chemistry/pp_linoz_mam4_resus_mom_soag/mo_sim_dat.F90#L117)
  // NOTE: [2] (https://github.com/eagles-project/e3sm_mam4_refactor/blob/refactor-maint-2.0/components/eam/src/chemistry/pp_linoz_mam4_resus_mom_soag/mo_sim_dat.F90#L99)

  // populate lng_indexer (see https://github.com/eagles-project/e3sm_mam4_refactor/blob/refactor-maint-2.0/components/eam/src/chemistry/mozart/mo_jlong.F90#L180)
  static const char *var_names[4] = {"jh2o2", "usr_HO2_HO2", "usr_SO2_OH", "usr_DMS_OH"};
  for (int m = 0; m < mam4::mo_photo::phtcnt; ++m) {
    int var_id;
    int result = nc_inq_varid(nc_id, var_names[m], &var_id);
    EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't fetch ID for variable '"
      << var_names[m] << "' from NetCDF file '" << file << "'\n");
    lng_indexer(m) = var_id;
  }

  // set pht_alias_mult_1 to 1
  Kokkos::deep_copy(pht_alias_mult_1, 1.0);
}

// ON HOST (MPI root only), populates the etfphot view using rebinned
// solar data from our solar_data_file
void populate_etfphot(HostView1D we, HostView1D etfphot) {
  // FIXME: It looks like EAM is relying on a piece of infrastructure that
  // FIXME: we just don't have in EAMxx (eam/src/chemistry/utils/solar_data.F90).
  // FIXME: I have no idea whether EAMxx has a plan for supporting this
  // FIXME: solar irradiance / photon flux data, and I'm not going to recreate
  // FIXME: that capability here. So this is an unplugged hole.
  // FIXME:
  // FIXME: If we are going to do this the way EAM does it, the relevant logic
  // FIXME: is the call to rebin() in eam/src/chemistry/mozart/mo_jlong.F90,
  // FIXME: around line 104.

  // FIXME: zero the photon flux for now
  Kokkos::deep_copy(etfphot, 0);
}
// This version uses scream_scorpio_interface to read netcdf files.
mam4::mo_photo::PhotoTableData read_photo_table(const std::string& rsf_file,
                      const std::string& xs_long_file) {


// set up the lng_indexer and pht_alias_mult_1 views based on our
// (hardwired) chemical mechanism
HostViewInt1D lng_indexer_h("lng_indexer(host)", mam4::mo_photo::phtcnt);



int nw, nump, numsza, numcolo3, numalb, nt, np_xs; // table dimensions
scorpio::register_file(rsf_file,scorpio::Read);
// read and broadcast dimension data
nump     = scorpio::get_dimlen(rsf_file,"numz");
numsza   = scorpio::get_dimlen(rsf_file, "numsza");
numalb   = scorpio::get_dimlen(rsf_file, "numalb");
numcolo3 = scorpio::get_dimlen(rsf_file, "numcolo3fact");

scorpio::register_file(xs_long_file,scorpio::Read);
nt       = scorpio::get_dimlen(xs_long_file, "numtemp");
nw       = scorpio::get_dimlen(xs_long_file, "numwl");
np_xs    = scorpio::get_dimlen(xs_long_file, "numprs");

//FIXME: hard-coded for only one photo reaction.
std::string rxt_names[1] = {"jh2o2"};
int numj = 1;
lng_indexer_h(0)=0;
// allocate the photolysis table
auto table = mam4::mo_photo::create_photo_table_data(nw, nt, np_xs, numj,
                                                       nump, numsza, numcolo3,
                                                       numalb);

// allocate host views for table data
auto rsf_tab_h = Kokkos::create_mirror_view(table.rsf_tab);
auto xsqy_h = Kokkos::create_mirror_view(table.xsqy);
auto sza_h = Kokkos::create_mirror_view(table.sza);
auto alb_h = Kokkos::create_mirror_view(table.alb);
auto press_h = Kokkos::create_mirror_view(table.press);
auto colo3_h = Kokkos::create_mirror_view(table.colo3);
auto o3rat_h = Kokkos::create_mirror_view(table.o3rat);
auto etfphot_h = Kokkos::create_mirror_view(table.etfphot);
auto prs_h = Kokkos::create_mirror_view(table.prs);

// read file data into our host views
scorpio::read_var(rsf_file, "pm", press_h.data());
scorpio::read_var(rsf_file, "sza", sza_h.data());
scorpio::read_var(rsf_file, "alb", alb_h.data());
scorpio::read_var(rsf_file, "colo3fact", o3rat_h.data());
scorpio::read_var(rsf_file, "colo3", colo3_h.data());
// it produces an error.
scorpio::read_var(rsf_file, "RSF", rsf_tab_h.data());
scorpio::read_var(xs_long_file, "pressure", prs_h.data());

// read xsqy data (using lng_indexer_h for the first index)
//FIXME: hard-coded for only one photo reaction.
for (int m = 0; m < mam4::mo_photo::phtcnt; ++m) {
  auto xsqy_ndx_h = ekat::subview(xsqy_h, m);
  scorpio::read_var(xs_long_file, rxt_names[m], xsqy_h.data());
}

// populate etfphot by rebinning solar data
HostView1D wc_h("wc", nw), wlintv_h("wlintv", nw), we_h("we", nw+1);

scorpio::read_var(rsf_file, "wc", wc_h.data());
scorpio::read_var(rsf_file, "wlintv", wlintv_h.data());
for (int i = 0; i < nw; ++i) {
      we_h(i) = wc_h(i) - 0.5 * wlintv_h(i);
}
we_h(nw) = wc_h(nw-1) - 0.5 * wlintv_h(nw-1);
populate_etfphot(we_h, etfphot_h);
scorpio::release_file(rsf_file);
scorpio::release_file(xs_long_file);

// copy host photolysis table into place on device
Kokkos::deep_copy(table.rsf_tab,          rsf_tab_h);
Kokkos::deep_copy(table.xsqy,             xsqy_h);
Kokkos::deep_copy(table.sza,              sza_h);
Kokkos::deep_copy(table.alb,              alb_h);
Kokkos::deep_copy(table.press,            press_h);
Kokkos::deep_copy(table.colo3,            colo3_h);
Kokkos::deep_copy(table.o3rat,            o3rat_h);
Kokkos::deep_copy(table.etfphot,          etfphot_h);
Kokkos::deep_copy(table.prs,              prs_h);
// set pht_alias_mult_1 to 1
Kokkos::deep_copy(table.pht_alias_mult_1, 1.0);
Kokkos::deep_copy(table.lng_indexer,      lng_indexer_h);

  // compute gradients (on device)
  Kokkos::parallel_for("del_p", nump-1, KOKKOS_LAMBDA(int i) {
    table.del_p(i) = 1.0/::abs(table.press(i)- table.press(i+1));
  });
  Kokkos::parallel_for("del_sza", numsza-1, KOKKOS_LAMBDA(int i) {
    table.del_sza(i) = 1.0/(table.sza(i+1) - table.sza(i));
  });
  Kokkos::parallel_for("del_alb", numalb-1, KOKKOS_LAMBDA(int i) {
    table.del_alb(i) = 1.0/(table.alb(i+1) - table.alb(i));
  });
  Kokkos::parallel_for("del_o3rat", numcolo3-1, KOKKOS_LAMBDA(int i) {
    table.del_o3rat(i) = 1.0/(table.o3rat(i+1) - table.o3rat(i));
  });
  Kokkos::parallel_for("dprs", np_xs-1, KOKKOS_LAMBDA(int i) {
    table.dprs(i) = 1.0/(table.prs(i) - table.prs(i+1));
  });

return table;
}

#if 0
// ON HOST, reads the photolysis table (used for gas phase chemistry) from the
// files with the given names
mam4::mo_photo::PhotoTableData read_photo_table(const ekat::Comm& comm,
                                                const char *rsf_file,
                                                const char* xs_long_file) {
  // NOTE: at the time of development, SCREAM's SCORPIO interface seems intended
  // NOTE: for domain-decomposed grid data. The files we're reading here are not
  // NOTE: spatial data, and should be the same everywhere, so we read them
  // NOTE: using serial NetCDF calls on MPI rank 0 and broadcast to other ranks.
  const int mpi_root = 0;
  int rsf_id, xs_long_id; // NetCDF file IDs (used only on MPI root)
  int nw, nump, numsza, numcolo3, numalb, nt, np_xs; // table dimensions
  if (comm.rank() == mpi_root) { // read dimension data from files and broadcast
    // open files
    int result = nc_open(rsf_file, NC_NOWRITE, &rsf_id);
    EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't open rsf_file '" << rsf_file << "'\n");
    result = nc_open(xs_long_file, NC_NOWRITE, &xs_long_id);
    EKAT_REQUIRE_MSG(result == 0, "Error! Couldn't open xs_long_file '" << xs_long_file << "'\n");

    // read and broadcast dimension data
    nump     = nc_dimension(rsf_file, rsf_id, "numz");
    numsza   = nc_dimension(rsf_file, rsf_id, "numsza");
    numalb   = nc_dimension(rsf_file, rsf_id, "numalb");
    numcolo3 = nc_dimension(rsf_file, rsf_id, "numcolo3fact");
    nt       = nc_dimension(xs_long_file, xs_long_id, "numtemp");
    nw       = nc_dimension(xs_long_file, xs_long_id, "numwl");
    np_xs    = nc_dimension(xs_long_file, xs_long_id, "numprs");
    std::cout <<nump<<" : nump\n";

    int dim_data[7] = {nump, numsza, numcolo3, numalb, nt, nw, np_xs};
    comm.broadcast(dim_data, 7, mpi_root);
  } else { // receive broadcasted dimension data from root rank
    int dim_data[7];
    comm.broadcast(dim_data, 7, mpi_root);
    nump     = dim_data[0];
    numsza   = dim_data[1];
    numcolo3 = dim_data[2];
    numalb   = dim_data[3];
    nt       = dim_data[4];
    nw       = dim_data[5];
    np_xs    = dim_data[6];
  }

  // set up the lng_indexer and pht_alias_mult_1 views based on our
  // (hardwired) chemical mechanism
  HostViewInt1D lng_indexer_h("lng_indexer(host)", mam4::mo_photo::phtcnt);
  HostView1D pht_alias_mult_1_h("pht_alias_mult_1(host)", 2);
  if (comm.rank() == mpi_root) {
    set_lng_indexer_and_pht_alias_mult_1(xs_long_file, xs_long_id,
                                         lng_indexer_h, pht_alias_mult_1_h);
  }
  comm.broadcast(lng_indexer_h.data(),      mam4::mo_photo::phtcnt, mpi_root);
  comm.broadcast(pht_alias_mult_1_h.data(), 2,                      mpi_root);

  // compute the size of the foremost dimension of xsqy using lng_indexer
  int numj = 0;
  for (int m = 0; m < mam4::mo_photo::phtcnt; ++m) {
    if (lng_indexer_h(m) > 0) {
      for (int mm = 0; mm < m; ++mm) {
        if (lng_indexer_h(mm) == lng_indexer_h(m)) {
          break;
        }
        ++numj;
      }
    }
  }

  // allocate the photolysis table
  auto table = mam4::mo_photo::create_photo_table_data(nw, nt, np_xs, numj,
                                                       nump, numsza, numcolo3,
                                                       numalb);

  // allocate host views for table data
  auto rsf_tab_h = Kokkos::create_mirror_view(table.rsf_tab);
  auto xsqy_h = Kokkos::create_mirror_view(table.xsqy);
  auto sza_h = Kokkos::create_mirror_view(table.sza);
  auto alb_h = Kokkos::create_mirror_view(table.alb);
  auto press_h = Kokkos::create_mirror_view(table.press);
  auto colo3_h = Kokkos::create_mirror_view(table.colo3);
  auto o3rat_h = Kokkos::create_mirror_view(table.o3rat);
  auto etfphot_h = Kokkos::create_mirror_view(table.etfphot);
  auto prs_h = Kokkos::create_mirror_view(table.prs);

  if (comm.rank() == mpi_root) { // read data from files and broadcast
    // read file data into our host views
    read_nc_var(rsf_file, rsf_id, "pm", press_h);
    read_nc_var(rsf_file, rsf_id, "sza", sza_h);
    read_nc_var(rsf_file, rsf_id, "alb", alb_h);
    read_nc_var(rsf_file, rsf_id, "colo3fact", o3rat_h);
    read_nc_var(rsf_file, rsf_id, "colo3", colo3_h);
    read_nc_var(rsf_file, rsf_id, "RSF", rsf_tab_h);

    read_nc_var(xs_long_file, xs_long_id, "pressure", prs_h);

    // read xsqy data (using lng_indexer_h for the first index)
    int ndx = 0;
    for (int m = 0; m < mam4::mo_photo::phtcnt; ++m) {
      if (lng_indexer_h(m) > 0) {
        auto xsqy_ndx_h = ekat::subview(xsqy_h, ndx);
        read_nc_var(xs_long_file, xs_long_id, lng_indexer_h(m), xsqy_ndx_h);
        ++ndx;
      }
    }

    // populate etfphot by rebinning solar data
    HostView1D wc_h("wc", nw), wlintv_h("wlintv", nw), we_h("we", nw+1);
    read_nc_var(rsf_file, rsf_id, "wc", wc_h);
    read_nc_var(rsf_file, rsf_id, "wlintv", wlintv_h);
    for (int i = 0; i < nw; ++i) {
      we_h(i) = wc_h(i) - 0.5 * wlintv_h(i);
    }
    we_h(nw) = wc_h(nw-1) - 0.5 * wlintv_h(nw-1);
    populate_etfphot(we_h, etfphot_h);

    // close the files
    nc_close(rsf_id);
    nc_close(xs_long_id);
  }

  // broadcast host views from MPI root to others
  comm.broadcast(rsf_tab_h.data(), nw*numalb*numcolo3*numsza*nump, mpi_root);
  comm.broadcast(xsqy_h.data(),    numj*nw*nt*np_xs,               mpi_root);
  comm.broadcast(sza_h.data(),     numsza,                         mpi_root);
  comm.broadcast(alb_h.data(),     numalb,                         mpi_root);
  comm.broadcast(press_h.data(),   nump,                           mpi_root);
  comm.broadcast(o3rat_h.data(),   numcolo3,                       mpi_root);
  comm.broadcast(colo3_h.data(),   nump,                           mpi_root);
  comm.broadcast(etfphot_h.data(), nw,                             mpi_root);
  comm.broadcast(prs_h.data(),     np_xs,                          mpi_root);

  // copy host photolysis table into place on device
  Kokkos::deep_copy(table.rsf_tab,          rsf_tab_h);
  Kokkos::deep_copy(table.xsqy,             xsqy_h);
  Kokkos::deep_copy(table.sza,              sza_h);
  Kokkos::deep_copy(table.alb,              alb_h);
  Kokkos::deep_copy(table.press,            press_h);
  Kokkos::deep_copy(table.colo3,            colo3_h);
  Kokkos::deep_copy(table.o3rat,            o3rat_h);
  Kokkos::deep_copy(table.etfphot,          etfphot_h);
  Kokkos::deep_copy(table.prs,              prs_h);
  Kokkos::deep_copy(table.pht_alias_mult_1, pht_alias_mult_1_h);
  Kokkos::deep_copy(table.lng_indexer,      lng_indexer_h);

  // compute gradients (on device)
  Kokkos::parallel_for("del_p", nump-1, KOKKOS_LAMBDA(int i) {
    table.del_p(i) = 1.0/::abs(table.press(i)- table.press(i+1));
  });
  Kokkos::parallel_for("del_sza", numsza-1, KOKKOS_LAMBDA(int i) {
    table.del_sza(i) = 1.0/(table.sza(i+1) - table.sza(i));
  });
  Kokkos::parallel_for("del_alb", numalb-1, KOKKOS_LAMBDA(int i) {
    table.del_alb(i) = 1.0/(table.alb(i+1) - table.alb(i));
  });
  Kokkos::parallel_for("del_o3rat", numcolo3-1, KOKKOS_LAMBDA(int i) {
    table.del_o3rat(i) = 1.0/(table.o3rat(i+1) - table.o3rat(i));
  });
  Kokkos::parallel_for("dprs", np_xs-1, KOKKOS_LAMBDA(int i) {
    table.dprs(i) = 1.0/(table.prs(i) - table.prs(i+1));
  });

  return table;
}
#endif
// performs gas phase chemistry calculations on a single level of a single
// atmospheric column
KOKKOS_INLINE_FUNCTION
void gas_phase_chemistry(Real zm, Real zi, Real phis, Real temp, Real pmid, Real pdel, Real dt,
                         const Real photo_rates[mam4::mo_photo::phtcnt], // in
                         Real q[mam4::gas_chemistry::gas_pcnst], // VMRs, inout
                         Real invariants[mam4::gas_chemistry::nfs]) { // out
  // constexpr Real rga = 1.0/haero::Constants::gravity;
  // constexpr Real m2km = 0.01; // converts m -> km

  // The following things are chemical mechanism dependent! See mam4xx/src/mam4xx/gas_chem_mechanism.hpp)
  constexpr int gas_pcnst = mam4::gas_chemistry::gas_pcnst; // number of gas phase species
  constexpr int rxntot = mam4::gas_chemistry::rxntot;       // number of chemical reactions
  constexpr int indexm = mam4::gas_chemistry::indexm;       // index of total atm density in invariant array 

  constexpr int phtcnt = mam4::mo_photo::phtcnt; // number of photolysis reactions

  constexpr int itermax = mam4::gas_chemistry::itermax;
  constexpr int clscnt4 = mam4::gas_chemistry::clscnt4;
  constexpr int nfs = mam4::gas_chemistry::nfs;

  // NOTE: vvv these arrays were copied from mam4xx/gas_chem_mechanism.hpp vvv
  constexpr int permute_4[gas_pcnst] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                        10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                        20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
  constexpr int clsmap_4[gas_pcnst] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

  // These indices for species are fixed by the chemical mechanism
  // std::string solsym[] = {"O3", "H2O2", "H2SO4", "SO2", "DMS", "SOAG",
  //                         "so4_a1", "pom_a1", "soa_a1", "bc_a1", "dst_a1",
  //                         "ncl_a1", "mom_a1", "num_a1", "so4_a2", "soa_a2",
  //                         "ncl_a2", "mom_a2", "num_a2", "dst_a3", "ncl_a3",
  //                         "so4_a3", "bc_a3", "pom_a3", "soa_a3", "mom_a3",
  //                         "num_a3", "pom_a4", "bc_a4", "mom_a4", "num_a4"};
  constexpr int ndx_h2so4 = 2;
  // std::string extfrc_list[] = {"SO2", "so4_a1", "so4_a2", "pom_a4", "bc_a4",
  //                              "num_a1", "num_a2", "num_a3", "num_a4", "SOAG"};
  constexpr int synoz_ndx = -1;

  // fetch the zenith angle (not its cosine!) in degrees for this column.
  // FIXME: For now, we fix the zenith angle. At length, we need to compute it
  // FIXME: from EAMxx's current set of orbital parameters, which requires some
  // FIXME: conversation with the EAMxx team.

  // xform geopotential height from m to km and pressure from Pa to mb
  // Real zsurf = rga * phis;
  // Real zmid = m2km * (zm + zsurf);

  // ... compute the column's invariants
  // Real h2ovmr = q[0];
  // setinv(invariants, temp, h2ovmr, q, pmid); FIXME: not ported yet
  for (int i = 0; i < nfs; ++i) {
    invariants[i] = 0.1;
  }

  // ... set rates for "tabular" and user specified reactions
  Real reaction_rates[rxntot];
  mam4::gas_chemistry::setrxt(reaction_rates, temp);

  // set reaction rates based on chemical invariants
  // (indices (ndxes?) are taken from mam4 validation data and translated from
  // 1-based indices to 0-based indices)
  int usr_HO2_HO2_ndx = 1, usr_DMS_OH_ndx = 5,
      usr_SO2_OH_ndx = 3, inv_h2o_ndx = 3;
  mam4::gas_chemistry::usrrxt(reaction_rates, temp, invariants, invariants[indexm],
                              usr_HO2_HO2_ndx, usr_DMS_OH_ndx,
                              usr_SO2_OH_ndx, inv_h2o_ndx);
  mam4::gas_chemistry::adjrxt(reaction_rates, invariants, invariants[indexm]);

  //===================================
  // Photolysis rates at time = t(n+1)
  //===================================


  // ... Form the washout rates
  Real het_rates[gas_pcnst];
  // FIXME: not ported yet
  //sethet(het_rates, pmid, zmid, phis, temp, cmfdqr, prain, nevapr, delt,
  //       invariants[indexm], q);


  // save h2so4 before gas phase chem (for later new particle nucleation)
  Real del_h2so4_gasprod = q[ndx_h2so4];

  //===========================
  // Class solution algorithms
  //===========================

  // copy photolysis rates into reaction_rates (assumes photolysis rates come first)
  for (int i = 0; i < phtcnt; ++i) {
    reaction_rates[i] = photo_rates[i];
  }

  // ... solve for "Implicit" species
  bool factor[itermax];
  for (int i = 0; i < itermax; ++i) {
    factor[i] = true;
  }

  // initialize error tolerances
  Real epsilon[clscnt4];
  mam4::gas_chemistry::imp_slv_inti(epsilon);

  // solve chemical system implicitly
  Real prod_out[clscnt4], loss_out[clscnt4];
  mam4::gas_chemistry::imp_sol(q, reaction_rates, het_rates, dt,
    permute_4, clsmap_4, factor, epsilon, prod_out, loss_out);

  // save h2so4 change by gas phase chem (for later new particle nucleation)
  if (ndx_h2so4 > 0) {
    del_h2so4_gasprod = q[ndx_h2so4] - del_h2so4_gasprod;
  }
}

} // namespace scream::impl
