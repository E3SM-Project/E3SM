#include <mam4xx/mam4.hpp>

namespace scream::impl {

// number of photolysis reactions
using mam4::mo_photo::phtcnt;

using HostView1D    = mam4::DeviceType::view_1d<Real>::host_mirror_type;
using HostView5D   = mam4::DeviceType::view<Real*****>::host_mirror_type;
using HostView3D   = mam4::DeviceType::view<Real***>::host_mirror_type;
using HostViewInt1D = mam4::DeviceType::view_1d<int>::host_mirror_type;

//-------------------------------------------------------------------------
//                    Reading the photolysis table
//-------------------------------------------------------------------------

// ON HOST (MPI root only), sets the lng_indexer and pht_alias_mult_1 host views
// according to parameters in our (hardwired) chemical mechanism

std::vector<Real> populate_etfphot_from_e3sm_case() {
  // We obtained these values from an e3sm simulations.
  // We should only use this function on Host.
  std::vector<Real> etfphot_data = {
      0.75691227453241626E+012, 0.86525904597344678E+012, 0.10355748678445210E+013, 0.11846288453143215E+013, 0.21524405047611838E+013,
      0.32362583636383438E+013, 0.37289849127086353E+013, 0.44204330229059023E+013, 0.46835350139683008E+013, 0.61217728454045146E+013,
      0.45575051094967529E+013, 0.53491446243876533E+013, 0.47016062694342764E+013, 0.54281722298247529E+013, 0.45023968313414365E+013,
      0.68931981401230361E+013, 0.62012647462481055E+013, 0.61430770669364131E+013, 0.57820384729408037E+013, 0.76770646262530391E+013,
      0.13966508541416857E+014, 0.12105347510143980E+014, 0.28588979654418141E+014, 0.32160820948665508E+014, 0.24978065543030500E+014,
      0.27825400776036188E+014, 0.23276451219415352E+014, 0.36343683716296695E+014, 0.61787885646314477E+014, 0.78009914475741344E+014,
      0.76440824240882500E+014, 0.76291457600771391E+014, 0.94645085080390984E+014, 0.10124627769922270E+015, 0.10354111421691689E+015,
      0.10999649606948711E+015, 0.10889946060495367E+015, 0.11381912455165878E+015, 0.13490042469475880E+015, 0.15941519351184984E+015,
      0.14983265369952531E+015, 0.15184267258496494E+015, 0.15991419729740088E+015, 0.16976696691694741E+015, 0.18771840486614825E+015,
      0.16434366552645634E+015, 0.18371960453616509E+015, 0.21966368981040753E+015, 0.19617878628663241E+015, 0.22399700059898819E+015,
      0.18429911731380941E+015, 0.20129735694980109E+015, 0.20541588491339825E+015, 0.24334961879677731E+015, 0.35077121778312700E+015,
      0.34517894220011569E+015, 0.35749668154179594E+015, 0.36624304237331069E+015, 0.34975112547690056E+015, 0.35566025203681831E+015,
      0.42825273260963562E+015, 0.48406375456076200E+015, 0.49511158653410975E+015, 0.52695367706176038E+015, 0.52401610578239200E+015,
      0.50877746346978994E+015, 0.48780852943692825E+015};
  return etfphot_data;
}

// This version uses eamxx_scorpio_interface to read netcdf files.
mam4::mo_photo::PhotoTableData read_photo_table(
    const std::string &rsf_file, const std::string &xs_long_file,
    const std::vector<std::string> &rxt_names, const int numj,
    const HostViewInt1D &lng_indexer_h) {

  EKAT_REQUIRE_MSG(numj > 0, "Error: read_photo_table requires numj > 0.\n");
  EKAT_REQUIRE_MSG(lng_indexer_h.extent_int(0) == phtcnt,
                   "Error: read_photo_table requires lng_indexer_h sized by phtcnt.\n");


  int nw, nump, numsza, numcolo3, numalb, nt, np_xs;  // table dimensions
  scorpio::register_file(rsf_file, scorpio::Read);
  // read and broadcast dimension data
  nump     = scorpio::get_dimlen(rsf_file, "numz");
  numsza   = scorpio::get_dimlen(rsf_file, "numsza");
  numalb   = scorpio::get_dimlen(rsf_file, "numalb");
  numcolo3 = scorpio::get_dimlen(rsf_file, "numcolo3fact");

  scorpio::register_file(xs_long_file, scorpio::Read);
  nt    = scorpio::get_dimlen(xs_long_file, "numtemp");
  nw    = scorpio::get_dimlen(xs_long_file, "numwl");
  np_xs = scorpio::get_dimlen(xs_long_file, "numprs");

  // allocate the photolysis table
  auto table = mam4::mo_photo::create_photo_table_data(
      nw, nt, np_xs, numj, nump, numsza, numcolo3, numalb);

  // allocate host views for table data
  HostView5D l_rsf_tab_h("rsf_tab_h",numalb,numcolo3,numsza,nump,nw);
  HostView3D l_xsqy_h("xsqy_h",np_xs,nt,nw);
  auto rsf_tab_h = Kokkos::create_mirror_view(table.rsf_tab);
  auto xsqy_h    = Kokkos::create_mirror_view(table.xsqy);
  auto sza_h     = Kokkos::create_mirror_view(table.sza);
  auto alb_h     = Kokkos::create_mirror_view(table.alb);
  auto press_h   = Kokkos::create_mirror_view(table.press);
  auto colo3_h   = Kokkos::create_mirror_view(table.colo3);
  auto o3rat_h   = Kokkos::create_mirror_view(table.o3rat);
  auto prs_h = Kokkos::create_mirror_view(table.prs);
  

  // read file data into our host views
  scorpio::read_var(rsf_file, "pm", press_h.data());
  scorpio::read_var(rsf_file, "sza", sza_h.data());
  scorpio::read_var(rsf_file, "alb", alb_h.data());
  scorpio::read_var(rsf_file, "colo3fact", o3rat_h.data());
  scorpio::read_var(rsf_file, "colo3", colo3_h.data());
  scorpio::read_var(rsf_file, "RSF", l_rsf_tab_h.data());
  scorpio::read_var(xs_long_file, "pressure", prs_h.data());

  // read xsqy data (using lng_indexer_h for the first index)
  using policy_t3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>, Kokkos::DefaultHostExecutionSpace>;
  for(int m = 0; m < numj; ++m) {
    scorpio::read_var(xs_long_file, rxt_names[m], l_xsqy_h.data());
    Kokkos::parallel_for("xsqy_h", 
    policy_t3({0, 0, 0}, {xsqy_h.extent(1), xsqy_h.extent(2), xsqy_h.extent(3)}),
    [&](const int i, const int j, const int k) {
        xsqy_h(m, i, j, k) = l_xsqy_h(k,j,i);
  });
  }
  // populate etfphot by rebinning solar data
  HostView1D wc_h("wc", nw), wlintv_h("wlintv", nw), we_h("we", nw + 1);

  scorpio::read_var(rsf_file, "wc", wc_h.data());
  scorpio::read_var(rsf_file, "wlintv", wlintv_h.data());
  for(int i = 0; i < nw; ++i) {
    we_h(i) = wc_h(i) - 0.5 * wlintv_h(i);
  }
  we_h(nw) = wc_h(nw - 1) - 0.5 * wlintv_h(nw - 1);
  // populate_etfphot(we_h, etfphot_h);
  // FIXME: etfphot_data is hard-coded.
  auto etfphot_data = populate_etfphot_from_e3sm_case();
  auto etfphot_h    = HostView1D((Real *)etfphot_data.data(), nw);

  using policy_t = Kokkos::MDRangePolicy<Kokkos::Rank<4>, Kokkos::DefaultHostExecutionSpace>;
  
  Kokkos::parallel_for("scale_rsf_tab", 
    policy_t({0, 0, 0, 0}, {rsf_tab_h.extent(1), rsf_tab_h.extent(2), rsf_tab_h.extent(3), rsf_tab_h.extent(4)}),
    [&](const int l, const int i, const int j, const int k) {
      for (int w = 0; w < nw; ++w) {
        rsf_tab_h(w,l, i, j, k) = l_rsf_tab_h(k,j,i,l,w)*wlintv_h(w);
      }
  });

  scorpio::release_file(rsf_file);
  scorpio::release_file(xs_long_file);

  // copy host photolysis table into place on device
  Kokkos::deep_copy(table.rsf_tab, rsf_tab_h);
  Kokkos::deep_copy(table.xsqy, xsqy_h);
  Kokkos::deep_copy(table.sza, sza_h);
  Kokkos::deep_copy(table.alb, alb_h);
  Kokkos::deep_copy(table.press, press_h);
  Kokkos::deep_copy(table.colo3, colo3_h);
  Kokkos::deep_copy(table.o3rat, o3rat_h);
  Kokkos::deep_copy(table.etfphot, etfphot_h);
  Kokkos::deep_copy(table.prs, prs_h);
  // set pht_alias_mult_1 to 1
  Kokkos::deep_copy(table.pht_alias_mult_1, 1.0);
  Kokkos::deep_copy(table.lng_indexer, lng_indexer_h);

  // compute gradients (on device)
  Kokkos::parallel_for(
      "del_p", nump - 1, KOKKOS_LAMBDA(int i) {
        table.del_p(i) = 1.0 / mam4::abs(table.press(i) - table.press(i + 1));
      });
  Kokkos::parallel_for(
      "del_sza", numsza - 1, KOKKOS_LAMBDA(int i) {
        table.del_sza(i) = 1.0 / (table.sza(i + 1) - table.sza(i));
      });
  Kokkos::parallel_for(
      "del_alb", numalb - 1, KOKKOS_LAMBDA(int i) {
        table.del_alb(i) = 1.0 / (table.alb(i + 1) - table.alb(i));
      });
  Kokkos::parallel_for(
      "del_o3rat", numcolo3 - 1, KOKKOS_LAMBDA(int i) {
        table.del_o3rat(i) = 1.0 / (table.o3rat(i + 1) - table.o3rat(i));
      });
  Kokkos::parallel_for(
      "dprs", np_xs - 1, KOKKOS_LAMBDA(int i) {
        table.dprs(i) = 1.0 / (table.prs(i) - table.prs(i + 1));
      });

  return table;
}

// MAM4xx E3SM v2 photolysis table reader.
mam4::mo_photo::PhotoTableData read_photo_table(
    const std::string &rsf_file, const std::string &xs_long_file) {
  
  HostViewInt1D lng_indexer_h("lng_indexer", phtcnt);
  std::vector<std::string> rxt_names = {"jh2o2"};
  int numj                 = 1;
  lng_indexer_h(0)         = 0;
  return read_photo_table(rsf_file, xs_long_file, rxt_names, numj, lng_indexer_h);
}

}  // namespace scream::impl
