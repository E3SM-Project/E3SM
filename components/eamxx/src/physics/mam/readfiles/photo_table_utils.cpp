#include <mam4xx/mam4.hpp>

namespace scream::impl {

// number of photolysis reactions
using mam4::mo_photo::phtcnt;

using HostView1D    = haero::DeviceType::view_1d<Real>::HostMirror;
using HostViewInt1D = haero::DeviceType::view_1d<int>::HostMirror;

//-------------------------------------------------------------------------
//                    Reading the photolysis table
//-------------------------------------------------------------------------

// ON HOST (MPI root only), sets the lng_indexer and pht_alias_mult_1 host views
// according to parameters in our (hardwired) chemical mechanism

std::vector<Real> populate_etfphot_from_e3sm_case() {
  // We obtained these values from an e3sm simulations.
  // We should only use this function on Host.
  std::vector<Real> etfphot_data = {
      7.5691227E+11, 8.6525905E+11, 1.0355749E+12, 1.1846288E+12, 2.1524405E+12,
      3.2362584E+12, 3.7289849E+12, 4.4204330E+12, 4.6835350E+12, 6.1217728E+12,
      4.5575051E+12, 5.3491446E+12, 4.7016063E+12, 5.4281722E+12, 4.5023968E+12,
      6.8931981E+12, 6.2012647E+12, 6.1430771E+12, 5.7820385E+12, 7.6770646E+12,
      1.3966509E+13, 1.2105348E+13, 2.8588980E+13, 3.2160821E+13, 2.4978066E+13,
      2.7825401E+13, 2.3276451E+13, 3.6343684E+13, 6.1787886E+13, 7.8009914E+13,
      7.6440824E+13, 7.6291458E+13, 9.4645085E+13, 1.0124628E+14, 1.0354111E+14,
      1.0999650E+14, 1.0889946E+14, 1.1381912E+14, 1.3490042E+14, 1.5941519E+14,
      1.4983265E+14, 1.5184267E+14, 1.5991420E+14, 1.6976697E+14, 1.8771840E+14,
      1.6434367E+14, 1.8371960E+14, 2.1966369E+14, 1.9617879E+14, 2.2399700E+14,
      1.8429912E+14, 2.0129736E+14, 2.0541588E+14, 2.4334962E+14, 3.5077122E+14,
      3.4517894E+14, 3.5749668E+14, 3.6624304E+14, 3.4975113E+14, 3.5566025E+14,
      4.2825273E+14, 4.8406375E+14, 4.9511159E+14, 5.2695368E+14, 5.2401611E+14,
      5.0877746E+14, 4.8780853E+14};
  return etfphot_data;
}

// This version uses eamxx_scorpio_interface to read netcdf files.
mam4::mo_photo::PhotoTableData read_photo_table(
    const std::string &rsf_file, const std::string &xs_long_file) {
  // set up the lng_indexer and pht_alias_mult_1 views based on our
  // (hardwired) chemical mechanism
  HostViewInt1D lng_indexer_h("lng_indexer", phtcnt);

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

  // FIXME: hard-coded for only one photo reaction.
  std::string rxt_names[1] = {"jh2o2"};
  int numj                 = 1;
  lng_indexer_h(0)         = 0;
  // allocate the photolysis table
  auto table = mam4::mo_photo::create_photo_table_data(
      nw, nt, np_xs, numj, nump, numsza, numcolo3, numalb);

  // allocate host views for table data
  auto rsf_tab_h = Kokkos::create_mirror_view(table.rsf_tab);
  auto xsqy_h    = Kokkos::create_mirror_view(table.xsqy);
  auto sza_h     = Kokkos::create_mirror_view(table.sza);
  auto alb_h     = Kokkos::create_mirror_view(table.alb);
  auto press_h   = Kokkos::create_mirror_view(table.press);
  auto colo3_h   = Kokkos::create_mirror_view(table.colo3);
  auto o3rat_h   = Kokkos::create_mirror_view(table.o3rat);
  // auto etfphot_h = Kokkos::create_mirror_view(table.etfphot);
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
  // FIXME: hard-coded for only one photo reaction.
  for(int m = 0; m < phtcnt; ++m) {
    auto xsqy_ndx_h = ekat::subview(xsqy_h, m);
    scorpio::read_var(xs_long_file, rxt_names[m], xsqy_h.data());
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
        table.del_p(i) = 1.0 / haero::abs(table.press(i) - table.press(i + 1));
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

}  // namespace scream::impl
