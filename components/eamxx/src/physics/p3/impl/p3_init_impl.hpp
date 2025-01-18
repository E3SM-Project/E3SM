#ifndef P3_INIT_IMPL_HPP
#define P3_INIT_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

#include "ekat/util/ekat_file_utils.hpp"

#include <fstream>

namespace scream {
namespace p3 {

namespace {

template <typename S, typename IceT, typename CollT>
void read_ice_lookup_tables(const bool masterproc, const char* p3_lookup_base, const char* p3_version, IceT& ice_table_vals, CollT& collect_table_vals, int densize, int rimsize, int isize, int rcollsize)
{
  using DeviceIcetable = typename IceT::non_const_type;
  using DeviceColtable = typename CollT::non_const_type;

  const auto ice_table_vals_d     = DeviceIcetable("ice_table_vals");
  const auto collect_table_vals_d = DeviceColtable("collect_table_vals");

  const auto ice_table_vals_h    = Kokkos::create_mirror_view(ice_table_vals_d);
  const auto collect_table_vals_h = Kokkos::create_mirror_view(collect_table_vals_d);

  //
  // read in ice microphysics table into host views. We always read these as doubles.
  //

  std::string filename = std::string(p3_lookup_base) + std::string(p3_version);

  if (masterproc) {
    std::cout << "Reading ice lookup tables in file: " << filename << std::endl;
  }

  std::ifstream in(filename);

  // read header
  std::string version, version_val;
  in >> version >> version_val;
  EKAT_REQUIRE_MSG(version == "VERSION", "Bad " << filename << ", expected VERSION X.Y.Z header");
  EKAT_REQUIRE_MSG(version_val == p3_version, "Bad " << filename << ", expected version " << p3_version << ", but got " << version_val);

  // read tables
  double dum_s; int dum_i; // dum_s needs to be double to stream correctly
  for (int jj = 0; jj < densize; ++jj) {
    for (int ii = 0; ii < rimsize; ++ii) {
      for (int i = 0; i < isize; ++i) {
        in >> dum_i >> dum_i;
        int j_idx = 0;
        for (int j = 0; j < 15; ++j) {
          in >> dum_s;
          if (j > 1 && j != 10) {
            ice_table_vals_h(jj, ii, i, j_idx++) = dum_s;
          }
        }
      }

      for (int i = 0; i < isize; ++i) {
        for (int j = 0; j < rcollsize; ++j) {
          in >> dum_i >> dum_i;
          int k_idx = 0;
          for (int k = 0; k < 6; ++k) {
            in >> dum_s;
            if (k == 3 || k == 4) {
              collect_table_vals_h(jj, ii, i, j, k_idx++) = std::log10(dum_s);
            }
          }
        }
      }
    }
  }

  // deep copy to device
  Kokkos::deep_copy(ice_table_vals_d, ice_table_vals_h);
  Kokkos::deep_copy(collect_table_vals_d, collect_table_vals_h);
  ice_table_vals    = ice_table_vals_d;
  collect_table_vals = collect_table_vals_d;
}

template <typename S, typename C, typename MuRT, typename VNT, typename VMT, typename RevapT>
void compute_tables(const bool masterproc, MuRT& mu_r_table_vals, VNT& vn_table_vals, VMT& vm_table_vals, RevapT& revap_table_vals)
{
  using c = scream::physics::Constants<S>;

  int ii,jj,kk;
  S lamr,mu_r,dm,dum1,dum2,dum3,dum4,dum5,dd,amg,vt,dia;

  using MuRT_NC   = typename MuRT::non_const_type;
  using VNT_NC    = typename VNT::non_const_type;
  using VMT_NC    = typename VMT::non_const_type;
  using RevapT_NC = typename RevapT::non_const_type;

  MuRT_NC   mu_r_table_vals_nc("mu_r_table_vals");
  VNT_NC    vn_table_vals_nc("vn_table_vals");
  VMT_NC    vm_table_vals_nc("vm_table_vals");
  RevapT_NC revap_table_vals_nc("revap_table_vals");

  // Get host views
  auto mu_r_table_vals_h  = Kokkos::create_mirror_view(mu_r_table_vals_nc);
  auto revap_table_vals_h = Kokkos::create_mirror_view(revap_table_vals_nc);
  auto vn_table_vals_h    = Kokkos::create_mirror_view(vn_table_vals_nc);
  auto vm_table_vals_h    = Kokkos::create_mirror_view(vm_table_vals_nc);

  if (masterproc) {
    std::cout << "Recomputing lookup (non-ice) tables" << std::endl;
  }

  // ------------------------------------------------------------------------------------------

  // Generate lookup table for rain shape parameter mu_r
  // this is very fast so it can be generated at the start of each run
  // make a 150x1 1D lookup table, this is done in parameter
  // space of a scaled mean size proportional qr/Nr -- initlamr

  // write(iulog,*) '   Generating rain lookup-table ...'

  // AaronDonahue: Switching to table ver 4 means switching to a constand mu_r,
  // so this section is commented out.
  Kokkos::deep_copy(mu_r_table_vals_h, 1); // mu_r_constant =1. In other places, this is runtime_options.constant_mu_rain

  static constexpr S thrd = 1./3;
  static constexpr S small = 1.e-30;

  //.......................................................................
  // Generate lookup table for rain fallspeed and ventilation parameters
  // the lookup table is two dimensional as a function of number-weighted mean size
  // proportional to qr/Nr and shape parameter mu_r
  for (ii = 1; ii <= 10; ++ii) {
    mu_r = 1; // mu_r_constant = 1

    // loop over number-weighted mean size
    for (jj = 1; jj <= 300; ++jj) {
      if (jj <= 20) {
        dm = (jj*10 - 5)*1.e-6; // mean size [m]
      }
      else {
        dm = ((jj-20)*30 + 195)*1.e-6; // mean size [m]
      }

      lamr = (mu_r + 1)/dm;

      // do numerical integration over PSD

      dum1 = 0; // numerator,   number-weighted fallspeed
      dum2 = 0; // denominator, number-weighted fallspeed
      dum3 = 0; // numerator,   mass-weighted fallspeed
      dum4 = 0; // denominator, mass-weighted fallspeed
      dum5 = 0; // term for ventilation factor in evap
      dd   = 2;

      // loop over PSD to numerically integrate number and mass-weighted mean fallspeeds
      for (kk = 1; kk <= 10000; ++kk) {

        dia = (kk*dd - dd*0.5)*1.e-6;     // size bin [m]
        amg = c::PIOV6*997 * std::pow(dia, 3); // mass [kg]
        amg = amg*1000;                   // convert [kg] to [g]

        // get fallspeed as a function of size [m s-1]
        if (dia*1.e+6 <= 134.43) {
          vt = 4.5795e+3 * std::pow(amg, 2*thrd);
        }
        else if (dia*1.e+6 < 1511.64) {
          vt = 4.962e+1 * std::pow(amg, thrd);
        }
        else if (dia*1.e+6 < 3477.84) {
          vt = 1.732e+1 * std::pow(amg, c::SXTH);
        }
        else {
          vt = 9.17;
        }

        // note: factor of 4.*mu_r is non-answer changing and only needed to
        //       prevent underflow/overflow errors, same with 3.*mu_r for dum5
        dum1 += vt * std::pow(10, mu_r*std::log10(dia) + 4*mu_r) * std::exp(-lamr*dia) * dd * 1.e-6;
        dum2 += std::pow(10, mu_r*std::log10(dia) + 4*mu_r) * std::exp(-lamr*dia) * dd * 1.e-6;
        dum3 += vt * std::pow(10, (mu_r+3)*std::log10(dia) + 4*mu_r) * std::exp(-lamr*dia) * dd * 1.e-6;
        dum4 += std::pow(10, (mu_r+3)*std::log10(dia) + 4*mu_r) * std::exp(-lamr*dia) * dd * 1.e-6;
        dum5 += std::pow(vt*dia, 0.5) * std::pow(10, (mu_r+1)*std::log10(dia) + 3*mu_r) * std::exp(-lamr*dia) * dd * 1.e-6;
      }

      dum2 = std::max(dum2, small); // to prevent divide-by-zero below
      dum4 = std::max(dum4, small); // to prevent divide-by-zero below
      dum5 = std::max(dum5, small); // to prevent log10-of-zero below

      vn_table_vals_h(jj-1,ii-1)    = dum1/dum2;
      vm_table_vals_h(jj-1,ii-1)    = dum3/dum4;
      revap_table_vals_h(jj-1,ii-1) = std::pow(10, std::log10(dum5) + (mu_r+1)*std::log10(lamr) - (3*mu_r));
    }
  }

  Kokkos::deep_copy(mu_r_table_vals_nc, mu_r_table_vals_h);
  Kokkos::deep_copy(revap_table_vals_nc, revap_table_vals_h);
  Kokkos::deep_copy(vn_table_vals_nc, vn_table_vals_h);
  Kokkos::deep_copy(vm_table_vals_nc, vm_table_vals_h);

  mu_r_table_vals = mu_r_table_vals_nc;
  vn_table_vals = vn_table_vals_nc;
  vm_table_vals = vm_table_vals_nc;
  revap_table_vals = revap_table_vals_nc;
}

template <bool IsRead, typename S>
static void action(const ekat::FILEPtr& fid, S* data, const size_t size)
{
  if constexpr (IsRead) {
    ekat::read(data, size, fid);
  }
  else {
    ekat::write(data, size, fid);
  }
}

template <bool IsRead, typename MuRT, typename VNT, typename VMT, typename RevapT>
void io_impl(const bool masterproc, const char* dir, MuRT& mu_r_table_vals, VNT& vn_table_vals, VMT& vm_table_vals, RevapT& revap_table_vals)
{
  if (masterproc) {
    std::cout << (IsRead ? "Reading" : "Writing") << " lookup (non-ice) tables in dir " << dir << std::endl;
  }

  std::string extension =
#ifdef SCREAM_DOUBLE_PRECISION
    "8"
#else
    "4"
#endif
    ;

  const char* rw_flag = IsRead ? "r" : "w";

  // Get host views
  auto mu_r_table_vals_h  = Kokkos::create_mirror_view(mu_r_table_vals);
  auto revap_table_vals_h = Kokkos::create_mirror_view(revap_table_vals);
  auto vn_table_vals_h    = Kokkos::create_mirror_view(vn_table_vals);
  auto vm_table_vals_h    = Kokkos::create_mirror_view(vm_table_vals);

  // Add v2 because these tables are not identical to v1 due to roundoff differences
  // caused by doing the math in C++ instead of f90.
  std::string mu_r_filename  = std::string(dir) + "/mu_r_table_vals_v2.dat" + extension;
  std::string revap_filename = std::string(dir) + "/revap_table_vals_v2.dat" + extension;
  std::string vn_filename    = std::string(dir) + "/vn_table_vals_v2.dat" + extension;
  std::string vm_filename    = std::string(dir) + "/vm_table_vals_v2.dat" + extension;

  ekat::FILEPtr mu_r_file(fopen(mu_r_filename.c_str(), rw_flag));
  ekat::FILEPtr revap_file(fopen(revap_filename.c_str(), rw_flag));
  ekat::FILEPtr vn_file(fopen(vn_filename.c_str(), rw_flag));
  ekat::FILEPtr vm_file(fopen(vm_filename.c_str(), rw_flag));

  // Read files
  action<IsRead>(mu_r_file, mu_r_table_vals_h.data(), mu_r_table_vals.size());
  action<IsRead>(revap_file, revap_table_vals_h.data(), revap_table_vals.size());
  action<IsRead>(vn_file, vn_table_vals_h.data(), vn_table_vals.size());
  action<IsRead>(vm_file, vm_table_vals_h.data(), vm_table_vals.size());

  // Copy back to device
  if constexpr (IsRead) {
    Kokkos::deep_copy(mu_r_table_vals, mu_r_table_vals_h);
    Kokkos::deep_copy(revap_table_vals, revap_table_vals_h);
    Kokkos::deep_copy(vn_table_vals, vn_table_vals_h);
    Kokkos::deep_copy(vm_table_vals, vm_table_vals_h);
  }
}

template <typename MuRT, typename VNT, typename VMT, typename RevapT>
void read_computed_tables(const bool masterproc, const char* dir, MuRT& mu_r_table_vals, VNT& vn_table_vals, VMT& vm_table_vals, RevapT& revap_table_vals)
{
  using MuRT_NC   = typename MuRT::non_const_type;
  using VNT_NC    = typename VNT::non_const_type;
  using VMT_NC    = typename VMT::non_const_type;
  using RevapT_NC = typename RevapT::non_const_type;

  MuRT_NC   mu_r_table_vals_nc("mu_r_table_vals");
  VNT_NC    vn_table_vals_nc("vn_table_vals");
  VMT_NC    vm_table_vals_nc("vm_table_vals");
  RevapT_NC revap_table_vals_nc("revap_table_vals");

  io_impl<true>(masterproc, dir, mu_r_table_vals_nc, vn_table_vals_nc, vm_table_vals_nc, revap_table_vals_nc);

  mu_r_table_vals = mu_r_table_vals_nc;
  vn_table_vals = vn_table_vals_nc;
  vm_table_vals = vm_table_vals_nc;
  revap_table_vals = revap_table_vals_nc;
}

template <typename MuRT, typename VNT, typename VMT, typename RevapT>
void write_computed_tables(const bool masterproc, const char* dir, const MuRT& mu_r_table_vals, const VNT& vn_table_vals, const VMT& vm_table_vals, const RevapT& revap_table_vals)
{
  io_impl<false>(masterproc, dir, mu_r_table_vals, vn_table_vals, vm_table_vals, revap_table_vals);
}

template <typename S, typename DnuT>
void compute_dnu(DnuT& dnu_table_vals)
{
  typename DnuT::non_const_type dnu_table_vals_non_const("dnu_table_vals");
  const auto dnu_table_h   = Kokkos::create_mirror_view(dnu_table_vals_non_const);
  dnu_table_h(0)  =  0.000;
  dnu_table_h(1)  = -0.557;
  dnu_table_h(2)  = -0.430;
  dnu_table_h(3)  = -0.307;
  dnu_table_h(4)  = -0.186;
  dnu_table_h(5)  = -0.067;
  dnu_table_h(6)  = -0.050;
  dnu_table_h(7)  = -0.167;
  dnu_table_h(8)  = -0.282;
  dnu_table_h(9)  = -0.397;
  dnu_table_h(10) = -0.512;
  dnu_table_h(11) = -0.626;
  dnu_table_h(12) = -0.739;
  dnu_table_h(13) = -0.853;
  dnu_table_h(14) = -0.966;
  dnu_table_h(15) = -0.966;
  Kokkos::deep_copy(dnu_table_vals_non_const, dnu_table_h);
  dnu_table_vals = DnuT(dnu_table_vals_non_const);
}

}

/*
 * Implementation of p3 init. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */
template <typename S, typename D>
typename Functions<S,D>::P3LookupTables Functions<S,D>
::p3_init (const bool write_tables, const bool masterproc) {
  P3LookupTables lookup_tables; // This struct could be our global singleton
  auto version = P3C::p3_version;
  auto p3_lookup_base = P3C::p3_lookup_base;
  static const char* dir = SCREAM_DATA_DIR "/tables";
  // p3_init_a (reads ice_table, collect_table)
  read_ice_lookup_tables<S>(masterproc, p3_lookup_base, version, lookup_tables.ice_table_vals, lookup_tables.collect_table_vals, P3C::densize, P3C::rimsize, P3C::isize, P3C::rcollsize);
  if (write_tables) {
    //p3_init_b (computes tables mu_r_table, revap_table, vn_table, vm_table)
    compute_tables<S, P3C>(masterproc, lookup_tables.mu_r_table_vals, lookup_tables.vn_table_vals, lookup_tables.vm_table_vals, lookup_tables.revap_table_vals);
    write_computed_tables(masterproc, dir, lookup_tables.mu_r_table_vals, lookup_tables.vn_table_vals, lookup_tables.vm_table_vals, lookup_tables.revap_table_vals);
  }
  else {
    read_computed_tables(masterproc, dir, lookup_tables.mu_r_table_vals, lookup_tables.vn_table_vals, lookup_tables.vm_table_vals, lookup_tables.revap_table_vals);
  }
  // dnu is always computed/hardcoded
  compute_dnu<S>(lookup_tables.dnu_table_vals);

  return lookup_tables;
}

} // namespace p3
} // namespace scream

#endif
