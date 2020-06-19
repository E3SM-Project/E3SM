#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestBackToCellAverage {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  // Generate n test structs, each populated with random data (values within
  // [0,1]) by the default constructor.
  BackToCellAverageData back_to_cell_average_data[Spack::n];
  for (Int i = 0; i < Spack::n; ++i) {
    back_to_cell_average_data[i].randomize();
  }

  // Sync to device.
  view_1d<BackToCellAverageData> device_data("back_to_cell_average", Spack::n);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&back_to_cell_average_data[0], &back_to_cell_average_data[0] + Spack::n,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < Spack::n; ++i) {
    back_to_cell_average(back_to_cell_average_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack lcldm, rcldm, icldm, qcacc, qrevp, qcaut, ncacc, ncslf, ncautc, nrslf,
      nrevp, ncautr, qisub, nrshdr, qcheti, qrcol, qcshd, qimlt,
      qccol, qrheti, nimlt, nccol, ncshdc, ncheti, nrcol, nislf, qidep, nrheti,
      nisub, qinuc, ninuc, qiberg;
    for (Int s = 0; s < Spack::n; ++s) {
      lcldm[s] = device_data[s].lcldm;
      rcldm[s] = device_data[s].rcldm;
      icldm[s] = device_data[s].icldm;
      qcacc[s] = device_data[s].qcacc;
      qrevp[s] = device_data[s].qrevp;
      qcaut[s] = device_data[s].qcaut;
      ncacc[s] = device_data[s].ncacc;
      ncslf[s] = device_data[s].ncslf;
      ncautc[s] = device_data[s].ncautc;
      nrslf[s] = device_data[s].nrslf;
      nrevp[s] = device_data[s].nrevp;
      ncautr[s] = device_data[s].ncautr;
      qisub[s] = device_data[s].qisub;
      nrshdr[s] = device_data[s].nrshdr;
      qcheti[s] = device_data[s].qcheti;
      qrcol[s] = device_data[s].qrcol;
      qcshd[s] = device_data[s].qcshd;
      qimlt[s] = device_data[s].qimlt;
      qccol[s] = device_data[s].qccol;
      qrheti[s] = device_data[s].qrheti;
      nimlt[s] = device_data[s].nimlt;
      nccol[s] = device_data[s].nccol;
      ncshdc[s] = device_data[s].ncshdc;
      ncheti[s] = device_data[s].ncheti;
      nrcol[s] = device_data[s].nrcol;
      nislf[s] = device_data[s].nislf;
      qidep[s] = device_data[s].qidep;
      nrheti[s] = device_data[s].nrheti;
      nisub[s] = device_data[s].nisub;
      qinuc[s] = device_data[s].qinuc;
      ninuc[s] = device_data[s].ninuc;
      qiberg[s] = device_data[s].qiberg;
    }

    Functions::back_to_cell_average(lcldm, rcldm, icldm, qcacc, qrevp, qcaut,
      ncacc, ncslf, ncautc, nrslf, nrevp, ncautr, qisub, nrshdr,
      qcheti, qrcol, qcshd, qimlt, qccol, qrheti, nimlt, nccol, ncshdc, ncheti,
      nrcol, nislf, qidep, nrheti, nisub, qinuc, ninuc, qiberg);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      device_data(s).qcacc = qcacc[s];
      device_data(s).qrevp = qrevp[s];
      device_data(s).qcaut = qcaut[s];
      device_data(s).ncacc = ncacc[s];
      device_data(s).ncslf = ncslf[s];
      device_data(s).ncautc = ncautc[s];
      device_data(s).nrslf = nrslf[s];
      device_data(s).nrevp = nrevp[s];
      device_data(s).ncautr = ncautr[s];
      device_data(s).qisub = qisub[s];
      device_data(s).nrshdr = nrshdr[s];
      device_data(s).qcheti = qcheti[s];
      device_data(s).qrcol = qrcol[s];
      device_data(s).qcshd = qcshd[s];
      device_data(s).qimlt = qimlt[s];
      device_data(s).qccol = qccol[s];
      device_data(s).qrheti = qrheti[s];
      device_data(s).nimlt = nimlt[s];
      device_data(s).nccol = nccol[s];
      device_data(s).ncshdc = ncshdc[s];
      device_data(s).ncheti = ncheti[s];
      device_data(s).nrcol = nrcol[s];
      device_data(s).nislf = nislf[s];
      device_data(s).qidep = qidep[s];
      device_data(s).nrheti = nrheti[s];
      device_data(s).nisub = nisub[s];
      device_data(s).qinuc = qinuc[s];
      device_data(s).ninuc = ninuc[s];
      device_data(s).qiberg = qiberg[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < Spack::n; ++s) {
    REQUIRE(back_to_cell_average_data[s].qcacc == host_data[s].qcacc);
    REQUIRE(back_to_cell_average_data[s].qrevp == host_data[s].qrevp);
    REQUIRE(back_to_cell_average_data[s].qcaut == host_data[s].qcaut);
    REQUIRE(back_to_cell_average_data[s].ncacc == host_data[s].ncacc);
    REQUIRE(back_to_cell_average_data[s].ncslf == host_data[s].ncslf);
    REQUIRE(back_to_cell_average_data[s].ncautc == host_data[s].ncautc);
    REQUIRE(back_to_cell_average_data[s].nrslf == host_data[s].nrslf);
    REQUIRE(back_to_cell_average_data[s].nrevp == host_data[s].nrevp);
    REQUIRE(back_to_cell_average_data[s].ncautr == host_data[s].ncautr);
    REQUIRE(back_to_cell_average_data[s].qisub == host_data[s].qisub);
    REQUIRE(back_to_cell_average_data[s].nrshdr == host_data[s].nrshdr);
    REQUIRE(back_to_cell_average_data[s].qcheti == host_data[s].qcheti);
    REQUIRE(back_to_cell_average_data[s].qrcol == host_data[s].qrcol);
    REQUIRE(back_to_cell_average_data[s].qcshd == host_data[s].qcshd);
    REQUIRE(back_to_cell_average_data[s].qimlt == host_data[s].qimlt);
    REQUIRE(back_to_cell_average_data[s].qccol == host_data[s].qccol);
    REQUIRE(back_to_cell_average_data[s].qrheti == host_data[s].qrheti);
    REQUIRE(back_to_cell_average_data[s].nimlt == host_data[s].nimlt);
    REQUIRE(back_to_cell_average_data[s].nccol == host_data[s].nccol);
    REQUIRE(back_to_cell_average_data[s].ncshdc == host_data[s].ncshdc);
    REQUIRE(back_to_cell_average_data[s].ncheti == host_data[s].ncheti);
    REQUIRE(back_to_cell_average_data[s].nrcol == host_data[s].nrcol);
    REQUIRE(back_to_cell_average_data[s].nislf == host_data[s].nislf);
    REQUIRE(back_to_cell_average_data[s].qidep == host_data[s].qidep);
    REQUIRE(back_to_cell_average_data[s].nrheti == host_data[s].nrheti);
    REQUIRE(back_to_cell_average_data[s].nisub == host_data[s].nisub);
    REQUIRE(back_to_cell_average_data[s].qinuc == host_data[s].qinuc);
    REQUIRE(back_to_cell_average_data[s].ninuc == host_data[s].ninuc);
    REQUIRE(back_to_cell_average_data[s].qiberg == host_data[s].qiberg);
  }

}

};

}
}
}

namespace {

TEST_CASE("p3_back_to_cell_average", "[p3_functions]")
{
  using TRIF = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestBackToCellAverage;

  TRIF::run_phys();
  TRIF::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
