#include <catch2/catch.hpp>
#include <numeric>

#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "share/grid/point_grid.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("utils") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using kt = KokkosTypes<DefaultDevice>;
  using kt_host = KokkosTypes<HostDevice>;

  using P8 = ekat::Pack<Real,8>;

  ekat::Comm comm(MPI_COMM_WORLD);
  bool parallel_test = comm.size()>1;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {3,24};

  FieldIdentifier fid ("field_1", {tags,dims}, m/s,"some_grid");
  Field f1(fid);
  f1.get_header().get_alloc_properties().request_allocation(P8::n);
  f1.allocate_view();

  SECTION ("compare") {

    Field f2(fid);
    f2.allocate_view();

    auto v1 = f1.get_view<P8**>();
    auto v2 = f2.get_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto am_i_root = comm.am_i_root();
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v2(i,j) = i*dim1+j;

      int jpack = j / P8::n;
      int jvec = j % P8::n;
      v1(i,jpack)[jvec] = i*dim1+j;
      if (parallel_test and not am_i_root) {
        v1(i,jpack)[jvec] *= -1;
      }
    });
    Kokkos::fence();

    if (not parallel_test) {
      // The views were filled the same way, so they should test equal
      REQUIRE(views_are_equal(f1,f2));
    } else {
      // On rank 0, all is equal, but on other ranks all is different
      if (comm.am_i_root()) {
        REQUIRE(views_are_equal(f1,f2));
      } else {
        REQUIRE(not views_are_equal(f1,f2));
      }

      // On all ranks, the *global* comparison should fail
      REQUIRE(not views_are_equal(f1,f2,&comm));
    }
  }

  SECTION ("sum") {

    auto v1 = f1.get_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx + 1;
    });
    Kokkos::fence();

    Real lsum = offset*lsize + lsize*(lsize+1) / 2.0;
    Real gsum = gsize*(gsize+1) / 2.0;

    REQUIRE(field_sum<Real>(f1)==lsum);
    REQUIRE(field_sum<Real>(f1,&comm)==gsum);
  }

  SECTION ("frobenius") {

    auto v1 = f1.get_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx + 1;
    });
    Kokkos::fence();

    // Summing the squares from n=a+1 to n=a+N gives
    // (a+1)^2+(a+2)^2+...
    // which ultimately gives
    // N*a^2 + 2a*(1+2+...+N)
    Real lsum = offset*offset*lsize + 2*offset*lsize*(lsize+1)/2 + lsize*(lsize+1)*(2*lsize+1) / 6.0;
    Real gsum = gsize*(gsize+1)*(2*gsize+1) / 6.0;

    REQUIRE(frobenius_norm<Real>(f1)==std::sqrt(lsum));
    REQUIRE(frobenius_norm<Real>(f1,&comm)==std::sqrt(gsum));
  }

  SECTION ("max") {

    auto v1 = f1.get_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx+1;
    });
    Kokkos::fence();

    Real lmax = offset + lsize;
    Real gmax = gsize;

    REQUIRE(field_max<Real>(f1)==lmax);
    REQUIRE(field_max<Real>(f1,&comm)==gmax);
  }

  SECTION ("min") {

    auto v1 = f1.get_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx+1;
    });
    Kokkos::fence();

    Real lmin = offset + 1;
    Real gmin = 1;

    REQUIRE(field_min<Real>(f1)==lmin);
    REQUIRE(field_min<Real>(f1,&comm)==gmin);
  }

  SECTION ("perturb") {
    using namespace ShortFieldTagsNames;
    using RPDF = std::uniform_real_distribution<Real>;
    using IPDF = std::uniform_int_distribution<int>;
    auto engine = setup_random_test ();

    const int ncols = 6;
    const int ncmps = 2;
    const int nlevs = IPDF(3,9)(engine); // between 3-9 levels

    // Create 1d, 2d, 3d fields with a level dimension, and set all to 1
    FieldIdentifier fid1 ("f_1d",   FieldLayout({LEV},           {nlevs}),               Units::nondimensional(), "");
    FieldIdentifier fid2a("f_2d_a", FieldLayout({CMP, LEV},      {ncmps, nlevs}),        Units::nondimensional(), "");
    FieldIdentifier fid2b("f_2d_b", FieldLayout({COL, LEV},      {ncols, nlevs}),        Units::nondimensional(), "");
    FieldIdentifier fid3 ("f_3d",   FieldLayout({COL, CMP, LEV}, {ncols, ncmps, nlevs}), Units::nondimensional(), "");
    Field f1(fid1), f2a(fid2a), f2b(fid2b), f3(fid3);
    f1.allocate_view(), f2a.allocate_view(), f2b.allocate_view(), f3.allocate_view();
    f1.deep_copy(1), f2a.deep_copy(1), f2b.deep_copy(1), f3.deep_copy(1);

    // We need GIDs for fields with COL component. This test is not over
    // multiple ranks, so just set as [0, ncols-1].
    Field gids(FieldIdentifier("gids", FieldLayout({COL}, {ncols}), Units::nondimensional(), "", DataType::IntType));
    gids.allocate_view();
    auto gids_data = gids.get_internal_view_data<int,Host>();
    std::iota(gids_data, gids_data+ncols, 0);
    gids.sync_to_dev();

    // Create masks s.t. only last 3 levels are perturbed. For variety,
    // 1d and 2d fields will use lambda mask and 3 field will use a view.
    auto mask_lambda = [&nlevs] (const int& i0) {
      return i0 >= nlevs-3;
    };
    kt_host::view_1d<bool> mask_view("mask_view", nlevs);
    Kokkos::deep_copy(mask_view, false);
    for (int ilev=0; ilev<nlevs; ++ilev) {
      if (ilev >= nlevs-3) mask_view(ilev) = true;
    }

    // Compute random perturbation between [2, 3]
    RPDF pdf(2, 3);
    int base_seed = 0;
    perturb(f1,  engine, pdf, base_seed, mask_lambda);
    perturb(f2a, engine, pdf, base_seed, mask_lambda);
    perturb(f2b, engine, pdf, base_seed, mask_lambda, gids);
    perturb(f3,  engine, pdf, base_seed, mask_view,   gids);

    // Sync to host for checks
    f1.sync_to_host(), f2a.sync_to_host(), f2b.sync_to_host(), f3.sync_to_host();
    const auto v1  = f1.get_view <Real*,   Host>();
    const auto v2a = f2a.get_view<Real**,  Host>();
    const auto v2b = f2b.get_view<Real**,  Host>();
    const auto v3  = f3.get_view <Real***, Host>();

    // Check that all field values are 1 for all but last 3 levels and between [2,3] otherwise.
    auto check_level = [&] (const int ilev, const Real val) {
      if (ilev < nlevs-3) REQUIRE(val == 1);
      else REQUIRE((2 <= val && val <= 3));
    };
    for (int icol=0; icol<ncols; ++icol) {
      for (int icmp=0; icmp<ncmps; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          if (icol==0 && icmp==0) check_level(ilev, v1(ilev));
          if (icol==0) check_level(ilev, v2a(icmp,ilev));
          if (icmp==0) check_level(ilev, v2b(icol,ilev));
          check_level(ilev, v3(icol,icmp,ilev));
    }}}

    // Check that using a different seed gives different values
    auto f1_alt = f1.clone(); f1_alt.deep_copy(1.0);
    auto f3_alt = f3.clone(); f3_alt.deep_copy(1.0);
    int base_seed_alt = 100;
    perturb(f1_alt, engine, pdf, base_seed_alt, mask_lambda);
    perturb(f3_alt, engine, pdf, base_seed_alt, mask_lambda, gids);
    f1_alt.sync_to_host(), f3_alt.sync_to_host();

    const auto v1_alt = f1_alt.get_view<Real*,   Host>();
    const auto v3_alt = f3_alt.get_view<Real***, Host>();

    auto check_diff = [&] (const int ilev, const Real val1, const Real val2) {
      if (ilev < nlevs-3) REQUIRE(val1==val2);
      else                REQUIRE(val1!=val2);
    };
    for (int icol=0; icol<ncols; ++icol) {
      for (int icmp=0; icmp<ncmps; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          if (icol==0 && icmp==0) check_diff(ilev, v1(ilev), v1_alt(ilev));
          check_diff(ilev, v3(icol,icmp,ilev), v3_alt(icol,icmp,ilev));
    }}}

    // Finally check that the original seed gives same result
    f1_alt.deep_copy(1.0), f3_alt.deep_copy(1.0);
    perturb(f1_alt, engine, pdf, base_seed, mask_lambda);
    perturb(f3_alt, engine, pdf, base_seed, mask_lambda, gids);
    REQUIRE(views_are_equal(f1, f1_alt));
    REQUIRE(views_are_equal(f3, f3_alt));
  }

  SECTION ("wrong_st") {
    using wrong_real =
      typename std::conditional<std::is_same<Real,double>::value,
                                float,double>::type;
    REQUIRE_THROWS(field_min<int>(f1));
    REQUIRE_THROWS(field_max<int>(f1));
    REQUIRE_THROWS(field_sum<wrong_real>(f1));
    REQUIRE_THROWS(frobenius_norm<wrong_real>(f1));
  }
}

TEST_CASE ("print_field_hyperslab") {
  using namespace scream;
  using namespace ekat::units;

  using namespace ShortFieldTagsNames;
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;

  // Setup random number generation
  ekat::Comm comm(MPI_COMM_WORLD);
  auto engine = setup_random_test ();
  RPDF pdf(0,1);

  const int nel = 2;
  const int ncmp = 2;
  const int ngp  = 4;
  const int nlev = 12;
  const int iel  = IPDF(0,nel-1)(engine);
  const int icmp = IPDF(0,ncmp-1)(engine);
  const int igp  = IPDF(0,ngp-1)(engine);
  const int jgp  = IPDF(0,ngp-1)(engine);
  const int ilev = IPDF(0,nlev-1)(engine);

  // WARNING: make sure this value matches the one used in print_field_hyperslab
  constexpr int max_per_line = 5;

  // Create field (if available, use packs, to ensure we don't print garbage)
  std::vector<FieldTag> tags = {EL, CMP, GP, GP, LEV};
  std::vector<int>      dims = {nel,ncmp,ngp,ngp,nlev};
  int pack_size = std::min(SCREAM_PACK_SIZE,4);

  FieldIdentifier fid ("f", {tags,dims}, kg, "some_grid");
  Field f (fid);
  f.get_header().get_alloc_properties().request_allocation(pack_size);
  f.allocate_view();
  randomize (f,engine,pdf);

  auto v = f.get_view<const Real*****,Host>();

  SECTION ("slice_0") {
    std::vector<FieldTag> loc_tags = {EL,CMP};
    std::vector<int>      loc_idxs = {iel,icmp};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << to_string(fid.get_layout()) << "\n\n";
      for (int gp1=0; gp1<ngp; ++gp1) {
        for (int gp2=0; gp2<ngp; ++gp2) {
          expected <<  "  f(" << iel << "," << icmp << "," << gp1 << "," << gp2 << ",:)";
          for (int lev=0; lev<nlev; ++lev) {
            if (lev % max_per_line == 0) {
              expected << "\n    ";
            }
            expected << v(iel,icmp,gp1,gp2,lev) << ", ";
          }
          expected << "\n";
        }
      }

    REQUIRE (out.str()==expected.str());
  }
  SECTION ("slice_0234") {
    std::vector<FieldTag> loc_tags = {EL,GP,GP,LEV};
    std::vector<int>      loc_idxs = {iel,igp,jgp,ilev};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << to_string(fid.get_layout()) << "\n\n";
      expected <<  "  f(" << iel << ",:," << igp << "," << jgp << "," << ilev << ")";
      for (int cmp=0; cmp<ncmp; ++cmp) {
        if (cmp % max_per_line == 0) {
          expected << "\n    ";
        }
        expected << v(iel,cmp,igp,jgp,ilev) << ", ";
      }
      expected << "\n";

    REQUIRE (out.str()==expected.str());
  }
  SECTION ("slice_01234") {
    std::vector<FieldTag> loc_tags = {EL,CMP,GP,GP,LEV};
    std::vector<int>      loc_idxs = {iel,icmp,igp,jgp,ilev};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << to_string(fid.get_layout()) << "\n\n";
      expected <<  "  f(" << ekat::join(loc_idxs,",") << ")\n";
      expected << "    " << v(iel,icmp,igp,jgp,ilev) << ", \n";

    REQUIRE (out.str()==expected.str());
  }
}

} // anonymous namespace
