#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"

#include <ekat_pack.hpp>
#include <ekat_subview_utils.hpp>

namespace {

TEST_CASE("field_utils") {
  using namespace scream;
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using kt = KokkosTypes<DefaultDevice>;

  using P8 = ekat::Pack<Real,8>;

  ekat::Comm comm(MPI_COMM_WORLD);
  bool parallel_test = comm.size()>1;

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {3,24};

  FieldIdentifier fid ("field_1", {tags,dims}, m/s,"some_grid");
  Field f1(fid);
  f1.get_header().get_alloc_properties().request_allocation(P8::n);
  f1.allocate_view();

  SECTION("compare-rank-0") {
    // create two fields with rank-0 and get their views
    std::vector<FieldTag> tags_0 = {};
    std::vector<int> dims_0      = {};
    FieldIdentifier fid_01("field_01", {tags_0, dims_0}, m / s, "some_grid");
    FieldIdentifier fid_02("field_02", {tags_0, dims_0}, m / s, "some_grid");
    Field f01(fid_01);
    Field f02(fid_02);
    f01.allocate_view();
    f02.allocate_view();
    auto f01v = f01.get_view<Real>();
    auto f02v = f02.get_view<Real>();
    // fill the views with the same values
    Real val = 54321;
    Kokkos::deep_copy(f01v, val);
    Kokkos::deep_copy(f02v, val);
    // check that the views are equal
    REQUIRE(views_are_equal(f01, f02));

    // fill the views with different values
    Kokkos::deep_copy(f02v, 1 / val);
    // check that the views are not equal
    REQUIRE(not views_are_equal(f01, f02));
  }

  SECTION ("compare") {

    Field f2(fid);
    f2.allocate_view();

    auto v1 = f1.get_strided_view<P8**>();
    auto v2 = f2.get_strided_view<Real**>();
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

  SECTION ("perturb") {
    using namespace ShortFieldTagsNames;
    using IPDF = std::uniform_int_distribution<int>;

    auto nondim = Units::nondimensional();
    int seed = get_random_test_seed();
    std::mt19937_64 engine(seed);

    const int ncols = 6;
    const int ncmps = 2;
    const int nlevs = IPDF(3,9)(engine); // between 3-9 levels

    // Create 1d, 2d, 3d fields with a level dimension, and set all to 1
    FieldIdentifier fid1 ("f_1d",   FieldLayout({LEV},           {nlevs}),               nondim, "");
    FieldIdentifier fid2a("f_2d_a", FieldLayout({CMP, LEV},      {ncmps, nlevs}),        nondim, "");
    FieldIdentifier fid2b("f_2d_b", FieldLayout({COL, LEV},      {ncols, nlevs}),        nondim, "");
    FieldIdentifier fid3 ("f_3d",   FieldLayout({COL, CMP, LEV}, {ncols, ncmps, nlevs}), nondim, "");
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
    FieldIdentifier pmask_fid("pmask",fid1.get_layout(),nondim,"",DataType::IntType);
    Field pmask(pmask_fid,true);
    auto pmask_h = pmask.get_view<int*,Host>();
    for (int i=0; i<nlevs; ++i) {
      pmask_h(i) = static_cast<int>(i>=nlevs-3);
    }
    pmask.sync_to_dev();

    // Compute random perturbation of magnitude 0.1
    int base_seed = 0;
    perturb(f1,  0.1, base_seed, pmask);
    perturb(f2a, 0.1, base_seed, pmask);
    perturb(f2b, 0.1, base_seed, pmask, gids);
    perturb(f3,  0.1, base_seed, pmask, gids);

    // Sync to host for checks
    f1.sync_to_host(), f2a.sync_to_host(), f2b.sync_to_host(), f3.sync_to_host();
    const auto v1  = f1.get_strided_view <Real*,   Host>();
    const auto v2a = f2a.get_strided_view<Real**,  Host>();
    const auto v2b = f2b.get_strided_view<Real**,  Host>();
    const auto v3  = f3.get_strided_view <Real***, Host>();

    // Check that all field values are 1 for all but last 3 levels and between [0.9,1.1] otherwise.
    int num_perturbed = 0;
    auto check_level = [&] (const int ilev, const Real val) {
      if (ilev < nlevs-3) {
        REQUIRE(val == 1);
      } else {
        REQUIRE((0.9 <= val && val <= 1.1));
        if (val!=1)
          ++num_perturbed;
      }
    };
    for (int icol=0; icol<ncols; ++icol) {
      for (int icmp=0; icmp<ncmps; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          if (icol==0 && icmp==0) check_level(ilev, v1(ilev));
          if (icol==0) check_level(ilev, v2a(icmp,ilev));
          if (icmp==0) check_level(ilev, v2b(icol,ilev));
          check_level(ilev, v3(icol,icmp,ilev));
    }}}

    // While it is technically possible that a uniform pdf(-0.1,0.1) returns 0,
    // we DEFINITELY don't expect it to happen more than once. So num_perturb should be AT LEAST >0
    REQUIRE (num_perturbed>0);

    // Check that using a different seed gives different values
    auto f1_alt = f1.clone(); f1_alt.deep_copy(1.0);
    auto f3_alt = f3.clone(); f3_alt.deep_copy(1.0);
    int base_seed_alt = 100;
    perturb(f1_alt, 0.1, base_seed_alt, pmask);
    perturb(f3_alt, 0.1, base_seed_alt, pmask, gids);
    f1_alt.sync_to_host(), f3_alt.sync_to_host();

    const auto v1_alt = f1_alt.get_strided_view<Real*,   Host>();
    const auto v3_alt = f3_alt.get_strided_view<Real***, Host>();

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
    perturb(f1_alt, 0.1, base_seed, pmask);
    perturb(f3_alt, 0.1, base_seed, pmask, gids);
    REQUIRE(views_are_equal(f1, f1_alt));
    REQUIRE(views_are_equal(f3, f3_alt));
  }
}

TEST_CASE ("print_field_hyperslab") {
  using namespace scream;
  using namespace ekat::units;

  using namespace ShortFieldTagsNames;
  using IPDF = std::uniform_int_distribution<int>;

  // Setup random number generation
  ekat::Comm comm(MPI_COMM_WORLD);
  int seed = get_random_test_seed();

  std::mt19937_64 engine(seed);

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
  randomize_uniform (f, seed++);

  auto v = f.get_strided_view<const Real*****,Host>();

  SECTION ("slice_0") {
    std::vector<FieldTag> loc_tags = {EL,CMP};
    std::vector<int>      loc_idxs = {iel,icmp};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << fid.get_layout().to_string() << "\n\n";
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

    loc_idxs[0] = nel;
    REQUIRE_THROWS(print_field_hyperslab(f,loc_tags,loc_idxs,out));
    loc_idxs[0] = -1;
    REQUIRE_THROWS(print_field_hyperslab(f,loc_tags,loc_idxs,out));
  }
  SECTION ("slice_0234") {
    std::vector<FieldTag> loc_tags = {EL,GP,GP,LEV};
    std::vector<int>      loc_idxs = {iel,igp,jgp,ilev};
    std::stringstream out;
    print_field_hyperslab(f,loc_tags,loc_idxs,out);

    std::stringstream expected;
    expected << "     f" << fid.get_layout().to_string() << "\n\n";
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
    expected << "     f" << fid.get_layout().to_string() << "\n\n";
      expected <<  "  f(" << ekat::join(loc_idxs,",") << ")\n";
      expected << "    " << v(iel,icmp,igp,jgp,ilev) << ", \n";

    REQUIRE (out.str()==expected.str());
  }
}

TEST_CASE ("transpose") {
  using namespace scream;

  using namespace ShortFieldTagsNames;

  // Setup random number generation
  ekat::Comm comm(MPI_COMM_WORLD);
  int seed = get_random_test_seed();

  const int ncols = 3;
  const int nlevs = 10;
  const int ncmp1 = 2;
  const int ncmp2 = 4;
  const auto u1 = ekat::units::s;
  const auto u2 = ekat::units::m;

  std::vector<FieldTag> tags1d = {COL};
  std::vector<FieldTag> tags2d = {COL, LEV};
  std::vector<FieldTag> tags3d = {COL, CMP, LEV};
  std::vector<FieldTag> tags4d = {COL, CMP, CMP, LEV};

  std::vector<int>      dims1d = {ncols};
  std::vector<int>      dims2d = {ncols,nlevs};
  std::vector<int>      dims3d = {ncols,ncmp1,nlevs};
  std::vector<int>      dims4d = {ncols,ncmp1,ncmp2,nlevs};

  FieldIdentifier fid1d  ("foo", {tags1d,dims1d}, u1, "some_grid");
  FieldIdentifier fid2d  ("foo", {tags2d,dims2d}, u1, "some_grid");
  FieldIdentifier fid3d  ("foo", {tags3d,dims3d}, u1, "some_grid");
  FieldIdentifier fid3di ("foo", {tags3d,dims3d}, u1, "some_grid", DataType::IntType);
  FieldIdentifier fid4d  ("foo", {tags4d,dims4d}, u1, "some_grid");
  FieldIdentifier fid2du2("foo", {tags4d,dims4d}, u2, "some_grid");

  SECTION ("exceptions") {
    Field f2d  (fid2d);
    Field f2du2(fid2du2);
    Field f3d  (fid3d);
    Field f3di (fid3di);

    REQUIRE_THROWS(transpose(f2d)); // not allocated

    f2d.allocate_view();
    f2du2.allocate_view();
    f3d.allocate_view();
    f3di.allocate_view();

    REQUIRE_THROWS(transpose(f2d,f2du2));; // different units
    REQUIRE_THROWS(transpose(f2d,f3d));;   // different layout
    REQUIRE_THROWS(transpose(f3di,f3d));;  // different data type
  }

  SECTION ("1d") {
    Field f1d (fid1d);
    f1d.allocate_view();
    randomize_uniform(f1d, seed++);

    auto f1d_t = transpose(f1d);
    REQUIRE(views_are_equal(f1d,f1d_t));
  }

  SECTION ("2d") {
    Field f2d (fid2d);
    f2d.allocate_view();
    randomize_uniform(f2d, seed++);
    auto f2d_t = transpose(f2d);

    f2d.sync_to_host();
    f2d_t.sync_to_host();

    auto f2d_h = f2d.get_view<const Real**,Host>();
    auto f2d_t_h = f2d_t.get_view<const Real**,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        REQUIRE (f2d_h(icol,ilev)==f2d_t_h(ilev,icol));
      }
    }
  }

  SECTION ("3d") {
    Field f3d (fid3d);
    f3d.allocate_view();
    randomize_uniform(f3d, seed++);
    auto f3d_t = transpose(f3d);

    f3d.sync_to_host();
    f3d_t.sync_to_host();

    auto f3d_h = f3d.get_view<const Real***,Host>();
    auto f3d_t_h = f3d_t.get_view<const Real***,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      for (int icmp=0; icmp<ncmp1; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          REQUIRE (f3d_h(icol,icmp,ilev)==f3d_t_h(ilev,icmp,icol));
        }
      }
    }
  }

  SECTION ("4d") {
    Field f4d (fid4d);
    f4d.allocate_view();
    randomize_uniform(f4d, seed++);
    auto f4d_t = transpose(f4d);

    f4d.sync_to_host();
    f4d_t.sync_to_host();

    auto f4d_h = f4d.get_view<const Real****,Host>();
    auto f4d_t_h = f4d_t.get_view<const Real****,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      for (int icmp=0; icmp<ncmp1; ++icmp) {
        for (int jcmp=0; jcmp<ncmp2; ++jcmp) {
          for (int ilev=0; ilev<nlevs; ++ilev) {
            REQUIRE (f4d_h(icol,icmp,jcmp,ilev)==f4d_t_h(ilev,jcmp,icmp,icol));
          }
        }
      }
    }
  }
}

} // anonymous namespace
