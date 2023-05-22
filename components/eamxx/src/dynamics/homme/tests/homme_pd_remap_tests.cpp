#include <catch2/catch.hpp>

#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "share/field/field.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "mpi/BoundaryExchange.hpp"
#include "SimulationParams.hpp"

#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_grids_manager.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <memory>
#include <random>
#include <numeric>

extern "C" {
// These are specific C/F calls for these tests (i.e., not part of scream_homme_interface.hpp)
void init_test_params_f90 ();
void cleanup_test_f90 ();
}

namespace {

TEST_CASE("remap", "") {

  using namespace scream;
  using namespace ShortFieldTagsNames;

  // Some type defs
  using Remapper = PhysicsDynamicsRemapper;
  using IPDF = std::uniform_int_distribution<int>;
  using FID = FieldIdentifier;
  using FL  = FieldLayout;

  constexpr int pg_gll = 0;
  constexpr int PackSize = HOMMEXX_VECTOR_SIZE;

  // Create a comm
  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test(&comm);

  // Init homme context
  if (!is_parallel_inited_f90()) {
    auto comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
    init_parallel_f90(comm_f);
  }
  init_test_params_f90 ();

  // We'll use this extensively, so let's use a short ref name
  auto& c = Homme::Context::singleton();

  // Set a value for qsize that is not the full qsize_d
  auto& sp = c.create<Homme::SimulationParams>();
  sp.qsize = std::max(HOMMEXX_QSIZE_D/2,1);

  // Set parameters
  constexpr int ne = 2;
  set_homme_param("ne",ne);

  // Create the grids
  ekat::ParameterList params;
  params.set<std::string>("physics_grid_type","GLL");
  params.set<std::string>("vertical_coordinate_filename","NONE");
  HommeGridsManager gm(comm,params);
  gm.build_grids();

  // Local counters
  const int num_local_elems = get_num_local_elems_f90();
  const int num_local_cols = get_num_local_columns_f90(pg_gll);
  EKAT_REQUIRE_MSG(num_local_cols>0, "Internal test error! Fix homme_pd_remap_tests, please.\n");

  // Get physics and dynamics grids, and their dofs
  auto phys_grid = gm.get_grid("Physics GLL");
  auto dyn_grid  = std::dynamic_pointer_cast<const SEGrid>(gm.get_grid("Dynamics"));
  auto h_p_dofs = phys_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  auto h_d_dofs = dyn_grid->get_cg_dofs_gids().get_view<const gid_t*,Host>();
  auto h_d_lid2idx = dyn_grid->get_lid_to_idx_map().get_view<const int**,Host>();

  // Get some dimensions for Homme
  constexpr int np  = HOMMEXX_NP;
  constexpr int NVL = HOMMEXX_NUM_PHYSICAL_LEV;
  constexpr int NTL = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NQ  = HOMMEXX_QSIZE_D;
  const int nle = num_local_elems;
  const int nlc = num_local_cols;
  const auto units = ekat::units::m;  // Placeholder units (we don't care about units here)

  const int n0 = IPDF(0,NTL-1)(engine);

  // Note on prefixes: s=scalar, v=vector, ss=scalar state, vs=vector_state, tr=tracer array

  // Create tags and dimensions
  std::vector<FieldTag> s_2d_dyn_tags   = {EL,          GP, GP    };
  std::vector<FieldTag> v_2d_dyn_tags   = {EL,     CMP, GP, GP    };
  std::vector<FieldTag> s_3d_dyn_tags   = {EL,          GP, GP, LEV};
  std::vector<FieldTag> v_3d_dyn_tags   = {EL,     CMP, GP, GP, LEV};
  std::vector<FieldTag> ss_3d_dyn_tags  = {EL, TL,      GP, GP, LEV};
  std::vector<FieldTag> vs_3d_dyn_tags  = {EL, TL, CMP, GP, GP, LEV};

  std::vector<FieldTag> s_2d_phys_tags  = {COL         };
  std::vector<FieldTag> v_2d_phys_tags  = {COL, CMP    };
  std::vector<FieldTag> s_3d_phys_tags  = {COL,      LEV};
  std::vector<FieldTag> v_3d_phys_tags  = {COL, CMP, LEV};
  std::vector<FieldTag> vs_3d_phys_tags = {COL, CMP, LEV};
  std::vector<FieldTag> ss_3d_phys_tags = {COL,      LEV};

  std::vector<int> s_2d_dyn_dims   = {nle,          np, np     };
  std::vector<int> v_2d_dyn_dims   = {nle,       2, np, np     };
  std::vector<int> s_3d_dyn_dims   = {nle,          np, np, NVL};
  std::vector<int> v_3d_dyn_dims   = {nle,       2, np, np, NVL};
  std::vector<int> ss_3d_dyn_dims  = {nle, NTL,     np, np, NVL};
  std::vector<int> vs_3d_dyn_dims  = {nle, NTL,  2, np, np, NVL};
  std::vector<int> tr_3d_dyn_dims  = {nle,      NQ, np, np, NVL};

  std::vector<int> s_2d_phys_dims  = {nlc         };
  std::vector<int> v_2d_phys_dims  = {nlc,  2     };
  std::vector<int> s_3d_phys_dims  = {nlc,     NVL};
  std::vector<int> v_3d_phys_dims  = {nlc,  2, NVL};
  std::vector<int> ss_3d_phys_dims = {nlc,     NVL};
  std::vector<int> vs_3d_phys_dims = {nlc,  2, NVL};
  std::vector<int> tr_3d_phys_dims = {nlc, NQ, NVL};

  // Create identifiers
  const auto dgn = dyn_grid->name();
  const auto pgn = phys_grid->name();
  FID s_2d_dyn_fid  ("s_2d_dyn", FL(s_2d_dyn_tags,  s_2d_dyn_dims), units, dgn);
  FID v_2d_dyn_fid  ("v_2d_dyn", FL(v_2d_dyn_tags,  v_2d_dyn_dims), units, dgn);
  FID s_3d_dyn_fid  ("s_3d_dyn", FL(s_3d_dyn_tags,  s_3d_dyn_dims), units, dgn);
  FID v_3d_dyn_fid  ("v_3d_dyn", FL(v_3d_dyn_tags,  v_3d_dyn_dims), units, dgn);
  FID ss_3d_dyn_fid ("ss_3d_dyn", FL(ss_3d_dyn_tags, ss_3d_dyn_dims),units, dgn);
  FID vs_3d_dyn_fid ("vs_3d_dyn", FL(vs_3d_dyn_tags, vs_3d_dyn_dims),units, dgn);
  FID tr_3d_dyn_fid ("tr_3d_dyn", FL(v_3d_dyn_tags, tr_3d_dyn_dims),units, dgn);

  FID s_2d_phys_fid ("s_2d_phys",  FL(s_2d_phys_tags, s_2d_phys_dims),units, pgn);
  FID v_2d_phys_fid ("v_2d_phys",  FL(v_2d_phys_tags, v_2d_phys_dims),units, pgn);
  FID s_3d_phys_fid ("s_3d_phys",  FL(s_3d_phys_tags, s_3d_phys_dims),units, pgn);
  FID v_3d_phys_fid ("v_3d_phys",  FL(v_3d_phys_tags, v_3d_phys_dims),units, pgn);
  FID ss_3d_phys_fid ("ss_3d_phys", FL(ss_3d_phys_tags, ss_3d_phys_dims),units, pgn);
  FID vs_3d_phys_fid ("vs_3d_phys", FL(vs_3d_phys_tags, vs_3d_phys_dims),units, pgn);
  FID tr_3d_phys_fid ("tr_3d_phys", FL(v_3d_phys_tags, tr_3d_phys_dims),units, pgn);

  // Create fields
  Field s_2d_field_phys (s_2d_phys_fid);
  Field v_2d_field_phys (v_2d_phys_fid);
  Field s_3d_field_phys (s_3d_phys_fid);
  Field v_3d_field_phys (v_3d_phys_fid);
  Field ss_3d_field_phys (ss_3d_phys_fid);
  Field vs_3d_field_phys (vs_3d_phys_fid);
  Field tr_3d_field_phys (tr_3d_phys_fid);

  Field s_2d_field_dyn(s_2d_dyn_fid);
  Field v_2d_field_dyn(v_2d_dyn_fid);
  Field s_3d_field_dyn(s_3d_dyn_fid);
  Field v_3d_field_dyn(v_3d_dyn_fid);
  Field ss_3d_field_dyn (ss_3d_dyn_fid);
  Field vs_3d_field_dyn (vs_3d_dyn_fid);
  Field tr_3d_field_dyn (tr_3d_dyn_fid);

  // Request allocation to fit packs of reals for 3d views
  s_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  v_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  ss_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  vs_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  tr_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);

  s_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  v_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  ss_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  vs_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  tr_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);

  // Allocate view
  s_2d_field_phys.allocate_view();
  v_2d_field_phys.allocate_view();
  s_3d_field_phys.allocate_view();
  v_3d_field_phys.allocate_view();
  ss_3d_field_phys.allocate_view();
  vs_3d_field_phys.allocate_view();
  tr_3d_field_phys.allocate_view();

  s_2d_field_dyn.allocate_view();
  v_2d_field_dyn.allocate_view();
  s_3d_field_dyn.allocate_view();
  v_3d_field_dyn.allocate_view();
  ss_3d_field_dyn.allocate_view();
  vs_3d_field_dyn.allocate_view();
  tr_3d_field_dyn.allocate_view();

  // Build the remapper, and register the fields
  std::shared_ptr<Remapper> remapper(new Remapper(phys_grid,dyn_grid));
  remapper->registration_begins();
  remapper->register_field(s_2d_field_phys, s_2d_field_dyn);
  remapper->register_field(v_2d_field_phys, v_2d_field_dyn);
  remapper->register_field(s_3d_field_phys, s_3d_field_dyn);
  remapper->register_field(v_3d_field_phys, v_3d_field_dyn);
  remapper->register_field(ss_3d_field_phys, ss_3d_field_dyn.subfield("cmp",1,n0));
  remapper->register_field(vs_3d_field_phys, vs_3d_field_dyn.subfield("cmp",1,n0));
  remapper->register_field(tr_3d_field_phys, tr_3d_field_dyn);
  remapper->registration_ends();

  SECTION ("remap") {

    for (bool fwd : {true, false}) {
      if (comm.am_i_root()) {
        std::cout << " -> Remap " << (fwd ? " forward\n" : " backward\n");
      }

      // Note: for the dyn->phys test to run correctly, the dynamics input v must be synced,
      //       meaning that the values at the interface between two elements must match.
      //       To do this, we initialize each entry in the dynamic v with the id
      //       of the corresponding physics column.
      //       But since this approach makes checking answers much easier, we use it also for phys->dyn.

      if (fwd) {
        auto h_s_2d_view = s_2d_field_phys.get_view<Homme::Real*,Host>();
        auto h_v_2d_view = v_2d_field_phys.get_view<Homme::Real**,Host>();
        auto h_s_3d_view = s_3d_field_phys.get_view<Homme::Real**,Host>();
        auto h_v_3d_view = v_3d_field_phys.get_view<Homme::Real***,Host>();
        auto h_ss_3d_view = ss_3d_field_phys.get_view<Homme::Real**,Host>();
        auto h_vs_3d_view = vs_3d_field_phys.get_view<Homme::Real***,Host>();
        auto h_tr_3d_view = tr_3d_field_phys.get_view<Homme::Real***,Host>();

        for (int idof=0; idof<num_local_cols; ++idof) {
          auto gid = h_p_dofs(idof);
          h_s_2d_view(idof) = gid;
          h_v_2d_view(idof,0) = gid;
          h_v_2d_view(idof,1) = gid;
          for (int il=0; il<NVL; ++ il) {
            h_s_3d_view(idof,il) = gid;
            h_v_3d_view(idof,0,il) = gid;
            h_v_3d_view(idof,1,il) = gid;

            h_ss_3d_view(idof,il) = gid;
            h_vs_3d_view(idof,0,il) = gid;
            h_vs_3d_view(idof,1,il) = gid;
            for (int iq=0; iq<NQ; ++iq) {
              h_tr_3d_view(idof,iq,il) = gid;
            }
          }
        }
        s_2d_field_phys.sync_to_dev();
        v_2d_field_phys.sync_to_dev();
        s_3d_field_phys.sync_to_dev();
        v_3d_field_phys.sync_to_dev();

        ss_3d_field_phys.sync_to_dev();
        vs_3d_field_phys.sync_to_dev();
        tr_3d_field_phys.sync_to_dev();
      } else {
        auto h_s_2d_view = s_2d_field_dyn.get_view<Homme::Real***,Host>();
        auto h_v_2d_view = v_2d_field_dyn.get_view<Homme::Real****,Host>();
        auto h_s_3d_view = s_3d_field_dyn.get_view<Homme::Real****,Host>();
        auto h_v_3d_view = v_3d_field_dyn.get_view<Homme::Real*****,Host>();
        auto h_ss_3d_view = ss_3d_field_dyn.get_view<Homme::Real*****,Host>();
        auto h_vs_3d_view = vs_3d_field_dyn.get_view<Homme::Real******,Host>();
        auto h_tr_3d_view = tr_3d_field_dyn.get_view<Homme::Real*****,Host>();

        for (int ie=0; ie<num_local_elems; ++ie) {
          for (int ip=0; ip<NP; ++ip) {
            for (int jp=0; jp<NP; ++jp) {
              const int idof = ie*NP*NP + ip*NP + jp;
              auto gid = h_d_dofs(idof);
              h_s_2d_view(ie,ip,jp) = gid;
              h_v_2d_view(ie,0,ip,jp) = gid;
              h_v_2d_view(ie,1,ip,jp) = gid;
              for (int il=0; il<NVL; ++ il) {
                h_s_3d_view(ie,ip,jp,il) = gid;
                h_v_3d_view(ie,0,ip,jp,il) = gid;
                h_v_3d_view(ie,1,ip,jp,il) = gid;

                for (int itl=0; itl<NTL; ++itl) {
                  h_ss_3d_view(ie,itl,ip,jp,il) = gid;
                  h_vs_3d_view(ie,itl,0,ip,jp,il) = gid;
                  h_vs_3d_view(ie,itl,1,ip,jp,il) = gid;
                }
                for (int iq=0; iq<NQ; ++iq) {
                  h_tr_3d_view(ie,iq,ip,jp,il) = gid;
                }
              }
            }
          }
        }
        s_2d_field_dyn.sync_to_dev();
        v_2d_field_dyn.sync_to_dev();
        s_3d_field_dyn.sync_to_dev();
        v_3d_field_dyn.sync_to_dev();

        ss_3d_field_dyn.sync_to_dev();
        vs_3d_field_dyn.sync_to_dev();
        tr_3d_field_dyn.sync_to_dev();
      }

      // Remap
      remapper->remap(fwd);

      // Check
      {
        // 2d scalar
        s_2d_field_phys.sync_to_host();
        s_2d_field_dyn.sync_to_host();
        auto phys = s_2d_field_phys.get_view<Homme::Real*,Host>();
        auto dyn  = s_2d_field_dyn.get_view<Homme::Real***,Host>();

        if (fwd) {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            if (dyn(ie,ip,jp)!=h_d_dofs(idof)) {
                printf(" ** 2D Scalar ** \n");
                printf("d_out(%d,%d,%d): %2.16f\n",ie,ip,jp,dyn(ie,ip,jp));
                printf("expected: = %d\n",h_d_dofs(idof));
            }
            REQUIRE (dyn(ie,ip,jp)==h_d_dofs(idof));
          }
        } else {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            if (phys(idof)!=h_p_dofs(idof)) {
                printf(" ** 2D Scalar ** \n");
                printf("  p_out(%d) = %2.16f\n",idof,phys(idof));
                printf("  expected: = %d\n",h_p_dofs(idof));
            }
            REQUIRE (phys(idof)==h_p_dofs(idof));
          }
        }
      }

      {
        // 2d vector
        v_2d_field_phys.sync_to_host();
        v_2d_field_dyn.sync_to_host();
        auto phys = v_2d_field_phys.get_view<Homme::Real**,Host>();
        auto dyn  = v_2d_field_dyn.get_view<Homme::Real****,Host>();
        if (fwd) {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int icomp=0; icomp<2; ++icomp) {
              if (dyn(ie,icomp,ip,jp)!=h_d_dofs(idof)) {
                  printf(" ** 2D Vector ** \n");
                  printf("d_out(%d,%d,%d,%d): %2.16f\n",ie,ip,jp,icomp,dyn(ie,icomp,ip,jp));
                  printf("expected: = %d\n",h_d_dofs(idof));
              }
              REQUIRE (dyn(ie,icomp,ip,jp)==h_d_dofs(idof));
            }
          }
        } else {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int icomp=0; icomp<2; ++icomp) {
              if (phys(idof,icomp)!=h_p_dofs(idof)) {
                  printf(" ** 2D Vector ** \n");
                  printf("p_out(%d, %d) = %2.16f\n",idof,icomp,phys(idof,icomp));
                  printf("expected: = %d\n",h_p_dofs(idof));
              }
              REQUIRE (phys(idof,icomp)==h_p_dofs(idof));
            }
          }
        }
      }

      {
        // 3d scalar
        s_3d_field_phys.sync_to_host();
        s_3d_field_dyn.sync_to_host();
        auto phys = s_3d_field_phys.get_view<Homme::Real**,Host>();
        auto dyn  = s_3d_field_dyn.get_view<Homme::Real****,Host>();
        if (fwd) {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (dyn(ie,ip,jp,ilev)!=h_d_dofs(idof)) {
                  printf(" ** 3D Scalar ** \n");
                  printf("d_out(%d,%d,%d,%d): %2.16f\n",ie,ip,jp,ilev,dyn(ie,ip,jp,ilev));
                  printf("expected: = %d\n",h_d_dofs(idof));
              }
              REQUIRE (dyn(ie,ip,jp,ilev)==h_d_dofs(idof));
            }
          }
        } else {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (phys(idof,ilev)!=h_p_dofs(idof)) {
                  printf(" ** 3D Scalar ** \n");
                  printf("p_out(%d,%d) = %2.16f\n",idof,ilev,phys(idof,ilev));
                  printf("expected: = %d\n",h_p_dofs(idof));
              }
              REQUIRE (phys(idof,ilev)==h_p_dofs(idof));
            }
          }
        }
      }

      {
        // 3d vector
        v_3d_field_phys.sync_to_host();
        v_3d_field_dyn.sync_to_host();
        auto phys = v_3d_field_phys.get_view<Homme::Real***,Host>();
        auto dyn  = v_3d_field_dyn.get_view<Homme::Real*****,Host>();
        if (fwd) {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (dyn(ie,icomp,ip,jp,ilev)!=h_d_dofs(idof)) {
                    printf(" ** 3D Vector ** \n");
                    printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,icomp,ip,jp,ilev,dyn(ie,icomp,ip,jp,ilev));
                    printf("expected: = %d\n",h_d_dofs(idof));
                }
                REQUIRE (dyn(ie,icomp,ip,jp,ilev)==h_d_dofs(idof));
              }
            }
          }
        } else {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (phys(idof,icomp,ilev)!=h_p_dofs(idof)) {
                    printf(" ** 3D Vector ** \n");
                    printf("p_out(%d,%d,%d) = %2.16f\n",idof,icomp,ilev,phys(idof,icomp,ilev));
                    printf("expected: = %d\n",h_p_dofs(idof));
                }
                REQUIRE (phys(idof,icomp,ilev)==h_p_dofs(idof));
              }
            }
          }
        }
      }

      {
        // 3d scalar state
        ss_3d_field_phys.sync_to_host();
        ss_3d_field_dyn.sync_to_host();
        auto phys = ss_3d_field_phys.get_view<Homme::Real**,Host>();
        auto dyn  = ss_3d_field_dyn.get_view<Homme::Real*****,Host>();
        if (fwd) {
          for (int ilev=0; ilev<NVL; ++ilev) {
            for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
              int ie = h_d_lid2idx(idof,0);
              int ip = h_d_lid2idx(idof,1);
              int jp = h_d_lid2idx(idof,2);
              auto gid = h_d_dofs(idof);
              if (dyn(ie,n0,ip,jp,ilev)!=gid) {
                  printf(" ** 3D Scalar State ** \n");
                  printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,n0,ip,jp,ilev,dyn(ie,n0,ip,jp,ilev));
                  printf("expected: = %d\n",gid);
              }
              REQUIRE (dyn(ie,n0,ip,jp,ilev)==gid);
            }
          }
        } else {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (phys(idof,ilev)!=h_p_dofs(idof)) {
                  printf(" ** 3D Scalar State ** \n");
                  printf("p_out(%d,%d) = %2.16f\n",idof,ilev,phys(idof,ilev));
                  printf("expected: = %d\n",h_p_dofs(idof));
              }
              REQUIRE (phys(idof,ilev)==h_p_dofs(idof));
            }
          }
        }
      }

      {
        // 3d vector state
        vs_3d_field_phys.sync_to_host();
        vs_3d_field_dyn.sync_to_host();
        auto phys = vs_3d_field_phys.get_view<Homme::Real***,Host>();
        auto dyn  = vs_3d_field_dyn.get_view<Homme::Real******,Host>();
        if (fwd) {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (dyn(ie,n0,icomp,ip,jp,ilev)!=h_d_dofs(idof)) {
                    printf(" ** 3D Vector State ** \n");
                    printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,icomp,ip,jp,ilev,dyn(ie,n0,icomp,ip,jp,ilev));
                    printf("expected: = %d\n",h_d_dofs(idof));
                }
                REQUIRE (dyn(ie,n0,icomp,ip,jp,ilev)==h_d_dofs(idof));
              }
            }
          }
        } else {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (phys(idof,icomp,ilev)!=h_p_dofs(idof)) {
                    printf(" ** 3D Vector State ** \n");
                    printf("p_out(%d,%d,%d) = %2.16f\n",idof,icomp,ilev,phys(idof,icomp,ilev));
                    printf("expected: = %d\n",h_p_dofs(idof));
                }
                REQUIRE (phys(idof,icomp,ilev)==h_p_dofs(idof));
              }
            }
          }
        }
      }

      {
        // 3d tracers
        tr_3d_field_phys.sync_to_host();
        tr_3d_field_dyn.sync_to_host();
        auto phys = tr_3d_field_phys.get_view<Homme::Real***,Host>();
        auto dyn  = tr_3d_field_dyn.get_view<Homme::Real*****,Host>();
        if (fwd) {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int iq=0; iq<2; ++iq) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (dyn(ie,iq,ip,jp,ilev)!=h_d_dofs(idof)) {
                    printf(" ** 3D Tracer State ** \n");
                    printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,iq,ip,jp,ilev,dyn(ie,iq,ip,jp,ilev));
                    printf("expected: = %d\n",h_d_dofs(idof));
                }
                REQUIRE (dyn(ie,iq,ip,jp,ilev)==h_d_dofs(idof));
              }
            }
          }
        } else {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int iq=0; iq<2; ++iq) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (phys(idof,iq,ilev)!=h_p_dofs(idof)) {
                    printf(" ** 3D Tracer State ** \n");
                    printf("p_out(%d,%d,%d) = %2.16f\n",idof,iq,ilev,phys(idof,iq,ilev));
                    printf("expected: = %d\n",h_p_dofs(idof));
                }
                REQUIRE (phys(idof,iq,ilev)==h_p_dofs(idof));
              }
            }
          }
        }
      }
    }
  }

  // Delete remapper before finalizing the mpi context, since the remapper has some MPI stuff in it
  remapper = nullptr;

  // Finalize Homme::Context
  Homme::Context::finalize_singleton();

  // Cleanup f90 structures
  cleanup_test_f90();
}

TEST_CASE("combo_remap", "") {

  using namespace scream;
  using namespace ShortFieldTagsNames;

  // Some type defs
  using Remapper = PhysicsDynamicsRemapper;
  using IPDF = std::uniform_int_distribution<int>;
  using FID = FieldIdentifier;
  using FL  = FieldLayout;

  constexpr int pg_gll = 0;
  constexpr int PackSize = HOMMEXX_VECTOR_SIZE;

  // Create a comm
  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test(&comm);

  // Init homme context
  if (!is_parallel_inited_f90()) {
    auto comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
    init_parallel_f90(comm_f);
  }
  init_test_params_f90 ();

  // We'll use this extensively, so let's use a short ref name
  auto& c = Homme::Context::singleton();

  // Set a value for qsize that is not the full qsize_d
  auto& sp = c.create<Homme::SimulationParams>();
  sp.qsize = std::max(HOMMEXX_QSIZE_D/2,1);

  // Set parameters
  constexpr int ne = 2;
  set_homme_param("ne",ne);

  // Create the grids
  ekat::ParameterList params;
  params.set<std::string>("physics_grid_type","GLL");
  params.set<std::string>("vertical_coordinate_filename","NONE");
  HommeGridsManager gm(comm,params);
  gm.build_grids();

  // Local counters
  const int num_local_elems = get_num_local_elems_f90();
  const int num_local_cols = get_num_local_columns_f90(pg_gll);
  EKAT_REQUIRE_MSG(num_local_cols>0, "Internal test error! Fix homme_pd_remap_tests, please.\n");

  // Get physics and dynamics grids, and their dofs
  auto phys_grid = gm.get_grid("Physics GLL");
  auto dyn_grid  = std::dynamic_pointer_cast<const SEGrid>(gm.get_grid("Dynamics"));
  auto h_p_dofs = phys_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  auto h_d_dofs = dyn_grid->get_cg_dofs_gids().get_view<const gid_t*,Host>();
  auto h_d_lid2idx = dyn_grid->get_lid_to_idx_map().get_view<const int**,Host>();

  // Get some dimensions for Homme
  constexpr int np  = HOMMEXX_NP;
  constexpr int NVL = HOMMEXX_NUM_PHYSICAL_LEV;
  constexpr int NTL = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NQ  = HOMMEXX_QSIZE_D;
  const int nle = num_local_elems;
  const int nlc = num_local_cols;
  const auto units = ekat::units::m;  // Placeholder units (we don't care about units here)

  const int n0 = IPDF(0,NTL-1)(engine);

  // Note on prefixes: s=scalar, v=vector, ss=scalar state, vs=vector_state, tr=tracer array

  // Create tags and dimensions
  std::vector<FieldTag> s_2d_dyn_tags   = {EL,          GP, GP     };
  std::vector<FieldTag> v_2d_dyn_tags   = {EL,     CMP, GP, GP     };
  std::vector<FieldTag> s_3d_dyn_tags   = {EL,          GP, GP, LEV};
  std::vector<FieldTag> v_3d_dyn_tags   = {EL,     CMP, GP, GP, LEV};
  std::vector<FieldTag> ss_3d_dyn_tags  = {EL, TL,      GP, GP, LEV};
  std::vector<FieldTag> vs_3d_dyn_tags  = {EL, TL, CMP, GP, GP, LEV};
  std::vector<FieldTag> tr_3d_dyn_tags  = {EL,     CMP, GP, GP, LEV};

  std::vector<FieldTag> s_2d_phys_tags  = {COL          };
  std::vector<FieldTag> v_2d_phys_tags  = {COL, CMP     };
  std::vector<FieldTag> s_3d_phys_tags  = {COL,      LEV};
  std::vector<FieldTag> v_3d_phys_tags  = {COL, CMP, LEV};
  std::vector<FieldTag> vs_3d_phys_tags = {COL, CMP, LEV};
  std::vector<FieldTag> ss_3d_phys_tags = {COL,      LEV};
  std::vector<FieldTag> tr_3d_phys_tags = {COL, CMP, LEV};

  std::vector<int> s_2d_dyn_dims   = {nle,          np, np     };
  std::vector<int> v_2d_dyn_dims   = {nle,       2, np, np     };
  std::vector<int> s_3d_dyn_dims   = {nle,          np, np, NVL};
  std::vector<int> v_3d_dyn_dims   = {nle,       2, np, np, NVL};
  std::vector<int> ss_3d_dyn_dims  = {nle, NTL,     np, np, NVL};
  std::vector<int> vs_3d_dyn_dims  = {nle, NTL,  2, np, np, NVL};
  std::vector<int> tr_3d_dyn_dims  = {nle,      NQ, np, np, NVL};

  std::vector<int> s_2d_phys_dims  = {nlc         };
  std::vector<int> v_2d_phys_dims  = {nlc,  2     };
  std::vector<int> s_3d_phys_dims  = {nlc,     NVL};
  std::vector<int> v_3d_phys_dims  = {nlc,  2, NVL};
  std::vector<int> ss_3d_phys_dims = {nlc,     NVL};
  std::vector<int> vs_3d_phys_dims = {nlc,  2, NVL};
  std::vector<int> tr_3d_phys_dims = {nlc, NQ, NVL};

  // Create identifiers
  const auto dgn = dyn_grid->name();
  const auto pgn = phys_grid->name();
  FID s_2d_dyn_fid  ("s_2d_dyn", FL(s_2d_dyn_tags,  s_2d_dyn_dims), units, dgn);
  FID v_2d_dyn_fid  ("v_2d_dyn", FL(v_2d_dyn_tags,  v_2d_dyn_dims), units, dgn);
  FID s_3d_dyn_fid  ("s_3d_dyn", FL(s_3d_dyn_tags,  s_3d_dyn_dims), units, dgn);
  FID v_3d_dyn_fid  ("v_3d_dyn", FL(v_3d_dyn_tags,  v_3d_dyn_dims), units, dgn);
  FID ss_3d_dyn_fid ("ss_3d_dyn", FL(ss_3d_dyn_tags, ss_3d_dyn_dims),units, dgn);
  FID vs_3d_dyn_fid ("vs_3d_dyn", FL(vs_3d_dyn_tags, vs_3d_dyn_dims),units, dgn);
  FID tr_3d_dyn_fid ("tr_3d_dyn", FL(tr_3d_dyn_tags, tr_3d_dyn_dims),units, dgn);

  FID s_2d_phys_fid ("s_2d_phys",  FL(s_2d_phys_tags, s_2d_phys_dims),units, pgn);
  FID v_2d_phys_fid ("v_2d_phys",  FL(v_2d_phys_tags, v_2d_phys_dims),units, pgn);
  FID s_3d_phys_fid ("s_3d_phys",  FL(s_3d_phys_tags, s_3d_phys_dims),units, pgn);
  FID v_3d_phys_fid ("v_3d_phys",  FL(v_3d_phys_tags, v_3d_phys_dims),units, pgn);
  FID ss_3d_phys_fid ("ss_3d_phys", FL(ss_3d_phys_tags, ss_3d_phys_dims),units, pgn);
  FID vs_3d_phys_fid ("vs_3d_phys", FL(vs_3d_phys_tags, vs_3d_phys_dims),units, pgn);
  FID tr_3d_phys_fid ("tr_3d_phys", FL(tr_3d_phys_tags, tr_3d_phys_dims),units, pgn);

  // Create fields
  Field s_2d_field_phys (s_2d_phys_fid);
  Field v_2d_field_phys (v_2d_phys_fid);
  Field s_3d_field_phys (s_3d_phys_fid);
  Field v_3d_field_phys (v_3d_phys_fid);
  Field ss_3d_field_phys (ss_3d_phys_fid);
  Field vs_3d_field_phys (vs_3d_phys_fid);
  Field tr_3d_field_phys (tr_3d_phys_fid);

  Field s_2d_field_dyn(s_2d_dyn_fid);
  Field v_2d_field_dyn(v_2d_dyn_fid);
  Field s_3d_field_dyn(s_3d_dyn_fid);
  Field v_3d_field_dyn(v_3d_dyn_fid);
  Field ss_3d_field_dyn (ss_3d_dyn_fid);
  Field vs_3d_field_dyn (vs_3d_dyn_fid);
  Field tr_3d_field_dyn (tr_3d_dyn_fid);

  // Request allocation to fit packs of reals for 3d views
  s_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  v_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  ss_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  vs_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);
  tr_3d_field_phys.get_header().get_alloc_properties().request_allocation(PackSize);

  s_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  v_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  ss_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  vs_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);
  tr_3d_field_dyn.get_header().get_alloc_properties().request_allocation(PackSize);

  // Allocate view
  s_2d_field_phys.allocate_view();
  v_2d_field_phys.allocate_view();
  s_3d_field_phys.allocate_view();
  v_3d_field_phys.allocate_view();
  ss_3d_field_phys.allocate_view();
  vs_3d_field_phys.allocate_view();
  tr_3d_field_phys.allocate_view();

  s_2d_field_dyn.allocate_view();
  v_2d_field_dyn.allocate_view();
  s_3d_field_dyn.allocate_view();
  v_3d_field_dyn.allocate_view();
  ss_3d_field_dyn.allocate_view();
  vs_3d_field_dyn.allocate_view();
  tr_3d_field_dyn.allocate_view();

  // Build the remapper, and register the fields
  std::shared_ptr<Remapper> remapper(new Remapper(phys_grid,dyn_grid));
  remapper->registration_begins();
  remapper->register_field(s_2d_field_phys, s_2d_field_dyn);
  remapper->register_field(v_2d_field_phys, v_2d_field_dyn);
  remapper->register_field(s_3d_field_phys, s_3d_field_dyn);
  remapper->register_field(v_3d_field_phys, v_3d_field_dyn);
  remapper->register_field(ss_3d_field_phys, ss_3d_field_dyn.subfield("cmp",1,n0));
  remapper->register_field(vs_3d_field_phys, vs_3d_field_dyn.subfield("cmp",1,n0));
  remapper->register_field(tr_3d_field_phys, tr_3d_field_dyn);
  remapper->registration_ends();

  SECTION ("combo_test") {

    for (bool pdp : {true, false}) {
      if (comm.am_i_root()) {
        std::cout << " -> Remap " << (pdp ? " phys->dyn->pys\n" : " dyn->phys->dyn\n");
      }

      // Note: for the dyn->phys test to run correctly, the dynamics input v must be synced,
      //       meaning that the values at the interface between two elements must match.
      //       To do this, we initialize each entry in the dynamic v with the id
      //       of the corresponding column.
      //       But since this approach makes checking answers much easier, we use it also for phys->dyn.

      if (pdp) {
        auto h_s_2d_view = s_2d_field_phys.get_view<Homme::Real*,Host>();
        auto h_v_2d_view = v_2d_field_phys.get_view<Homme::Real**,Host>();
        auto h_s_3d_view = s_3d_field_phys.get_view<Homme::Real**,Host>();
        auto h_v_3d_view = v_3d_field_phys.get_view<Homme::Real***,Host>();
        auto h_ss_3d_view = ss_3d_field_phys.get_view<Homme::Real**,Host>();
        auto h_vs_3d_view = vs_3d_field_phys.get_view<Homme::Real***,Host>();
        auto h_tr_3d_view = tr_3d_field_phys.get_view<Homme::Real***,Host>();

        for (int idof=0; idof<num_local_cols; ++idof) {
          auto gid = h_p_dofs(idof);
          h_s_2d_view(idof) = gid;
          h_v_2d_view(idof,0) = gid;
          h_v_2d_view(idof,1) = gid;
          for (int il=0; il<NVL; ++ il) {
            h_s_3d_view(idof,il) = gid;
            h_v_3d_view(idof,0,il) = gid;
            h_v_3d_view(idof,1,il) = gid;

            h_ss_3d_view(idof,il) = gid;
            h_vs_3d_view(idof,0,il) = gid;
            h_vs_3d_view(idof,1,il) = gid;
            for (int iq=0; iq<NQ; ++iq) {
              h_tr_3d_view(idof,iq,il) = gid;
            }
          }
        }
        s_2d_field_phys.sync_to_dev();
        v_2d_field_phys.sync_to_dev();
        s_3d_field_phys.sync_to_dev();
        v_3d_field_phys.sync_to_dev();

        ss_3d_field_phys.sync_to_dev();
        vs_3d_field_phys.sync_to_dev();
        tr_3d_field_phys.sync_to_dev();
      } else {
        auto h_s_2d_view = s_2d_field_dyn.get_view<Homme::Real***,Host>();
        auto h_v_2d_view = v_2d_field_dyn.get_view<Homme::Real****,Host>();
        auto h_s_3d_view = s_3d_field_dyn.get_view<Homme::Real****,Host>();
        auto h_v_3d_view = v_3d_field_dyn.get_view<Homme::Real*****,Host>();
        auto h_ss_3d_view = ss_3d_field_dyn.get_view<Homme::Real*****,Host>();
        auto h_vs_3d_view = vs_3d_field_dyn.get_view<Homme::Real******,Host>();
        auto h_tr_3d_view = tr_3d_field_dyn.get_view<Homme::Real*****,Host>();

        for (int ie=0; ie<num_local_elems; ++ie) {
          for (int ip=0; ip<NP; ++ip) {
            for (int jp=0; jp<NP; ++jp) {
              const int idof = ie*NP*NP + ip*NP + jp;
              auto gid = h_d_dofs(idof);
              h_s_2d_view(ie,ip,jp) = gid;
              h_v_2d_view(ie,0,ip,jp) = gid;
              h_v_2d_view(ie,1,ip,jp) = gid;
              for (int il=0; il<NVL; ++ il) {
                h_s_3d_view(ie,ip,jp,il) = gid;
                h_v_3d_view(ie,0,ip,jp,il) = gid;
                h_v_3d_view(ie,1,ip,jp,il) = gid;

                for (int itl=0; itl<NTL; ++itl) {
                  h_ss_3d_view(ie,itl,ip,jp,il) = gid;
                  h_vs_3d_view(ie,itl,0,ip,jp,il) = gid;
                  h_vs_3d_view(ie,itl,1,ip,jp,il) = gid;
                }
                for (int iq=0; iq<NQ; ++iq) {
                  h_tr_3d_view(ie,iq,ip,jp,il) = gid;
                }
              }
            }
          }
        }
        s_2d_field_dyn.sync_to_dev();
        v_2d_field_dyn.sync_to_dev();
        s_3d_field_dyn.sync_to_dev();
        v_3d_field_dyn.sync_to_dev();

        ss_3d_field_dyn.sync_to_dev();
        vs_3d_field_dyn.sync_to_dev();
        tr_3d_field_dyn.sync_to_dev();
      }

      // Remap
      if (pdp) {
        remapper->remap(true);
        Kokkos::fence();
        remapper->remap(false);
      } else {
        remapper->remap(false);
        remapper->remap(true);
      }

      // Check
      {
        // 2d scalar
        s_2d_field_phys.sync_to_host();
        s_2d_field_dyn.sync_to_host();
        auto phys = s_2d_field_phys.get_view<Homme::Real*,Host>();
        auto dyn  = s_2d_field_dyn.get_view<Homme::Real***,Host>();

        if (pdp) {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            if (phys(idof)!=h_p_dofs(idof)) {
                printf(" ** 2D Scalar ** \n");
                printf("  p_out(%d) = %2.16f\n",idof,phys(idof));
                printf("  expected: = %d\n",h_p_dofs(idof));
            }
            REQUIRE (phys(idof)==h_p_dofs(idof));
          }
        } else {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            if (dyn(ie,ip,jp)!=h_d_dofs(idof)) {
                printf(" ** 2D Scalar ** \n");
                printf("d_out(%d,%d,%d): %2.16f\n",ie,ip,jp,dyn(ie,ip,jp));
                printf("expected: = %d\n",h_d_dofs(idof));
            }
            REQUIRE (dyn(ie,ip,jp)==h_d_dofs(idof));
          }
        }
      }

      {
        // 2d vector
        v_2d_field_phys.sync_to_host();
        v_2d_field_dyn.sync_to_host();
        auto phys = v_2d_field_phys.get_view<Homme::Real**,Host>();
        auto dyn  = v_2d_field_dyn.get_view<Homme::Real****,Host>();
        if (pdp) {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int icomp=0; icomp<2; ++icomp) {
              if (phys(idof,icomp)!=h_p_dofs(idof)) {
                  printf(" ** 2D Vector ** \n");
                  printf("p_out(%d, %d) = %2.16f\n",idof,icomp,phys(idof,icomp));
                  printf("expected: = %d\n",h_p_dofs(idof));
              }
              REQUIRE (phys(idof,icomp)==h_p_dofs(idof));
            }
          }
        } else {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int icomp=0; icomp<2; ++icomp) {
              if (dyn(ie,icomp,ip,jp)!=h_d_dofs(idof)) {
                  printf(" ** 2D Vector ** \n");
                  printf("d_out(%d,%d,%d,%d): %2.16f\n",ie,ip,jp,icomp,dyn(ie,icomp,ip,jp));
                  printf("expected: = %d\n",h_d_dofs(idof));
              }
              REQUIRE (dyn(ie,icomp,ip,jp)==h_d_dofs(idof));
            }
          }
        }
      }

      {
        // 3d scalar
        s_3d_field_phys.sync_to_host();
        s_3d_field_dyn.sync_to_host();
        auto phys = s_3d_field_phys.get_view<Homme::Real**,Host>();
        auto dyn  = s_3d_field_dyn.get_view<Homme::Real****,Host>();
        if (pdp) {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (phys(idof,ilev)!=h_p_dofs(idof)) {
                  printf(" ** 3D Scalar ** \n");
                  printf("p_out(%d,%d) = %2.16f\n",idof,ilev,phys(idof,ilev));
                  printf("expected: = %d\n",h_p_dofs(idof));
              }
              REQUIRE (phys(idof,ilev)==h_p_dofs(idof));
            }
          }
        } else {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (dyn(ie,ip,jp,ilev)!=h_d_dofs(idof)) {
                  printf(" ** 3D Scalar ** \n");
                  printf("d_out(%d,%d,%d,%d): %2.16f\n",ie,ip,jp,ilev,dyn(ie,ip,jp,ilev));
                  printf("expected: = %d\n",h_d_dofs(idof));
              }
              REQUIRE (dyn(ie,ip,jp,ilev)==h_d_dofs(idof));
            }
          }
        }
      }

      {
        // 3d vector
        v_3d_field_phys.sync_to_host();
        v_3d_field_dyn.sync_to_host();
        auto phys = v_3d_field_phys.get_view<Homme::Real***,Host>();
        auto dyn  = v_3d_field_dyn.get_view<Homme::Real*****,Host>();
        if (pdp) {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (phys(idof,icomp,ilev)!=h_p_dofs(idof)) {
                    printf(" ** 3D Vector ** \n");
                    printf("p_out(%d,%d,%d) = %2.16f\n",idof,icomp,ilev,phys(idof,icomp,ilev));
                    printf("expected: = %d\n",h_p_dofs(idof));
                }
                REQUIRE (phys(idof,icomp,ilev)==h_p_dofs(idof));
              }
            }
          }
        } else {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (dyn(ie,icomp,ip,jp,ilev)!=h_d_dofs(idof)) {
                    printf(" ** 3D Vector ** \n");
                    printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,icomp,ip,jp,ilev,dyn(ie,icomp,ip,jp,ilev));
                    printf("expected: = %d\n",h_d_dofs(idof));
                }
                REQUIRE (dyn(ie,icomp,ip,jp,ilev)==h_d_dofs(idof));
              }
            }
          }
        }
      }

      {
        // 3d scalar state
        ss_3d_field_phys.sync_to_host();
        ss_3d_field_dyn.sync_to_host();
        auto phys = ss_3d_field_phys.get_view<Homme::Real**,Host>();
        auto dyn  = ss_3d_field_dyn.get_view<Homme::Real*****,Host>();
        if (pdp) {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (phys(idof,ilev)!=h_p_dofs(idof)) {
                  printf(" ** 3D Scalar State ** \n");
                  printf("p_out(%d,%d) = %2.16f\n",idof,ilev,phys(idof,ilev));
                  printf("expected: = %d\n",h_p_dofs(idof));
              }
              REQUIRE (phys(idof,ilev)==h_p_dofs(idof));
            }
          }
        } else {
          for (int ilev=0; ilev<NVL; ++ilev) {
            for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
              int ie = h_d_lid2idx(idof,0);
              int ip = h_d_lid2idx(idof,1);
              int jp = h_d_lid2idx(idof,2);
              auto gid = h_d_dofs(idof);
              if (dyn(ie,n0,ip,jp,ilev)!=gid) {
                  printf(" ** 3D Scalar State ** \n");
                  printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,n0,ip,jp,ilev,dyn(ie,n0,ip,jp,ilev));
                  printf("expected: = %d\n",gid);
              }
              REQUIRE (dyn(ie,n0,ip,jp,ilev)==gid);
            }
          }
        }
      }

      {
        // 3d vector state
        vs_3d_field_phys.sync_to_host();
        vs_3d_field_dyn.sync_to_host();
        auto phys = vs_3d_field_phys.get_view<Homme::Real***,Host>();
        auto dyn  = vs_3d_field_dyn.get_view<Homme::Real******,Host>();
        if (pdp) {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (phys(idof,icomp,ilev)!=h_p_dofs(idof)) {
                    printf(" ** 3D Vector State ** \n");
                    printf("p_out(%d,%d,%d) = %2.16f\n",idof,icomp,ilev,phys(idof,icomp,ilev));
                    printf("expected: = %d\n",h_p_dofs(idof));
                }
                REQUIRE (phys(idof,icomp,ilev)==h_p_dofs(idof));
              }
            }
          }
        } else {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int icomp=0; icomp<2; ++icomp) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (dyn(ie,n0,icomp,ip,jp,ilev)!=h_d_dofs(idof)) {
                    printf(" ** 3D Vector State ** \n");
                    printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,icomp,ip,jp,ilev,dyn(ie,n0,icomp,ip,jp,ilev));
                    printf("expected: = %d\n",h_d_dofs(idof));
                }
                REQUIRE (dyn(ie,n0,icomp,ip,jp,ilev)==h_d_dofs(idof));
              }
            }
          }
        }
      }

      {
        // 3d tracers
        tr_3d_field_phys.sync_to_host();
        tr_3d_field_dyn.sync_to_host();
        auto phys = tr_3d_field_phys.get_view<Homme::Real***,Host>();
        auto dyn  = tr_3d_field_dyn.get_view<Homme::Real*****,Host>();
        if (pdp) {
          for (int idof=0; idof<phys_grid->get_num_local_dofs(); ++idof) {
            for (int iq=0; iq<2; ++iq) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (phys(idof,iq,ilev)!=h_p_dofs(idof)) {
                    printf(" ** 3D Tracer State ** \n");
                    printf("p_out(%d,%d,%d) = %2.16f\n",idof,iq,ilev,phys(idof,iq,ilev));
                    printf("expected: = %d\n",h_p_dofs(idof));
                }
                REQUIRE (phys(idof,iq,ilev)==h_p_dofs(idof));
              }
            }
          }
        } else {
          for (int idof=0; idof<dyn_grid->get_num_local_dofs(); ++idof) {
            int ie = h_d_lid2idx(idof,0);
            int ip = h_d_lid2idx(idof,1);
            int jp = h_d_lid2idx(idof,2);
            for (int iq=0; iq<2; ++iq) {
              for (int ilev=0; ilev<NVL; ++ilev) {
                if (dyn(ie,iq,ip,jp,ilev)!=h_d_dofs(idof)) {
                    printf(" ** 3D Tracer State ** \n");
                    printf("d_out(%d,%d,%d,%d,%d): %2.16f\n",ie,iq,ip,jp,ilev,dyn(ie,iq,ip,jp,ilev));
                    printf("expected: = %d\n",h_d_dofs(idof));
                }
                REQUIRE (dyn(ie,iq,ip,jp,ilev)==h_d_dofs(idof));
              }
            }
          }
        }
      }
    }
  }

  // Delete remapper before finalizing the mpi context, since the remapper has some MPI stuff in it
  remapper = nullptr;

  // Finalize Homme::Context
  Homme::Context::finalize_singleton();

  // Cleanup f90 structures
  cleanup_test_f90();
}

} // anonymous namespace
