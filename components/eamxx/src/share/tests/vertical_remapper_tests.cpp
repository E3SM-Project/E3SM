#include <catch2/catch.hpp>

#include "share/grid/remap/vertical_remapper.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace scream {

template<typename ViewT>
typename ViewT::HostMirror
cmvc (const ViewT& v) {
  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  return vh;
}

class VerticalRemapperTester : public VerticalRemapper {
public:
  VerticalRemapperTester (const grid_ptr_type& src_grid,
                          const std::string&   map_file,
                          const Field&         lev_prof,
                          const Field&         ilev_prof,
                          const Real           mask_val)
   : VerticalRemapper(src_grid, map_file, lev_prof, ilev_prof, mask_val)
  {
    // Nothing to do
  }
};

template<typename ViewT>
bool contains (const ViewT& v, const typename ViewT::traits::value_type& entry) {
  const auto vh = cmvc (v);
  const auto beg = vh.data();
  const auto end = vh.data() + vh.size();
  for (auto it=beg; it!=end; ++it) {
    if (*it == entry) {
      return true;
    }
  }
  return false;
}

void print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

// Helper function to create a grid given the number of dof's and a comm group.
std::shared_ptr<AbstractGrid>
build_src_grid(const ekat::Comm& comm, const int nldofs_src, const int nlevs_src) 
{
  auto src_grid = std::make_shared<PointGrid>("src",nldofs_src,nlevs_src,comm);

  auto src_dofs = src_grid->get_dofs_gids();
  auto src_dofs_h = src_dofs.get_view<gid_t*,Host>();
  std::iota(src_dofs_h.data(),src_dofs_h.data()+nldofs_src,nldofs_src*comm.rank());
  src_dofs.sync_to_dev();

  return src_grid;
}

// Helper function to create fields
Field
create_field(const std::string& name, const std::shared_ptr<const AbstractGrid>& grid, const bool twod, const bool vec, const bool mid = false, const int ps = 1)
{
  constexpr int vec_dim = 3;
  constexpr auto CMP = FieldTag::Component;
  constexpr auto units = ekat::units::Units::nondimensional();
  auto fl = twod
          ? (vec ? grid->get_2d_vector_layout (CMP,vec_dim)
                 : grid->get_2d_scalar_layout ())
          : (vec ? grid->get_3d_vector_layout (mid,CMP,vec_dim)
                 : grid->get_3d_scalar_layout (mid));
  FieldIdentifier fid(name,fl,units,grid->name());
  Field f(fid);
  f.get_header().get_alloc_properties().request_allocation(ps);
  f.allocate_view();
  return f;
}

// Helper function to set variable data
Real data_func(const int col, const int vec, const Real pres) {
  // Using a linear function dependent on 
  //   - col, the horizontal column,
  //   - vec, the vector dimension, and
  //   - pres, the current pressure
  // Should ensure that the interpolated values match exactly, since vertical interp is also a linear interpolator.
  // Note, we don't use the level, because here the vertical interpolation is over pressure, so it represents the level.
  return col*pres + vec*100.0;
}

// Helper function to create a remap file
void create_remap_file(const std::string& filename, const int nlevs, const std::vector<std::int64_t>& dofs_p, const std::vector<Real>& p_tgt) 
{

  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::register_dimension(filename,"nlevs","nlevs",nlevs, false);

  scorpio::register_variable(filename,"p_levs","p_levs","none",{"nlevs"},"real","real","Real-nlevs");

  scorpio::set_dof(filename,"p_levs",dofs_p.size(),dofs_p.data()); 
  
  scorpio::eam_pio_enddef(filename);

  scorpio::grid_write_data_array(filename,"p_levs",p_tgt.data(),nlevs);

  scorpio::eam_pio_closefile(filename);
}

TEST_CASE ("vertical_remap") {

  // -------------------------------------- //
  //           Init MPI and PIO             //
  // -------------------------------------- //

  ekat::Comm comm(MPI_COMM_WORLD);

  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());
  scorpio::eam_init_pio_subsystem(fcomm);

  // -------------------------------------- //
  //           Set grid/map sizes           //
  // -------------------------------------- //

  const int nlevs_src  = 2*SCREAM_PACK_SIZE + 2;  // Make sure we check what happens when the vertical extent is a little larger than the max PACK SIZE
  const int nlevs_tgt  = nlevs_src/2;
  const int nldofs_src = 10;
  const int nldofs_tgt =  5;
  const int ngdofs_src = nldofs_src*comm.size();
  const int ngdofs_tgt = nldofs_tgt*comm.size();
  const int nnz_local  = nldofs_src;
  const int nnz        = nnz_local*comm.size();

  // -------------------------------------- //
  //           Create a map file            //
  // -------------------------------------- //

  print (" -> creating map file ...\n",comm);


  std::string filename = "vertical_map_file_np" + std::to_string(comm.size()) + ".nc";

  // Create target pressure levels to be remapped onto
  const Real ptop_tgt = 1;
  const Real pbot_tgt = nlevs_src-2;
  const Real dp_tgt   = (pbot_tgt-ptop_tgt)/(nlevs_tgt-1);
  std::vector<std::int64_t> dofs_p(nlevs_tgt);
  std::iota(dofs_p.begin(),dofs_p.end(),0);
  std::vector<Real> p_tgt;
  for (int ii=0; ii<nlevs_tgt; ++ii) {
    p_tgt.push_back(ptop_tgt + dp_tgt*ii);
  }

  create_remap_file(filename, nlevs_tgt, dofs_p, p_tgt);
  print (" -> creating map file ... done!\n",comm);

  // -------------------------------------- //
  //      Build src grid and remapper       //
  // -------------------------------------- //

  print (" -> creating grid and remapper ...\n",comm);

  const Real mask_val = -99999.0;

  auto src_grid = build_src_grid(comm, nldofs_src, nlevs_src);

  // We need the source pressure level fields for both p_mid and p_int
  auto pmid_src   = create_field("p_mid",  src_grid, false, false, true,  SCREAM_PACK_SIZE);
  auto pint_src   = create_field("p_int",  src_grid, false, false, false, SCREAM_PACK_SIZE);
  // Set the source pressures
  {
    // By adding 1 to the pbot_tgt and subtrating 1 from ptop_tgt we ensure some masking, which 
    // we also want to check.
    const Real ptop_src = ptop_tgt-1;
    const Real pbot_src = pbot_tgt+1;
    const Real dp_src = (pbot_src-ptop_src)/(nlevs_src-1);
    auto pmid_v = pmid_src.get_view<Real**,Host>();
    auto pint_v = pint_src.get_view<Real**,Host>();
    for (int ii=0; ii<pmid_v.extent_int(0); ++ii) {
      pint_v(ii,0) = ptop_src;
      for (int kk=0; kk<nlevs_src; ++kk) {
        pint_v(ii,kk+1) = pint_v(ii,kk) + dp_src;
        pmid_v(ii,kk)   = 0.5*(pint_v(ii,kk) + pint_v(ii,kk+1));
      }
    }
  }
  pmid_src.sync_to_dev();
  pint_src.sync_to_dev();
  auto remap = std::make_shared<VerticalRemapperTester>(src_grid,filename,pmid_src,pint_src,mask_val);
  print (" -> creating grid and remapper ... done!\n",comm);

  // -------------------------------------- //
  //      Create src/tgt grid fields        //
  // -------------------------------------- //

  print (" -> creating fields ...\n",comm);
  constexpr int vec_dim = 3;

  auto tgt_grid = remap->get_tgt_grid();
  // Check that the target grid made by the remapper has the same number of columns as the source grid.
  // Also check that the number of levels matches the expectation.
  REQUIRE(tgt_grid->get_num_vertical_levels()==nlevs_tgt);
  REQUIRE(tgt_grid->get_num_global_dofs()==src_grid->get_num_global_dofs());

  auto src_s2d   = create_field("s2d",  src_grid,true,false);
  auto src_v2d   = create_field("v2d",  src_grid,true,true);
  // For now we can only support PS = SCREAM_PACK_SIZE because the source and target pressure levels will assume that packsize.
  // If we use a smaller packsize we throw an error in the property check step of the vertical interpolation scheme.
  // TODO: Fix that.
  auto src_s3d_m = create_field("s3d_m",src_grid,false,false,true, 1);
  auto src_s3d_i = create_field("s3d_i",src_grid,false,false,false,SCREAM_PACK_SIZE);
  auto src_v3d_m = create_field("v3d_m",src_grid,false,true ,true, 1);
  auto src_v3d_i = create_field("v3d_i",src_grid,false,true ,false,SCREAM_PACK_SIZE);

  auto tgt_s2d   = create_field("s2d",  tgt_grid,true,false);
  auto tgt_v2d   = create_field("v2d",  tgt_grid,true,true);
  auto tgt_s3d_m = create_field("s3d_m",tgt_grid,false,false,true, 1);
  auto tgt_s3d_i = create_field("s3d_i",tgt_grid,false,false,true, SCREAM_PACK_SIZE);
  auto tgt_v3d_m = create_field("v3d_m",tgt_grid,false,true ,true, 1);
  auto tgt_v3d_i = create_field("v3d_i",tgt_grid,false,true ,true, SCREAM_PACK_SIZE);

  std::vector<Field> src_f = {src_s2d,src_v2d,src_s3d_m,src_s3d_i,src_v3d_m,src_v3d_i};
  std::vector<Field> tgt_f = {tgt_s2d,tgt_v2d,tgt_s3d_m,tgt_s3d_i,tgt_v3d_m,tgt_v3d_i};

  const int nfields = src_f.size();

  print (" -> creating fields ... done!\n",comm);

  // -------------------------------------- //
  //     Register fields in the remapper    //
  // -------------------------------------- //

  print (" -> registering fields ...\n",comm);
  remap->registration_begins();
  remap->register_field(src_s2d,  tgt_s2d);
  remap->register_field(src_v2d,  tgt_v2d);
  remap->register_field(src_s3d_m,tgt_s3d_m);
  remap->register_field(src_s3d_i,tgt_s3d_i);
  remap->register_field(src_v3d_m,tgt_v3d_m);
  remap->register_field(src_v3d_i,tgt_v3d_i);
  remap->registration_ends();
  print (" -> registering fields ... done!\n",comm);

  // -------------------------------------- //
  //        Check remapper internals        //
  // -------------------------------------- //

  print (" -> Checking remapper internal state ...\n",comm);


  // Check the pressure levels read from map file

  // -------------------------------------- //
  //       Generate data for src fields     //
  // -------------------------------------- //

  print (" -> generate src fields data ...\n",comm);
  using namespace ShortFieldTagsNames;
  // Generate data in a deterministic way, so that when we check results,
  // we know a priori what the input data that generated the tgt field's
  // values was, even if that data was off rank.
  auto pmid_v = pmid_src.get_view<Real**,Host>();
  auto pint_v = pint_src.get_view<Real**,Host>();
  auto src_gids    = remap->get_src_grid()->get_dofs_gids().get_view<const gid_t*,Host>();
  for (const auto& f : src_f) {
    const auto& l = f.get_header().get_identifier().get_layout();
    switch (get_layout_type(l.tags())) {
      case LayoutType::Scalar2D:
      {
        const auto v_src = f.get_view<Real*,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          v_src(i) = data_func(i,0,pmid_v(i,0));
        }
      } break;
      case LayoutType::Vector2D:
      {
        const auto v_src = f.get_view<Real**,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<vec_dim; ++j) {
            v_src(i,j) = data_func(i,j+1,pmid_v(i,0));
        }}
      } break;
      case LayoutType::Scalar3D:
      {
        const int  nlevs = l.dims().back();
        const auto p_v   = l.has_tag(LEV) ? pmid_v : pint_v;
        const auto v_src = f.get_view<Real**,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<nlevs; ++j) {
            v_src(i,j) = data_func(i,0,p_v(i,j));
        }}
      } break;
      case LayoutType::Vector3D:
      {
        const int nlevs  = l.dims().back();
        const auto p_v   = l.has_tag(LEV) ? pmid_v : pint_v;
        const auto v_src = f.get_view<Real***,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<vec_dim; ++j) {
            for (int k=0; k<nlevs; ++k) {
              v_src(i,j,k) = data_func(i,j+1,p_v(i,k));
        }}}
      } break;
      default:
        EKAT_ERROR_MSG ("Unexpected layout.\n");
    }
    f.sync_to_dev();
  }
  print (" -> generate src fields data ... done!\n",comm);

  // No bwd remap
  REQUIRE_THROWS(remap->remap(false));

  for (int irun=0; irun<5; ++irun) {
    print (" -> run remap ...\n",comm);
    remap->remap(true);
    print (" -> run remap ... done!\n",comm);

    // -------------------------------------- //
    //          Check remapped fields         //
    // -------------------------------------- //

    print (" -> check tgt fields ...\n",comm);
    const auto tgt_gids = tgt_grid->get_dofs_gids().get_view<const gid_t*,Host>();
    const int ntgt_gids = tgt_gids.size();
    for (size_t ifield=0; ifield<tgt_f.size(); ++ifield) {
      const auto& f     = tgt_f[ifield];
      const auto& fsrc  = src_f[ifield];
      const auto& lsrc  = fsrc.get_header().get_identifier().get_layout();
      const auto p_v    = lsrc.has_tag(LEV) ? pmid_v : pint_v;
      const int nlevs_p = lsrc.has_tag(LEV) ? nlevs_src : nlevs_src+1;
      const auto ls     = to_string(lsrc);
      std::string dots (25-ls.size(),'.');
      print ("   -> Checking field with source layout " + to_string(lsrc) +" " + dots + "\n",comm);

      f.sync_to_host();

      switch (get_layout_type(lsrc.tags())) {
        case LayoutType::Scalar2D:
        {
          // This is a flat array w/ no LEV tag so the interpolated value for source and target should match.
          const auto v_src = fsrc.get_view<const Real*,Host>();
          const auto v_tgt = f.get_view<const Real*,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            REQUIRE ( v_tgt(i) == v_src(i) );
          }
        } break;
        case LayoutType::Vector2D:
        {
          // This is a flat array w/ no LEV tag so the interpolated value for source and target should match.
          const auto v_tgt = f.get_view<const Real**,Host>();
          const auto v_src = fsrc.get_view<const Real**,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            const auto gid = tgt_gids(i);
            for (int j=0; j<vec_dim; ++j) {
              REQUIRE ( v_tgt(i,j) == v_src(i,j) );
          }}
        } break;
        case LayoutType::Scalar3D:
        {
          const auto v_tgt = f.get_view<const Real**,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            const auto gid = tgt_gids(i);
            for (int j=0; j<nlevs_tgt; ++j) {
              if (p_tgt[j]>p_v(i,nlevs_p-1) || p_tgt[j]<p_v(i,0)) {
                REQUIRE ( v_tgt(i,j) == mask_val );
              } else {
                REQUIRE ( v_tgt(i,j) == data_func(i,0,p_tgt[j]) );
              }
          }}
        } break;
        case LayoutType::Vector3D:
        {
          const auto v_tgt = f.get_view<const Real***,Host>();
          for (int i=0; i<ntgt_gids; ++i) {
            const auto gid = tgt_gids(i);
            for (int j=0; j<vec_dim; ++j) {
              for (int k=0; k<nlevs_tgt; ++k) {
                const auto term1 = gid*vec_dim*nlevs_tgt+j*nlevs_tgt+k;
                const auto term2 = (gid+ngdofs_tgt)*vec_dim*nlevs_tgt+j*nlevs_tgt+k;
                if (p_tgt[k]>p_v(i,nlevs_p-1) || p_tgt[k]<p_v(i,0)) {
                  REQUIRE ( v_tgt(i,j,k) == mask_val );
                } else {
                  REQUIRE ( v_tgt(i,j,k)== data_func(i,j+1,p_tgt[k]) );
                }
          }}}
        } break;
        default:
          EKAT_ERROR_MSG ("Unexpected layout.\n");
      }

      print ("   -> Checking field with source layout " + to_string(lsrc) + " " + dots + " OK!\n",comm);
    }
    print ("check tgt fields ... done!\n",comm);
  }

  // Clean up scorpio stuff
  scorpio::eam_pio_finalize();
}

} // namespace scream
