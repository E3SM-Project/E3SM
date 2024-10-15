#include "share/io/scream_output_manager.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"

namespace scream
{

constexpr int ngcols_data  = 12;
constexpr int nlevs_data   = 20;
constexpr int nsteps_data  = 5;
constexpr int dt_data      = 100;
constexpr int nlevs_filled = 2;
constexpr double fill_val  = 1e30;

util::TimeStamp get_t0 () {
  return  util::TimeStamp ({2000,1,1},{12,0,0});
}

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ngcols, const int nlevs) {

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_t{"Point Grid"});
  auto& pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", ngcols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

std::shared_ptr<FieldManager>
create_fm (const std::shared_ptr<const AbstractGrid>& grid)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using FR = FieldRequest;

  auto fm = std::make_shared<FieldManager>(grid);

  const std::string& gn = grid->name();

  auto scalar3d = grid->get_3d_scalar_layout(true);
  auto vector3d = grid->get_3d_vector_layout(true,2);

  FieldIdentifier fid1("p_mid",scalar3d,Pa,gn);
  FieldIdentifier fid2("horiz_winds",vector3d,m/s,gn);

  // Register fields with fm
  fm->registration_begins();
  fm->register_field(FR(fid1));
  fm->register_field(FR(fid2));
  fm->registration_ends();

  auto U = fm->get_field("horiz_winds").subfield("U",1,0);
  auto V = fm->get_field("horiz_winds").subfield("V",1,1);
  fm->add_field(U);
  fm->add_field(V);

  return fm;
}

void compute_field (Field f,
                    const util::TimeStamp& time,
                    const ekat::Comm& comm,
                    const int num_masked_levs = 0)
{
  const auto& fl = f.get_header().get_identifier().get_layout();

  const int ncols = fl.dim(0);
  const int nlevs = fl.dim(1);

  int offset = ncols;
  comm.scan(&offset,1,MPI_SUM);
  offset -= ncols;

  const int step = time.get_num_steps();
  const int lev_beg = num_masked_levs;
  const int lev_end = nlevs - num_masked_levs;

  if (fl.rank()==2) {
    const auto f_h = f.get_view<Real**, Host>();
    for (int icol=0;icol<ncols;++icol) {
      for (int ilev=0; ilev<lev_beg; ++ilev) {
        f_h(icol,ilev) = fill_val;
      }
      for (int ilev=lev_beg;ilev<lev_end;++ilev) {
        f_h(icol,ilev) = step + (offset+icol)*nlevs + ilev + 1;
      }
      for (int ilev=lev_end; ilev<nlevs; ++ilev) {
        f_h(icol,ilev) = fill_val;
      }
    }
  } else {
    const auto f_h = f.get_view<Real***, Host>();
    for (int icol=0;icol<ncols;++icol) {
      for (int ilev=0; ilev<lev_beg; ++ilev) {
        f_h(icol,0,ilev) = f_h(icol,1,ilev) = fill_val;
      }
      for (int ilev=lev_beg;ilev<lev_end;++ilev) {
        f_h(icol,0,ilev) = step + (offset+icol)*2*nlevs + ilev + 1;
        f_h(icol,1,ilev) = step + (offset+icol)*2*nlevs + ilev + 1;
      }
      for (int ilev=lev_end; ilev<nlevs; ++ilev) {
        f_h(icol,0,ilev) = f_h(icol,1,ilev) = fill_val;
      }
    }
  }
  f.sync_to_dev();
  f.get_header().get_tracking().update_time_stamp(time);
}

void compute_fields (const std::shared_ptr<FieldManager>& fm,
                     const util::TimeStamp& time,
                     const ekat::Comm& comm,
                     const int num_masked_levs = 0,
                     const bool update_p_mid = true)
{
  if (update_p_mid) {
    // Don't mask pressure
    compute_field(fm->get_field("p_mid"),time,comm,0);
  }
  compute_field(fm->get_field("U"),time,comm,num_masked_levs);
  compute_field(fm->get_field("V"),time,comm,num_masked_levs);

  // Not sure if we need it, since we don't handle horiz_winds directly, I think
  fm->get_field("horiz_winds").get_header().get_tracking().update_time_stamp(time);
}

std::shared_ptr<OutputManager>
create_om (const std::string& filename_prefix,
           const std::shared_ptr<FieldManager>& fm,
           const std::shared_ptr<GridsManager>& gm,
           const util::TimeStamp& t0,
           const ekat::Comm& comm)
{
  using strvec_t = std::vector<std::string>;

  // NOTE: ask "real" fp precision, so even when building in double precision
  //       we can retrieve exactly the nudging data (if no remapping happens)
  ekat::ParameterList params;
  params.set<std::string>("Averaging Type","INSTANT");
  params.set<std::string>("filename_prefix",filename_prefix);
  params.set<std::string>("Floating Point Precision","real");
  params.set("Field Names",strvec_t{"p_mid","U","V"});
  params.set("fill_value",fill_val);

  auto& ctrl_pl = params.sublist("output_control");
  ctrl_pl.set<std::string>("frequency_units","nsteps");
  ctrl_pl.set("Frequency",1);
  ctrl_pl.set("save_grid_data",false);

  auto om = std::make_shared<OutputManager>();
  om->initialize(comm,params,t0,false);
  om->setup(fm,gm);
  return om;
}

} // namespace scream
