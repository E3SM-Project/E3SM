#include "share/io/scorpio_scm_input.hpp"
#include "share/io/scorpio_input.hpp"

#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"

#include <ekat/util/ekat_string_utils.hpp>

#include <memory>

// Anonymous namespace, in case other TU's define the same type
namespace {
struct RealInt
{
  scream::Real val;
  int  idx;
};
}
// Specialize ekat's get_mpi_type for RealInt (needed for the all_reduce call)
namespace ekat {
  template<>
  MPI_Datatype get_mpi_type<RealInt> () {
#ifdef SCREAM_DOUBLE_PRECISION
    return MPI_DOUBLE_INT;
#else
    return MPI_FLOAT_INT;
#endif
  }
} // namespace ekat

namespace scream
{

SCMInput::
SCMInput (const std::string& filename,
          const double lat, const double lon,
          const std::vector<Field>& fields,
          const ekat::Comm& comm)
 : m_comm(comm)
 , m_filename (filename)
{
  auto iotype = scorpio::str2iotype("default");
  scorpio::register_file(m_filename,scorpio::Read,iotype);

  // Some input files have the "time" dimension as non-unlimited. This messes up our
  // scorpio interface, which stores a pointer to a "time" dim to be used to read/write
  // slices. This ptr is automatically inited to the unlimited dim in the file. If there is
  // no unlim dim, this ptr remains inited.
  if (not scorpio::has_time_dim(m_filename) and scorpio::has_dim(m_filename,"time")) {
    scorpio::mark_dim_as_time(m_filename,"time");
  }

  create_io_grid ();
  auto ncols = m_io_grid->get_num_local_dofs();

  create_closest_col_info (lat, lon);

  // Init fields specs
  for (const auto& f : fields) {
    const auto& fh  = f.get_header();
    const auto& fid = fh.get_identifier();
    const auto& fl  = fid.get_layout();

    EKAT_REQUIRE_MSG (fl.tags()[0]==FieldTag::Column,
      "Error! SCMInput only works for physics-type layouts.\n"
      "  - field name: " + f.name() + "\n"
      "  - field layout: " + fl.to_string() + "\n");

    m_fields.push_back(f.subfield(0,0));
    FieldIdentifier fid_io(f.name(),fl.clone().reset_dim(0,ncols),fid.get_units(),m_io_grid->name());
    auto& f_io = m_io_fields.emplace_back(fid_io);
    f_io.allocate_view();
  }

  // Init scorpio internal structures
  init_scorpio_structures ();
}

SCMInput::
~SCMInput ()
{
  scorpio::release_file(m_filename);
}

void SCMInput::create_io_grid ()
{
  EKAT_REQUIRE_MSG (scorpio::has_dim(m_filename,"ncol"),
      "Error! Dimension 'ncol' not found in input file.\n"
      "  - filename: " + m_filename + "\n");
  const int ncols = scorpio::get_dimlen(m_filename,"ncol");
  const int nlevs = scorpio::has_dim(m_filename,"lev") ? scorpio::get_dimlen(m_filename,"lev") : 1;

  m_io_grid = create_point_grid("scm_io_grid",ncols,nlevs,m_comm);
}

void SCMInput::create_closest_col_info (double target_lat, double target_lon)
{
  // Read lat/lon fields
  const auto ncols = m_io_grid->get_num_local_dofs();

  auto nondim = ekat::units::Units::nondimensional();
  auto lat = m_io_grid->create_geometry_data("lat",m_io_grid->get_2d_scalar_layout(),nondim);
  auto lon = m_io_grid->create_geometry_data("lon",m_io_grid->get_2d_scalar_layout(),nondim);

  // Read from file
  AtmosphereInput file_reader(m_filename, m_io_grid, {lat,lon});
  file_reader.read_variables();
  file_reader.finalize();

  // Find column index of closest lat/lon to target_lat/lon params
  auto lat_d = lat.get_view<Real*>();
  auto lon_d = lon.get_view<Real*>();
  using minloc_t = Kokkos::MinLoc<Real,int>;
  using minloc_value_t = typename minloc_t::value_type;
  minloc_value_t minloc;
  Kokkos::parallel_reduce(ncols, KOKKOS_LAMBDA (int icol, minloc_value_t& result) {
    auto dist = std::abs(lat_d(icol)-target_lat)+std::abs(lon_d(icol)-target_lon);
    if(dist<result.val) {
      result.val = dist;
      result.loc = icol;
    }
  }, minloc_t(minloc));

  // Find processor with closest lat/lon match
  const auto my_rank = m_comm.rank();
  RealInt min_dist_and_rank = {minloc.val, my_rank};
  m_comm.all_reduce(&min_dist_and_rank, 1, MPI_MINLOC);

  // Set local col idx to -1 for mpi ranks not containing minimum lat/lon distance
  m_closest_col_info.mpi_rank = min_dist_and_rank.idx;
  m_closest_col_info.col_lid  = my_rank==min_dist_and_rank.idx ? minloc.loc : -1;
}

void SCMInput::read_variables (const int time_index)
{
  auto func_start = std::chrono::steady_clock::now();
  auto fname = [](const Field& f) { return f.name(); };
  if (m_atm_logger) {
    m_atm_logger->info("[EAMxx::scorpio_scm_input] Reading variables from file");
    m_atm_logger->info("  file name: " + m_filename);
    m_atm_logger->info("  var names: " + ekat::join(m_fields,fname,", "));
    if (time_index!=-1) {
      m_atm_logger->info("  time idx : " + std::to_string(time_index));
    }
  }

  // MPI rank with closest column index store column data
  const int n = m_fields.size();
  for (int i=0; i<n; ++i) {
    auto& f_io = m_io_fields[i];
    const auto& name = f_io.name();

    // Read the data
    scorpio::read_var(m_filename,name,f_io.get_internal_view_data<Real,Host>(),time_index);

    auto& f = m_fields[i];
    if (m_comm.rank() == m_closest_col_info.mpi_rank) {
      // This rank read the column we need, so copy it in the output field
      f.deep_copy<Host>(f_io.subfield(0,m_closest_col_info.col_lid));
    }

    // Broadcast column data to all other ranks
    const auto col_size = f.get_header().get_identifier().get_layout().size();
    m_comm.broadcast(f.get_internal_view_data<Real,Host>(), col_size, m_closest_col_info.mpi_rank);

    // Sync fields to device
    f.sync_to_dev();
  }

  auto func_finish = std::chrono::steady_clock::now();
  if (m_atm_logger) {
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start)/1000.0;
    m_atm_logger->info("  Done! Elapsed time: " + std::to_string(duration.count()) +" seconds");
  }
}

void SCMInput::init_scorpio_structures()
{
  using namespace ShortFieldTagsNames;

  // Check variables are in the input file
  for (const auto& f : m_io_fields) {
    const auto& layout = f.get_header().get_identifier().get_layout();
    auto dim_names = layout.names();
    for (int i=0; i<layout.rank(); ++i) {
      // If t==CMP, and the name stored in the layout is the default ("dim"),
      // we append also the extent, to allow different vector dims in the file
      if (dim_names[i]=="dim") {
        dim_names[i] += std::to_string(layout.dim(i));
      }
    }

    // Check that the variable is in the file.
    EKAT_REQUIRE_MSG (scorpio::has_var(m_filename,f.name()),
        "Error! Input file does not store a required variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + f.name() + "\n");

    const auto& var = scorpio::get_var(m_filename,f.name());
    EKAT_REQUIRE_MSG (var.dim_names()==dim_names,
        "Error! Layout mismatch for input file variable.\n"
        " - filename: " + m_filename + "\n"
        " - varname : " + f.name() + "\n"
        " - expected dims : " + ekat::join(dim_names,",") + "\n"
        " - dims from file: " + ekat::join(var.dim_names(),",") + "\n");

    // Check that all dims for this var match the ones on file
    for (int i=0; i<layout.rank(); ++i) {
      const int file_len  = scorpio::get_dimlen(m_filename,dim_names[i]);
      const int eamxx_len = i==0 ? m_io_grid->get_partitioned_dim_global_size() : layout.dim(i);
      EKAT_REQUIRE_MSG (eamxx_len==file_len,
          "Error! Dimension mismatch for input file variable.\n"
        " - filename : " + m_filename + "\n"
        " - varname  : " + f.name() + "\n"
        " - var dims : " + ekat::join(dim_names,",") + "\n"
        " - dim name : " + dim_names[i] + "\n"
        " - expected extent : " + std::to_string(eamxx_len) + "\n"
        " - extent from file: " + std::to_string(file_len) + "\n");
    }

    // Ensure that we can read the var using Real data type
    scorpio::change_var_dtype (m_filename,f.name(),"real");
  }

  // Set decompositions for the variables
  set_decompositions();
}

/* ---------------------------------------------------------- */
void SCMInput::set_decompositions()
{
  using namespace ShortFieldTagsNames;

  // First, check if any of the vars is indeed partitioned
  const auto decomp_tag  = m_io_grid->get_partitioned_dim_tag();

  bool has_decomposed_layouts = false;
  for (const auto& f : m_io_fields) {
    const auto& layout = f.get_header().get_identifier().get_layout();
    if (layout.has_tag(decomp_tag)) {
      has_decomposed_layouts = true;
      break;
    }
  }
  if (not has_decomposed_layouts) {
    // If none of the input vars are decomposed on this grid,
    // then there's nothing to do here
    return;
  }

  // Set the decomposition for the partitioned dimension
  const int local_dim = m_io_grid->get_partitioned_dim_local_size();
  std::string decomp_dim = m_io_grid->has_special_tag_name(decomp_tag)
                         ? m_io_grid->get_special_tag_name(decomp_tag)
                         : e2str(decomp_tag);

  auto gids_f = m_io_grid->get_partitioned_dim_gids();
  auto gids_h = gids_f.get_view<const AbstractGrid::gid_type*,Host>();
  auto min_gid = m_io_grid->get_global_min_partitioned_dim_gid();
  std::vector<scorpio::offset_t> offsets(local_dim);
  for (int idof=0; idof<local_dim; ++idof) {
    offsets[idof] = gids_h[idof] - min_gid;
  }
  scorpio::set_dim_decomp(m_filename,decomp_dim,offsets);
}

} // namespace scream
