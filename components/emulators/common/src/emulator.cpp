/**
 * @file emulator.cpp
 * @brief Implementation of the Emulator base class.
 */

#include "emulator.hpp"
#include "include/emulator.hpp"

#include <memory_resource>
#include <moab/iMOAB.h>

#include <stdexcept>
#include <iostream>

namespace emulator {

Emulator::Emulator(EmulatorType type, MPI_Comm comm, int id, const std::string &name)
    : m_type(type), m_comm(comm), m_component_id(id), m_moab_app_id(-1), m_name(name), m_initialized(false),
      m_step_count(0) {}

void Emulator::initialize() {
  if (m_initialized) {
    throw std::runtime_error("Emulator already initialized");
  }
  init_impl();
  m_initialized = true;
}

void Emulator::run(int dt) {
  if (!m_initialized) {
    throw std::runtime_error("Emulator::run() called before initialize()");
  }
  run_impl(dt);
  m_step_count++;
}

void Emulator::finalize() {
  if (!m_initialized) {
    return; // Already finalized or never initialized
  }
  final_impl();
  m_initialized = false;
}


static std::string to_string(EmulatorType t) {
  switch (t) {
  case EmulatorType::ATM_COMP: return "ATM_COMP";
  case EmulatorType::OCN_COMP: return "OCN_COMP";
  case EmulatorType::ICE_COMP: return "ICE_COMP";
  case EmulatorType::LND_COMP: return "LND_COMP";
  default:                     return "UNKNOWN";
  }
}

void Emulator::print_info(std::ostream& os) const {
  os << "Emulator '" << m_name << "'\n";
  os << "  type          : " << to_string(m_type) << "\n";
  os << "  id            : " << m_id << "\n";
  os << "  initialized   : " << std::boolalpha << m_initialized << "\n";
  os << "  step_count    : " << m_step_count << "\n";

  // Grid summary via virtual getters
  int nx    = get_nx();
  int ny    = get_ny();
  int nloc  = get_num_local_cols();
  int nglob = get_num_global_cols();

  os << "  grid          : nx=" << nx
     << " ny=" << ny
     << " num_local_cols=" << nloc
     << " num_global_cols=" << nglob << "\n";

  // Hook for derived classes to print config / component-specific info
  print_extra_info(os);
}

namespace { // stuff for registering with MOAB

// TODO: replace error checking with e.g. EKAT_REQUIRE (throws exception)
// TODO: abstract away some boilerplate?

void create_global_id_tag(const Emulator &e) {
  int num_local_cols = e.get_num_local_cols();

  int type = 0; // dense, integer
  int stride = 1;
  int index;

  ErrCode err = iMOAB_DefineTagStorage(e.moab_app_id(), "GLOBAL_ID", &type, &stride, &index);
  if (err) {
    fprintf(stderr, "Error: failed to define GLOBAL_ID tag in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }

  std::vector<int> idata(num_local_columns);
  get_local_col_gids(idata.data());
  int ent_type = 0; // entity type (vertex)
  err = iMOAB_SetIntTagStorage(e.moab_app_id(), "GLOBAL_ID", num_local_cols, ent_type, idata);
  if (err) {
    fprintf(stderr, "Error: failed to create GLOBAL_ID tag in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }

  err = iMOAB_ResolveSharedEntities(e.moab_app_id(), num_local_cols, gids);
  if (err) {
    fprintf(stderr, "Error: failed to resolve shared entities in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }
  err = iMOAB_UpdateMeshInfo(e.moab_app_id());
  if (err) {
    fprintf(stderr, "Error: failed to update mesh info in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }
}

void create_seq_flds_dom_fields_tag(const Emulator &e) {
  /* FIXME: currently, components create tags from seq_flds_dom_fields, which we don't have! :-(
  tag_type = 1; // dense, double
  err = iMOAB_DefineTagStorage(e.moab_app_id(), "seq_flds_dom_fields", &tag_type, &tag_size, tag_index);
  if (err) {
    fprintf(stderr, "Error: failed to create tags from seq_flds_dom_fields\n");
    MPI_abort(e.comm(), err);
  }
  */
}

//----------------------------- pardon the mess! -------------------------------------
template <typename T>
ErrCode set_tag_storage(const Emulator &e, const std::string &name, int ent_type,
                       const std::vector<T> &data) {}
template <>
ErrCode set_tag_storage<int>(const Emulator &e, const std::string &name, int ent_type,
                             const std::vector<T> &data) {
  int len = data.size();
  return iMOAB_SetIntTagStorage(e.moab_app_id(), name.c_str(), &len, ent_type, data);
}
template <>
ErrCode set_tag_storage<double>(const Emulator &e, const std::string &name, int ent_type,
                                const std::vector<T> &data) {
  int len = data.size();
  return iMOAB_SetDoubleTagStorage(e.moab_app_id(), name.c_str(), &len, &ent_type, data.data());
}
//----------------------------- //////////////// -------------------------------------

template <typename T>
void create_tag(const Emulator &e, const std::string &name, int type, const std::vector<T> &data) {
  ErrCode err = set_tag_storage(e, name, type, data); // ^^^^
  if (ierr) {
    fprintf(stderr, "Error: failed to create %s tag\n", name.c_str());
    MPI_abort(e.comm(), err);
  }
}

void create_seq_flds_a2x_fields(const Emulator &e) {
  int type = 1; // dense, double
  int stride = 1;  // one value per vertex / entity
  int index;
  // FIXME: remove quotes around name when we have this data
  err = iMOAB_DefineTagStorage(mphaid, "seq_flds_a2x_fields", &type, &stride, &index);
  if (err) {
    fprintf(stderr, "Error: failed to define seq_flds_a2x_fields\n");
    MPI_abort(e.comm(), err);
  }

  // make sure this is defined too; it could have the same fields, but in different order, or
  // really different fields; need to make sure we have them
  // FIXME: remove quotes around name when we have this data
  err = iMOAB_DefineTagStorage(mphaid, "seq_flds_x2a_fields", &type, &stride, &index);
  if (err) {
    fprintf(stderr, "Error: failed to define seq_flds_x2a_fields\n");
    MPI_abort(e.comm(), err);
  }
}

void create_moab_vertices(const Emulator & e,
                          const std::vector<double> &lat,
                          const std::vector<double> &lon) {

  int num_local_cols = e.get_num_local_cols();
  std::vector<double> moab_vertex_coords(3 * num_local_cols);
  for (std::size_t i = 0; i < num_local_cols; ++i) {
    double lat_v = lat[i] * M_PI/180.0;
    double lon_v = lon[i] * M_PI/180.0;
    moab_vertex_coords[3*i]   = math::cos(latv) * math::cos(lonv);
    moab_vertex_coords[3*i+1] = math::cos(latv) * math::sin(lonv);
    moab_vertex_coords[3*i+2] = math::sin(latv);
  }

  int dimension = 3;
  int num_coords = dimension * num_local_cols;
  ErrCode err = iMOAB_CreateVertices(e.moab_app_id(), &num_coords, &dimension,
                                     moab_vertex_coords.data());
  if (err) {
    fprintf(stderr, "Error! Couldn't create MOAB vertices in %s model", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }
}

void Emulator::register_with_moab() {
  ErrCode err;

  err = iMOAB_RegisterApplication(m_name.c_str(), m_comm,
                                  &m_component_id, &m_moab_app_id);
  if (err) {
    MPI_Abort(m_comm, err);
  }
  int my_rank;
  MPI_Comm_rank(m_comm, &my_rank);
  if (my_rank == 0) {
    printf("Registered MOAB app %d (component ID %d).\n", m_moab_app_id, m_component_id);
  }

  ErrCode err;
  int tag_type, tag_stride;
  string tag_name;

  int num_local_cols = e.get_num_local_cols();
  std::vector<double> data1(num_local_cols), data2(num_local_cols);

  create_global_id_tag(e); // "GLOBAL_ID"
  create_seq_flds_dom_fields_tag(e);

  // fetch and store latitudes (data1) and longitudes (data2)
  e.get_cols_latlon(data1.data(), data2.data());

  create_tag(e, "lat", num_local_cols, 0, data1);
  create_tag(e, "lon", num_local_cols, 0, data2);

  // pause here to create MOAB's vertices
  create_moab_vertices(e, data1, data2);

  e.get_cols_area(data1.data());
  create_tag(e, "area", num_local_cols, 0, data1);

  // Mask and frac are both exactly 1
  std::fill(data1.begin(), data1.end(), 1.0);
  create_tag(e, "mask", num_local_cols, 0, data1);
  create_tag(e, "frac", num_local_cols, 0, data1);

  //!call mct_gGrid_importRAttr(dom_atm,"mask",data1,lsize)
  //!call mct_gGrid_importRAttr(dom_atm,"frac",data1,lsize)

  // Aream is computed by mct, so give invalid initial value
  std::fill(data1.begin(), data1.end(), -9999.0);
  create_tag(e, "aream", num_local_cols, 0, data1);

#ifdef MOABDEBUG
    outfile = 'AtmPhys.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mphaid, outfile, wopts)
    if (ierr > 0 )  then
      print *, "Error: fail to write PhysAtm mesh "
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
#endif

  // FIXME: uncomment this line when we have this data
  //create_seq_flds_a2x_fields(e);
}

} // namespace emulator
