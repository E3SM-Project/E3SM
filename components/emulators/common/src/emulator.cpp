/**
 * @file emulator.cpp
 * @brief Implementation of the Emulator base class.
 */

#include "emulator.hpp"
#include "include/emulator.hpp"

//#include <memory_resource>
#include <moab/iMOAB.h>

#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace std;

namespace emulator {

Emulator::Emulator(EmulatorType type, MPI_Comm comm, int id, const string &name)
    : m_type(type), m_comm(comm), m_component_id(id), m_moab_app_id(-1), m_name(name), m_initialized(false),
      m_step_count(0) {}

void Emulator::initialize() {
  if (m_initialized) {
    throw runtime_error("Emulator already initialized");
  }
  init_impl();
  m_initialized = true;
}

void Emulator::run(int dt) {
  if (!m_initialized) {
    throw runtime_error("Emulator::run() called before initialize()");
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


static string to_string(EmulatorType t) {
  switch (t) {
  case EmulatorType::ATM: return "ATM";
  case EmulatorType::OCN: return "OCN";
  case EmulatorType::ICE: return "ICE";
  case EmulatorType::LND: return "LND";
  default: return "UNKNOWN";
  }
}

void Emulator::print_info(ostream& os) const {
  os << "Emulator '" << m_name << "'\n";
  os << "  type          : " << to_string(m_type) << "\n";
  os << "  id            : " << m_component_id << "\n";
  os << "  initialized   : " << boolalpha << m_initialized << "\n";
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

void create_global_id_tag(const Emulator &e) {
  int num_local_cols = e.get_num_local_cols();

  int type = 0; // dense, integer
  int stride = 1;
  int index;

  int app_id = e.moab_app_id();
  ErrCode err = iMOAB_DefineTagStorage(&app_id, "GLOBAL_ID", &type, &stride, &index);
  if (err) {
    fprintf(stderr, "Error: failed to define GLOBAL_ID tag in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }

  vector<int> gids(num_local_cols);
  e.get_local_col_gids(gids.data());
  int ent_type = 0; // entity type (vertex)
  err = iMOAB_SetIntTagStorage(&app_id, "GLOBAL_ID", &num_local_cols, &ent_type, gids.data());
  if (err) {
    fprintf(stderr, "Error: failed to create GLOBAL_ID tag in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }

  err = iMOAB_ResolveSharedEntities(&app_id, &num_local_cols, gids.data());
  if (err) {
    fprintf(stderr, "Error: failed to resolve shared entities in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }
  err = iMOAB_UpdateMeshInfo(&app_id);
  if (err) {
    fprintf(stderr, "Error: failed to update mesh info in %s model\n", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }
}

void create_seq_flds_dom_fields_tag(const Emulator &e) {
  /* FIXME: currently, components create tags from seq_flds_dom_fields, which we don't have! :-(
  int type = 1; // dense, double
  int size = 1;
  int app_id = e.moab_app_id();
  err = iMOAB_DefineTagStorage(&app_id, "seq_flds_dom_fields", &type, &size, tag_index);
  if (err) {
    fprintf(stderr, "Error: failed to create tags from seq_flds_dom_fields\n");
    MPI_Abort(e.comm(), err);
  }
  */
}

//----------------------------- pardon the mess! -------------------------------------
template <typename T>
ErrCode set_tag_storage(const Emulator &e, const string &name, int ent_type,
                        vector<T> &data) { return 1; }
template <>
ErrCode set_tag_storage<int>(const Emulator &e, const string &name, int ent_type,
                             vector<int> &data) {
  int len = data.size();
  int app_id = e.moab_app_id();
  return iMOAB_SetIntTagStorage(&app_id, name.c_str(), &len, &ent_type, data.data());
}
template <>
ErrCode set_tag_storage<double>(const Emulator &e, const string &name, int ent_type,
                                vector<double> &data) {
  int len = data.size();
  int app_id = e.moab_app_id();
  return iMOAB_SetDoubleTagStorage(&app_id, name.c_str(), &len, &ent_type, data.data());
}
//----------------------------- //////////////// -------------------------------------

template <typename T>
void create_tag(const Emulator &e, const string &name, int type, vector<T> &data) {
  ErrCode err = set_tag_storage(e, name, type, data); // ^^^^
  if (err) {
    fprintf(stderr, "Error: failed to create %s tag\n", name.c_str());
    MPI_Abort(e.comm(), err);
  }
}

void create_seq_flds_a2x_fields(const Emulator &e) {
  int type = 1; // dense, double
  int stride = 1;  // one value per vertex / entity
  int index;
  int app_id = e.moab_app_id();
  // FIXME: remove quotes around name when we have this data
  ErrCode err = iMOAB_DefineTagStorage(&app_id, "seq_flds_a2x_fields", &type, &stride, &index);
  if (err) {
    fprintf(stderr, "Error: failed to define seq_flds_a2x_fields\n");
    MPI_Abort(e.comm(), err);
  }

  // make sure this is defined too; it could have the same fields, but in different order, or
  // really different fields; need to make sure we have them
  // FIXME: remove quotes around name when we have this data
  err = iMOAB_DefineTagStorage(&app_id, "seq_flds_x2a_fields", &type, &stride, &index);
  if (err) {
    fprintf(stderr, "Error: failed to define seq_flds_x2a_fields\n");
    MPI_Abort(e.comm(), err);
  }
}

void create_moab_vertices(const Emulator & e,
                          const vector<double> &lat,
                          const vector<double> &lon) {

  int num_local_cols = e.get_num_local_cols();
  vector<double> moab_vertex_coords(3 * num_local_cols);
  for (size_t i = 0; i < num_local_cols; ++i) {
    double lat_v = lat[i] * M_PI/180.0;
    double lon_v = lon[i] * M_PI/180.0;
    moab_vertex_coords[3*i]   = std::cos(lat_v) * std::cos(lon_v);
    moab_vertex_coords[3*i+1] = std::cos(lat_v) * std::sin(lon_v);
    moab_vertex_coords[3*i+2] = std::sin(lat_v);
  }

  int dimension = 3;
  int num_coords = dimension * num_local_cols;
  int app_id = e.moab_app_id();;
  ErrCode err = iMOAB_CreateVertices(&app_id, &num_coords, &dimension,
                                     moab_vertex_coords.data());
  if (err) {
    fprintf(stderr, "Error! Couldn't create MOAB vertices in %s model", e.name().c_str());
    MPI_Abort(e.comm(), err);
  }
}

#ifdef MOABDEBUG
void write_mesh(const Emulator &e) {
  string filename = e.name() + ".h5m";
  string options = "PARALLEL=WRITE_PART";
  ErrCode err = iMOAB_WriteMesh(e.moab_app_id(), filename.c_str(), options.c_str());
  if (err) {
    fpritnf(stderr, "Error: fail to write mesh file %s\n", filename.c_str());
    MPI_Abort(e.comm(), err);
  }
}
#endif

} // anonymous namespace

void Emulator::register_with_moab() {
  ErrCode err;

  err = iMOAB_RegisterApplication(m_name.c_str(), &m_comm,
                                  &m_component_id, &m_moab_app_id);
  if (err) {
    MPI_Abort(m_comm, err);
  }
  int my_rank;
  MPI_Comm_rank(m_comm, &my_rank);
  if (my_rank == 0) {
    printf("Registered MOAB app %d (component ID %d).\n", m_moab_app_id, m_component_id);
  }

  int num_local_cols = get_num_local_cols();
  vector<double> data1(num_local_cols), data2(num_local_cols);

  create_global_id_tag(*this); // "GLOBAL_ID"
  create_seq_flds_dom_fields_tag(*this);

  // fetch and store latitudes (data1) and longitudes (data2)
  get_cols_latlon(data1.data(), data2.data());

  create_tag(*this, "lat", 0, data1);
  create_tag(*this, "lon", 0, data2);

  // pause here to create MOAB's vertices
  create_moab_vertices(*this, data1, data2);

  get_cols_area(data1.data());
  create_tag(*this, "area", 0, data1);

  // Mask and frac are both exactly 1
  fill(data1.begin(), data1.end(), 1.0);
  create_tag(*this, "mask", 0, data1);
  create_tag(*this, "frac", 0, data1);

  //!call mct_gGrid_importRAttr(dom_atm,"mask",data1,lsize)
  //!call mct_gGrid_importRAttr(dom_atm,"frac",data1,lsize)

  // Aream is computed by mct, so give invalid initial value
  fill(data1.begin(), data1.end(), -9999.0);
  create_tag(*this, "aream", 0, data1);

#ifdef MOABDEBUG
  write_mesh(e);
#endif

  // FIXME: uncomment this line when we have this data
  //create_seq_flds_a2x_fields(e);
}

} // namespace emulator
