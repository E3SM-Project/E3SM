/**
 * @file emulator.cpp
 * @brief Implementation of the Emulator base class.
 */

#include "moab/Types.hpp"
#include <coupler_types.hpp>
#include <emulator.hpp>
#include <emulator_config.hpp>
#include <field_registry.hpp>

// #include <memory_resource>
#include <limits>
#include <moab/iMOAB.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <ekat_yaml.hpp>
#include <iostream>
#include <process.hpp>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace emulator {

RegisteredField
Emulator::make_registered_field(const EmulatorFieldDesc& field) {
  const auto role = field.role == ModelFieldRole::Input ? FieldRole::Import
                                                        : FieldRole::Export;
  return RegisteredField{
      .role = role,
      .component = m_name,
      .attributes = field.config.attributes,
      .size = field.size,
  };
}

std::string to_string(const DimensionDesc& desc) {
  std::ostringstream out;
  out << "Dimension: " << desc.name
      << " Local/Global sizes: " << desc.local_extent << ", "
      << desc.global_extent << "\n";
  return out.str();
}

const DimensionDesc& Emulator::dimension(std::string_view name) const {
  const auto it =
      std::ranges::find(m_grid.dimensions, name, &DimensionDesc::name);
  if (it == m_grid.dimensions.end()) {
    throw std::runtime_error("Requested invalid dimension name: " +
                             std::string(name));
  }
  return *it;
}

const EmulatorFieldDesc&
Emulator::get_field_descriptor(std::string_view name) const {
  const auto it =
      std::ranges::find(field_descriptors_, name,
                        [](const EmulatorFieldDesc& field) -> std::string_view {
                          return field.config.name;
                        });
  if (it == field_descriptors_.end()) {
    throw std::runtime_error("Requested field not found: " + std::string(name));
  }
  return *it;
}

std::span<const EmulatorFieldDesc> Emulator::field_descriptors() const {
  return field_descriptors_;
}

std::span<const double> Emulator::get_field(std::string_view name) {
  const auto& desc = get_field_descriptor(name);
  return field_view(desc);
}

std::size_t Emulator::local_field_size(const FieldConfig& field) const {
  std::size_t size = 1;

  if (field.dimensions.empty()) {
    throw std::runtime_error("Field '" + field.name +
                             "' does not declare any dimensions");
  }

  for (const auto& dim_name : field.dimensions) {
    const auto& desc = dimension(dim_name);
    if (desc.local_extent == 0) {
      throw std::runtime_error(
          "Dimension '" + dim_name +
          "' has zero local extent for field: " + field.name);
    }
    // check that size*desc.local_exten < maximum value of size_t
    if (size > std::numeric_limits<std::size_t>::max() / desc.local_extent) {
      throw std::runtime_error("Inappropriate dimensions for field:\n" +
                               to_string(field) + to_string(desc));
    }
    size *= desc.local_extent;
  }
  return size;
}

std::size_t Emulator::append_fields(std::span<const FieldConfig> fields,
                                    ModelFieldRole role) {

  std::size_t next_offset = 0;

  for (const auto& field : fields) {
    const auto field_size = local_field_size(field);

    if (next_offset > std::numeric_limits<std::size_t>::max() - field_size) {
      throw std::overflow_error("Emulator field-buffer size overflow");
    }

    field_descriptors_.push_back(EmulatorFieldDesc{
        .role = role,
        .config = field,
        .offset = next_offset,
        .size = field_size,
    });

    next_offset += field_size;
  }

  // this should be total size needed for the flat buffer
  return next_offset;
}

std::span<double> Emulator::field_view(const EmulatorFieldDesc& field) {
  auto& buffer =
      field.role == ModelFieldRole::Input ? input_fields_ : output_fields_;

  return std::span<double>{buffer.data() + field.offset, field.size};
}

std::span<double> Emulator::field_view(std::string_view name) {
  const auto& desc = get_field_descriptor(name);
  return field_view(desc);
}

Emulator::Emulator(const CreateConfig cfg)
    : m_config_path(cfg.path), m_grid(cfg.grid), m_component_id(-1),
      m_moab_app_id(-1), m_initialized(false), m_step_count(0) {
  initialize();
}

void Emulator::initialize() {
  if (m_initialized) {
    throw std::runtime_error("Emulator already initialized");
  }

  const auto config = read_emulator_config(m_config_path);

  m_name = config.name;
  if (config.type == "process") {
    m_type = EmulatorType::Process;
  } else {
    m_type = EmulatorType::Component;
  }

  const auto input_size = append_fields(config.inputs, ModelFieldRole::Input);
  const auto output_size =
      append_fields(config.outputs, ModelFieldRole::Output);

  input_fields_.resize(input_size);
  output_fields_.resize(output_size);

  m_initialized = true;
}

void Emulator::populate_registry(FieldRegistry& registry) {

  if (!m_initialized) {
    throw std::logic_error(
        "Cannot populate registry before emulator initialization");
  }

  for (const auto& field : field_descriptors_) {
    m_coupling.register_coupled_field(registry, make_registered_field(field),
                                      field_view(field));
  }
}

void Emulator::run(const ProcessRunOpts opts) {
  if (!m_initialized) {
    throw std::runtime_error("Emulator::run() called before initialize()");
  }
  run_impl(opts);
  m_step_count++;
}

void Emulator::finalize() {
  if (!m_initialized) {
    return; // Already finalized or never initialized
  }
  m_initialized = false;
}

static std::string to_string(EmulatorType t) {
  switch (t) {
  case EmulatorType::ATM:
    return "ATM";
  case EmulatorType::OCN:
    return "OCN";
  case EmulatorType::ICE:
    return "ICE";
  case EmulatorType::LND:
    return "LND";
  default:
    return "UNKNOWN";
  }
}

void Emulator::print_info(std::ostream& os) const {
  os << "Emulator '" << m_name << "'\n";
  os << "  type          : " << to_string(m_type) << "\n";
  os << "  id            : " << m_component_id << "\n";
  os << "  initialized   : " << std::boolalpha << m_initialized << "\n";
  os << "  step_count    : " << m_step_count << "\n";

  // Grid summary via virtual getters
  int nx = get_nx();
  int ny = get_ny();
  int nloc = get_num_local_cols();
  int nglob = get_num_global_cols();

  os << "  grid          : nx=" << nx << " ny=" << ny
     << " num_local_cols=" << nloc << " num_global_cols=" << nglob << "\n";

  // Hook for derived classes to print config / component-specific info
  print_extra_info(os);
}

namespace { // stuff for registering with MOAB

// TODO: replace error checking with e.g. EKAT_REQUIRE (throws exception)

void create_global_id_tag(const Emulator& e) {
  int num_local_cols = e.get_num_local_cols();

  int type = 0; // dense, integer
  int stride = 1;
  int index;

  int app_id = e.moab_app_id();
  ErrCode err =
      iMOAB_DefineTagStorage(&app_id, "GLOBAL_ID", &type, &stride, &index);
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << e.name()
              << " model\n";
    MPI_Abort(e.comm(), err);
  }

  auto temp = e.get_local_col_gids();
  auto gids = std::vector<int>{temp.begin(), temp.end()};
  int ent_type = 0; // entity type (vertex)
  err = iMOAB_SetIntTagStorage(&app_id, "GLOBAL_ID", &num_local_cols, &ent_type,
                               gids.data());
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << e.name()
              << " model\n";
    MPI_Abort(e.comm(), err);
  }

  err = iMOAB_ResolveSharedEntities(&app_id, &num_local_cols, gids.data());
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << e.name()
              << " model\n";
    MPI_Abort(e.comm(), err);
  }
  err = iMOAB_UpdateMeshInfo(&app_id);
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << e.name()
              << " model\n";
    MPI_Abort(e.comm(), err);
  }
}

void create_seq_flds_dom_fields_tag(const Emulator& e) {
  /* FIXME: currently, components create tags from seq_flds_dom_fields, which we
  don't have! :-( int type = 1; // dense, double int size = 1; int app_id =
  e.moab_app_id(); err = iMOAB_DefineTagStorage(&app_id, "seq_flds_dom_fields",
  &type, &size, tag_index); if (err) { fprintf(stderr, "Error: failed to create
  tags from seq_flds_dom_fields\n"); MPI_Abort(e.comm(), err);
  }
  */
}

//----------------------------- pardon the mess!
//-------------------------------------
template <typename T>
ErrCode set_tag_storage(const Emulator& e, const std::string& name,
                        int ent_type, std::span<T> data) {
  return moab::ErrorCode::MB_SUCCESS;
}
template <>
ErrCode set_tag_storage<int>(const Emulator& e, const std::string& name,
                             int ent_type, std::span<int> data) {
  int len = data.size();
  int app_id = e.moab_app_id();
  return iMOAB_SetIntTagStorage(&app_id, name.c_str(), &len, &ent_type,
                                data.data());
}
template <>
ErrCode set_tag_storage<double>(const Emulator& e, const std::string& name,
                                int ent_type, std::span<double> data) {
  int len = data.size();
  int app_id = e.moab_app_id();
  return iMOAB_SetDoubleTagStorage(&app_id, name.c_str(), &len, &ent_type,
                                   data.data());
}
//----------------------------- ////////////////
//-------------------------------------

template <typename T>
void create_tag(const Emulator& e, const std::string& name, int type,
                std::span<T>& data) {
  ErrCode err = set_tag_storage(e, name, type, data); // ^^^^
  if (err) {
    fprintf(stderr, "Error: failed to create %s tag\n", name.c_str());
    MPI_Abort(e.comm(), err);
  }
}

void create_seq_flds_a2x_fields(const Emulator& e) {
  int type = 1;   // dense, double
  int stride = 1; // one value per vertex / entity
  int index;
  int app_id = e.moab_app_id();
  // FIXME: remove quotes around name when we have this data
  ErrCode err = iMOAB_DefineTagStorage(&app_id, "seq_flds_a2x_fields", &type,
                                       &stride, &index);
  if (err) {
    fprintf(stderr, "Error: failed to define seq_flds_a2x_fields\n");
    MPI_Abort(e.comm(), err);
  }

  // make sure this is defined too; it could have the same fields, but in
  // different order, or really different fields; need to make sure we have them
  // FIXME: remove quotes around name when we have this data
  err = iMOAB_DefineTagStorage(&app_id, "seq_flds_x2a_fields", &type, &stride,
                               &index);
  if (err) {
    fprintf(stderr, "Error: failed to define seq_flds_x2a_fields\n");
    MPI_Abort(e.comm(), err);
  }
}

void create_moab_vertices(const Emulator& e, std::span<const double> lat,
                          std::span<const double> lon) {

  int num_local_cols = e.get_num_local_cols();
  std::vector<double> moab_vertex_coords(3 * num_local_cols);
  for (size_t i = 0; i < num_local_cols; ++i) {
    double lat_v = lat[i] * M_PI / 180.0;
    double lon_v = lon[i] * M_PI / 180.0;
    moab_vertex_coords[3 * i] = std::cos(lat_v) * std::cos(lon_v);
    moab_vertex_coords[3 * i + 1] = std::cos(lat_v) * std::sin(lon_v);
    moab_vertex_coords[3 * i + 2] = std::sin(lat_v);
  }

  int dimension = 3;
  int num_coords = dimension * num_local_cols;
  int app_id = e.moab_app_id();
  ;
  ErrCode err = iMOAB_CreateVertices(&app_id, &num_coords, &dimension,
                                     moab_vertex_coords.data());
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << e.name()
              << " model\n";
    MPI_Abort(e.comm(), err);
  }
}

#ifdef MOABDEBUG
void write_mesh(const Emulator& e) {
  std::string filename = e.name() + ".h5m";
  std::string options = "PARALLEL=WRITE_PART";
  ErrCode err =
      iMOAB_WriteMesh(e.moab_app_id(), filename.c_str(), options.c_str());
  if (err) {
    fpritnf(stderr, "Error: fail to write mesh file %s\n", filename.c_str());
    MPI_Abort(e.comm(), err);
  }
}
#endif

} // anonymous namespace

void Emulator::register_with_moab() {
  ErrCode err;

  err = iMOAB_RegisterApplication(m_name.c_str(), &m_grid.comm, &m_component_id,
                                  &m_moab_app_id);
  if (err) {
    MPI_Abort(comm(), err);
  }
  int my_rank;
  MPI_Comm_rank(comm(), &my_rank);
  if (my_rank == 0) {
    printf("Registered MOAB app %d (component ID %d).\n", m_moab_app_id,
           m_component_id);
  }

  int num_local_cols = get_num_local_cols();
  // std::vector<double> data1(num_local_cols), data2(num_local_cols);
  auto data1 = get_latitudes();
  auto data2 = get_longitudes();

  create_global_id_tag(*this); // "GLOBAL_ID"
  create_seq_flds_dom_fields_tag(*this);

  create_tag(*this, "lat", 0, data1);
  create_tag(*this, "lon", 0, data2);

  // pause here to create MOAB's vertices
  create_moab_vertices(*this, data1, data2);

  auto areas = get_areas();
  create_tag(*this, "area", 0, areas);

  // Mask and frac are both exactly 1
  // fill(data1.begin(), data1.end(), 1.0);
  // create_tag(*this, "mask", 0, data1);
  // create_tag(*this, "frac", 0, data1);

  //! call mct_gGrid_importRAttr(dom_atm,"mask",data1,lsize)
  //! call mct_gGrid_importRAttr(dom_atm,"frac",data1,lsize)

  // Aream is computed by mct, so give invalid initial value
  // fill(data1.begin(), data1.end(), -9999.0);
  create_tag(*this, "aream", 0, data1);

#ifdef MOABDEBUG
  write_mesh(e);
#endif

  // FIXME: uncomment this line when we have this data
  // create_seq_flds_a2x_fields(e);
}

} // namespace emulator
