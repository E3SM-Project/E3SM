#include "moab_adapter.hpp"
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <limits>

#include <moab/Types.hpp>
#include <moab/iMOAB.h>
#include <stdexcept>
namespace e3sm::coupler {

namespace {
int checked_size_to_int(std::size_t size) {
  if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("Size excedes int length");
  }
  return static_cast<int>(size);
}

void checkMoabErr(ErrCode Err, const std::string& msg, MPI_Comm comm) {
  if (Err != moab::MB_SUCCESS) {
    std::cerr << msg << std::endl;
    MPI_Abort(comm, Err);
  }
}


} // namespace

MoabAdapter::MoabAdapter(std::string component_name, int id,
                         const MoabPartition& partition)
    : name_(std::move(component_name)), component_id_(id),
      comm_(partition.comm), local_vertex_count_(partition.global_ids.size()) {
  register_application(partition);
}

MoabAdapter::~MoabAdapter() {
  if (app_id_ >= 0) {
    ErrCode err = iMOAB_DeregisterApplication(&app_id_);
    checkMoabErr(err, "Failed to deregister MOAB application " + name_ + '\n',
                 comm_);
  }
}

//

void MoabAdapter::define_tag(const MoabTagDescriptor& descriptor) {
  if (descriptor.name.empty() ||
      descriptor.name.find(':') != std::string::npos) {
    throw std::invalid_argument(
        "Expected one nonempty MOAB tag name without ':'");
  }

  if (descriptor.values_per_vertex == 0) {
    throw std::invalid_argument(
        "MOAB tag must have at least one value per entity");
  }

  if (tags_.contains(descriptor.name)) {
    throw std::runtime_error("MOAB tag already registered: " + descriptor.name);
  }

  int type = DENSE_DOUBLE;
  int components = checked_size_to_int(descriptor.values_per_vertex);
  int index = -1;

  checkMoabErr(iMOAB_DefineTagStorage(&app_id_, descriptor.name.c_str(), &type,
                                      &components, &index),
               "defining tag " + descriptor.name, comm_);

  tags_.emplace(descriptor.name, TagInfo{
                                     .moab_index = index,
                                     .descriptor = descriptor,
                                 });
}

void MoabAdapter::write_tag(const std::string& name, std::span<double> data) {
  int len = checked_size_to_int(data.size());

  TagMap::iterator it = tags_.find(name);
  if (it == tags_.end()) {
    throw std::runtime_error("Invalid tag name " + name);
  }
  MoabTagDescriptor tag = it->second.descriptor;
  std::size_t expected_size = tag.values_per_vertex * local_vertex_count_;

  if (expected_size != data.size()) {
    std::stringstream os;
    os << "Wrong size for tag " << name << ": Expected " << expected_size
       << " Got " << data.size();
    throw std::runtime_error(os.str());
  }

  int type = TAG_VERTEX;
  ErrCode err = iMOAB_SetDoubleTagStorage(&app_id_, name.c_str(), &len, &type,
                                          data.data());
  if (err) {
    fprintf(stderr, "Error: failed to create %s tag\n", name_.c_str());
    MPI_Abort(comm_, err);
  }
}

void MoabAdapter::read_tag(const std::string& name,
                           std::span<double> destination)  {

  const auto it = tags_.find(name);
  if (it == tags_.end()) {
    throw std::runtime_error("Unknown MOAB tag: " + name);
  }

  const std::size_t components = it->second.descriptor.values_per_vertex;

  if (components == 0 ||
      local_vertex_count_ >
          std::numeric_limits<std::size_t>::max() / components) {
    throw std::runtime_error("Invalid MOAB tag size: " + name);
  }

  const std::size_t expected = local_vertex_count_ * components;
  if (destination.size() != expected) {
    throw std::runtime_error("Destination size mismatch for MOAB tag: " + name);
  }

  int length = checked_size_to_int(destination.size());
  int entity_type = TAG_VERTEX;
  int my_id = app_id();

  ErrCode err = iMOAB_GetDoubleTagStorage(
      &my_id, name.c_str(), &length, &entity_type, destination.data());

  if (err != moab::MB_SUCCESS) {
    throw std::runtime_error("Failed to read MOAB tag: " + name);
  }
}

void MoabAdapter::create_vertices(const MoabPartition& partition) {

  std::size_t num_local_cols = partition.global_ids.size();
  auto lat = partition.latitudes;
  auto lon = partition.longitudes;

  std::vector<double> moab_vertex_coords(3 * num_local_cols);
  for (std::size_t i = 0; i < num_local_cols; ++i) {
    double lat_v = lat[i] * M_PI / 180.0;
    double lon_v = lon[i] * M_PI / 180.0;
    moab_vertex_coords[3 * i] = std::cos(lat_v) * std::cos(lon_v);
    moab_vertex_coords[3 * i + 1] = std::cos(lat_v) * std::sin(lon_v);
    moab_vertex_coords[3 * i + 2] = std::sin(lat_v);
  }

  int dimension = 3;
  int num_coords = checked_size_to_int(dimension * num_local_cols);

  ErrCode err = iMOAB_CreateVertices(&app_id_, &num_coords, &dimension,
                                     moab_vertex_coords.data());
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << name_
              << " model\n";
    MPI_Abort(comm_, err);
  }
}

void MoabAdapter::register_application(const MoabPartition& partition) {

  ErrCode err;
  err = iMOAB_RegisterApplication(name_.c_str(), &comm_, &component_id_,
                                  &app_id_);

  if (err) {
    MPI_Abort(comm(), err);
  }
  create_vertices(partition);
  checkMoabErr(iMOAB_UpdateMeshInfo(&app_id_), "updating mesh information",
               comm_);
  create_global_id_tag(partition);
}

void MoabAdapter::create_global_id_tag(const MoabPartition& part) {

  int type = 0; // dense, integer
  int stride = 1;
  int index;
  ErrCode err =
      iMOAB_DefineTagStorage(&app_id_, "GLOBAL_ID", &type, &stride, &index);
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << name_
              << " model\n";
    MPI_Abort(comm(), err);
  }

  auto gids = std::vector<int>{part.global_ids.begin(), part.global_ids.end()};
  int ent_type = 0; // entity type (vertex)
  int num_local_cols = checked_size_to_int(part.global_ids.size());

  err = iMOAB_SetIntTagStorage(&app_id_, "GLOBAL_ID", &num_local_cols,
                               &ent_type, gids.data());
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << name_
              << " model\n";
    MPI_Abort(comm(), err);
  }

  err = iMOAB_ResolveSharedEntities(&app_id_, &num_local_cols, gids.data());
  if (err) {
    std::cerr << "Error: failed to define GLOBAL_ID tag in " << name_
              << " model\n";
    MPI_Abort(comm(), err);
  }
}
} // namespace e3sm::coupler
