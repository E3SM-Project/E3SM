#ifndef E3SM_EMULATOR_MOAB_ADAPTOR_HPP
#define E3SM_EMULATOR_MOAB_ADAPTOR_HPP

#include "field_registry.hpp"
#include "mpi.h"
#include <map>
#include <moab/Types.hpp>
#include <moab/iMOAB.h>
#include <span>
#include <string>

namespace e3sm::coupler {
/*
  std::string name;
  FieldRole role;
  std::size_t values_per_entity;
  */
struct MoabTagDescriptor {
  std::string name;
  FieldRole role;
  std::size_t values_per_vertex;
};

/*
  MPI_Comm comm = MPI_COMM_NULL;
  std::span<int> global_ids;
  std::span<double> latitudes;
  std::span<double> longitudes;
  */
struct MoabPartition {
  MPI_Comm comm = MPI_COMM_NULL;
  std::span<int> global_ids;
  std::span<double> latitudes;
  std::span<double> longitudes;
};

struct NamedConstField {
  std::string_view name;
  std::span<const double> data;
};

struct NamedMutableField {
  std::string_view name;
  std::span<double> data;
};

class MoabAdapter {
public:
  MoabAdapter(std::string component_name, int component_id,
              const MoabPartition& partition);

  ~MoabAdapter();

  MoabAdapter(const MoabAdapter&) = delete;
  MoabAdapter& operator=(const MoabAdapter&) = delete;
  MoabAdapter(MoabAdapter&&) = delete;
  MoabAdapter& operator=(MoabAdapter&&) = delete;

  std::string_view name() const { return name_; }
  int app_id() const { return app_id_; }
  int component_id() const { return component_id_; }
  MPI_Comm comm() const { return comm_; }

  void write_exports(std::span<const NamedConstField> fields);
  void read_imports(std::span<const NamedMutableField> fields);
  void register_application(const MoabPartition& partition);
  void define_tag(const MoabTagDescriptor& desc);
  void write_tag(const std::string& name, std::span<double> data);
  void read_tag(const std::string& name, std::span<double> destination);

private:
  struct TagInfo {
    int moab_index;
    MoabTagDescriptor descriptor;
  };

  void create_vertices(const MoabPartition& partition);
  void assign_global_ids(std::span<const int> global_ids);
  void create_global_id_tag(const MoabPartition& part);

  std::string name_;
  int component_id_;
  int app_id_ = -1;
  MPI_Comm comm_ = MPI_COMM_NULL;

  std::size_t local_vertex_count_ = 0;
  using TagMap = std::map<std::string, TagInfo>;
  TagMap tags_;
};
} // namespace e3sm::coupler

#endif
