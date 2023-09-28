#include "share/grid/abstract_grid.hpp"

#include "share/field/field_utils.hpp"

#include <algorithm>
#include <cstring>
#include <string>

namespace scream
{
// Constructor(s) & Destructor
AbstractGrid::
AbstractGrid (const std::string& name,
              const GridType type,
              const int num_local_dofs,
              const int num_vertical_lev,
              const ekat::Comm& comm)
 : m_type (type)
 , m_name (name)
 , m_num_local_dofs (num_local_dofs)
 , m_num_vert_levs  (num_vertical_lev)
 , m_comm (comm)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (m_num_local_dofs>=0, "Error! Number of local dofs must be non-negative.\n");

  m_comm.all_reduce(&m_num_local_dofs,&m_num_global_dofs,1,MPI_SUM);

  // This grid name is also an alias
  m_aliases.push_back(m_name);
}

void AbstractGrid::add_alias (const std::string& alias)
{
  if (not ekat::contains(m_aliases,alias) and alias!=m_name) {
    m_aliases.push_back(alias);
  }
}

FieldLayout AbstractGrid::
get_vertical_layout (const bool midpoints) const
{
  using namespace ShortFieldTagsNames;
  return midpoints ? FieldLayout ({ LEV},{m_num_vert_levs})
                   : FieldLayout ({ILEV},{m_num_vert_levs+1});

}

bool AbstractGrid::is_unique () const {
  // Get a copy of gids on host. CAREFUL: do not use the stored dofs,
  // since we need to sort dofs in order to call unique, and we don't
  // want to alter the order of gids in this grid.
  auto dofs = m_dofs_gids.clone();
  auto dofs_h = dofs.get_view<gid_type*,Host>();

  std::sort(dofs_h.data(),dofs_h.data()+m_num_local_dofs);
  auto unique_end = std::unique(dofs_h.data(),dofs_h.data()+m_num_local_dofs);

  int locally_unique = unique_end==(dofs_h.data()+m_num_local_dofs);
  int unique;
  m_comm.all_reduce(&locally_unique,&unique,1,MPI_PROD);
  if (unique==0) {
    return false;
  }

  // Each rank has unique gids locally. Now it's time to verify if they are also globally unique.
  int max_dofs;
  m_comm.all_reduce(&m_num_local_dofs,&max_dofs,1,MPI_MAX);
  std::vector<gid_type> gids(max_dofs,-1);
  int unique_gids = 1;

  for (int pid=0; pid<m_comm.size(); ++pid) {
    // Rank pid broadcasts its gids, everyone else checks if there are duplicates
    if (pid==m_comm.rank()) {
      auto start = dofs_h.data();
      auto end   = start + m_num_local_dofs;
      std::copy(start,end,gids.data());
    }

    int ndofs = m_num_local_dofs;
    m_comm.broadcast(&ndofs,1,pid);
    m_comm.broadcast(gids.data(),ndofs,pid);

    int my_unique_gids = 1;
    if (pid!=m_comm.rank()) {
      // Checking two sorted arrays of length m and n for elements in common is O(m+n) ops.
      int i=0, j=0;
      while (i<m_num_local_dofs && j<ndofs && my_unique_gids==1) {
        if (dofs_h[i]<gids[j]) {
          ++i;
        } else if (dofs_h[i]>gids[j]) {
          ++j;
        } else {
          // Found a match. We can stop here
          my_unique_gids = 0;
          break;
        }
      }
    }
    m_comm.all_reduce(&my_unique_gids,&unique_gids,1,MPI_PROD);
    if (unique_gids==0) {
      break;
    }
  }

  return unique_gids;
}

auto AbstractGrid::
get_global_min_dof_gid () const ->gid_type
{
  // Lazy calculation
  if (m_global_min_dof_gid==std::numeric_limits<gid_type>::max()) {
    m_global_min_dof_gid = field_min<gid_type>(m_dofs_gids,&get_comm());
  }
  return m_global_min_dof_gid;
}

auto AbstractGrid::
get_global_max_dof_gid () const ->gid_type
{
  // Lazy calculation
  if (m_global_max_dof_gid==-std::numeric_limits<gid_type>::max()) {
    m_global_max_dof_gid = field_max<gid_type>(m_dofs_gids,&get_comm());
  }
  return m_global_max_dof_gid;
}

Field
AbstractGrid::get_dofs_gids () const {
  return m_dofs_gids.get_const();
}

Field
AbstractGrid::get_dofs_gids () {
  return m_dofs_gids;
}

Field
AbstractGrid::get_lid_to_idx_map () const {
  return m_lid_to_idx.get_const();
}

Field
AbstractGrid::get_lid_to_idx_map () {
  return m_lid_to_idx;
}

Field
AbstractGrid::get_geometry_data (const std::string& name) const {
  EKAT_REQUIRE_MSG (has_geometry_data(name),
      "Error! Geometry data '" + name + "' not found.\n");

  return m_geo_fields.at(name).get_const();
}

Field
AbstractGrid::create_geometry_data (const FieldIdentifier& fid)
{
  const auto& name = fid.name();

  EKAT_REQUIRE_MSG (not has_geometry_data(name),
      "Error! Cannot create geometry data, since it already exists.\n"
      "  - grid name: " + this->name() + "\n"
      "  - geo data name: " + name + "\n"
      "  - geo data layout: " + to_string(m_geo_fields.at(name).get_header().get_identifier().get_layout()) + "\n"
      "  - input layout: " + to_string(fid.get_layout()) + "\n");

  // Create field and the read only copy as well
  auto& f = m_geo_fields[name] = Field(fid);
  f.allocate_view();
  return f;
}

void
AbstractGrid::set_geometry_data (const Field& f)
{
  EKAT_REQUIRE_MSG (not has_geometry_data(f.name()),
      "Error! Cannot set geometry data, since it already exists.\n"
      "  - grid name: " + this->name() + "\n"
      "  - geo data name: " + f.name() + "\n");

  m_geo_fields[f.name()] = f;
}

std::list<std::string>
AbstractGrid::get_geometry_data_names () const
{
  std::list<std::string> names;
  for (const auto& it : m_geo_fields) {
    names.push_back(it.first);
  }
  return names;
}

void AbstractGrid::reset_num_vertical_lev (const int num_vertical_lev) {
  m_num_vert_levs = num_vertical_lev;

  // TODO: when the PR storing geo data as Field goes in, you should
  //       invalidate all geo data whose FieldLayout contains LEV/ILEV
}

std::vector<AbstractGrid::gid_type>
AbstractGrid::get_unique_gids () const
{
  // Gather local sizes across all ranks
  std::vector<int> ngids (m_comm.size());
  ngids[m_comm.rank()] = get_num_local_dofs();
  m_comm.all_gather(ngids.data(),1);

  std::vector<int> offsets (m_comm.size()+1,0);
  for (int pid=1; pid<=m_comm.size(); ++pid) {
    offsets[pid] = offsets[pid-1] + ngids[pid-1];
  }
  EKAT_REQUIRE_MSG (offsets[m_comm.size()]==m_num_global_dofs,
      "Error! Something went wrong while computing offsets in AbstractGrid::get_unique_grid.\n");

  // Gather all dofs
  const auto mpi_gid_t = ekat::get_mpi_type<gid_type>();
  std::vector<gid_type> all_gids (m_num_global_dofs);
  auto dofs_gids_h = m_dofs_gids.get_view<const gid_type*,Host>();
  MPI_Allgatherv (dofs_gids_h.data(),m_num_local_dofs,mpi_gid_t,
                  all_gids.data(),ngids.data(),offsets.data(),
                  mpi_gid_t,m_comm.mpi_comm());

  // Figure out unique dofs
  std::vector<gid_type> unique_dofs;
  const auto all_gids_beg = all_gids.begin();
  const auto all_gids_end = all_gids.begin() + offsets[m_comm.rank()];
  const auto my_gids_beg = dofs_gids_h.data();
  const auto my_gids_end = dofs_gids_h.data() + m_num_local_dofs;
  for (auto it=my_gids_beg; it!=my_gids_end; ++it) {
    if (std::find(all_gids_beg,all_gids_end,*it)==all_gids_end) {
      unique_dofs.push_back(*it);
    }
  }

  return unique_dofs;
}

std::vector<int> AbstractGrid::
get_owners (const gid_view_h& gids) const
{
  const auto& comm = get_comm();

  // Init owners to in
  std::map<gid_type,int> owners;
  int num_gids_in = gids.size();
  for (int i=0; i<num_gids_in; ++i) {
    owners[gids[i]] = -1;
  }
  int num_found = 0;

  // Let each rank bcast its owned gids, so that other procs can
  // check against their input list
  // Note: we can't use a view of const, since the view ptr needs to be passed
  //       to MPI bcast routines, which expect pointer to nonconst. It's an
  //       innocuous issue though, since only the send rank will use the ptr
  //       from the view, and it's not writing in it.
  auto my_gids_h = m_dofs_gids.get_view<gid_type*,Host>();
  gid_type* data;
  std::vector<gid_type> pid_gids;
  for (int pid=0; pid<comm.size(); ++pid) {
    // Bcast gids count for this pid
    int num_gids_pid = my_gids_h.size();
    comm.broadcast(&num_gids_pid,1,pid);

    // Bcast gids
    if (pid==comm.rank()) {
      data = my_gids_h.data();
    } else {
      pid_gids.resize(num_gids_pid);
      data = pid_gids.data();
    }
    comm.broadcast(data,num_gids_pid,pid);

    // Checks if any of the input pids is in this list
    for (int i=0; i<num_gids_pid && num_found<num_gids_in; ++i) {
      auto it = owners.find(data[i]);
      if (it!=owners.end()) {
        EKAT_REQUIRE_MSG (it->second==-1,
            "Error! Found a GID with multiple owners.\n"
            "  - owner 1: " + std::to_string(it->second) + "\n"
            "  - owner 2: " + std::to_string(pid) + "\n");
        it->second = pid;
        ++num_found;
      }
    }
  }
  EKAT_REQUIRE_MSG (num_found==static_cast<int>(owners.size()),
      "Error! Could not locate the owner of one of the input GIDs.\n"
      "  - rank: " + std::to_string(comm.rank()) + "\n"
      "  - num found: " + std::to_string(num_found) + "\n"
      "  - num gids in: " + std::to_string(num_gids_in) + "\n");


  // Now create and fill output view
  std::vector<int> result(num_gids_in);
  for (int i=0; i<num_gids_in; ++i) {
    result[i] = owners.at(gids[i]);
  }

  return result;
}

void AbstractGrid::create_dof_fields (const int scalar2d_layout_rank)
{
  using namespace ShortFieldTagsNames;
  const auto units = ekat::units::Units::nondimensional();

  // The dof gids field is a 1d field, while lid2idx has rank 2.
  // For both, the 1st dim is the num of local dofs. The 2nd dime of
  // lid2idx is the rank of a 2d scalar layout.
  FieldLayout dof_layout({COL},{get_num_local_dofs()});
  FieldLayout lid2idx_layout({COL,CMP},{get_num_local_dofs(),scalar2d_layout_rank});
  m_dofs_gids = Field(FieldIdentifier("gids",dof_layout,units,m_name,DataType::IntType));
  m_lid_to_idx = Field(FieldIdentifier("lid2idx",lid2idx_layout,units,m_name,DataType::IntType));

  m_dofs_gids.allocate_view();
  m_lid_to_idx.allocate_view();
}

void AbstractGrid::copy_data (const AbstractGrid& src, const bool shallow)
{
  if (shallow) {
    m_dofs_gids = src.m_dofs_gids;
  } else {
    m_dofs_gids = src.m_dofs_gids.clone();
  }

  if (shallow) {
    m_lid_to_idx = src.m_lid_to_idx;
  } else {
    m_lid_to_idx = src.m_lid_to_idx.clone();
  }

  for (const auto& name : src.get_geometry_data_names()) {
    if (shallow) {
      m_geo_fields[name] = src.m_geo_fields.at(name);
    } else {
      m_geo_fields[name] = src.m_geo_fields.at(name).clone();
    }
  }
}

} // namespace scream
