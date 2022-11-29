#include "share/grid/abstract_grid.hpp"

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

bool AbstractGrid::is_unique () const {
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! We cannot establish if this grid is unique before the dofs GIDs are set.\n");

  // Get a copy of gids on host. CAREFUL: do not use create_mirror_view,
  // since it would create a shallow copy on CPU devices, but we need a
  // deep copy, to prevent altering order of gids.
  decltype(m_dofs_gids)::HostMirror dofs_h("",m_dofs_gids.size());
  Kokkos::deep_copy(dofs_h,m_dofs_gids);

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
  std::vector<int> gids(max_dofs);
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

void AbstractGrid::
set_dofs (const dofs_list_type& dofs)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_dofs_set, "Error! Dofs cannot be re-set, once set.\n");

  EKAT_REQUIRE_MSG (dofs.extent_int(0)==m_num_local_dofs,
      "Error! Wrong size for the input dofs list. It should match the number of local dofs stored in the mesh.\n"
      "       Expected gids list size: " + std::to_string(m_num_local_dofs) + "\n"
      "       Input gids list size   : " + std::to_string(dofs.size()) + "\n");

  m_dofs_gids  = dofs;
  m_dofs_set = true;

#ifndef NDEBUG
  EKAT_REQUIRE_MSG(this->valid_dofs_list(dofs), "Error! Invalid list of dofs gids.\n");
#endif

  m_dofs_gids_host = Kokkos::create_mirror_view(m_dofs_gids);
  Kokkos::deep_copy(m_dofs_gids_host,m_dofs_gids);

  set_global_min_dof_gid();
  set_global_max_dof_gid();
}

const AbstractGrid::dofs_list_type&
AbstractGrid::get_dofs_gids () const {
  // Sanity check
  EKAT_REQUIRE_MSG (m_dofs_set, "Error! You must call 'set_dofs' first.\n");

  return m_dofs_gids;
}

const AbstractGrid::dofs_list_h_type&
AbstractGrid::get_dofs_gids_host () const {
  // Sanity check
  EKAT_REQUIRE_MSG (m_dofs_set, "Error! You must call 'set_dofs' first.\n");

  return m_dofs_gids_host;
}


void AbstractGrid::
set_lid_to_idx_map (const lid_to_idx_map_type& lid_to_idx)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_lid_to_idx_set, "Error! The lid->idx map cannot be re-set, once set.\n");

  EKAT_REQUIRE_MSG (lid_to_idx.extent_int(0)==m_num_local_dofs &&
                    lid_to_idx.extent_int(1)==get_2d_scalar_layout().rank(),
      "Error! Wrong size(s) for the input lid_to_idx map. They should match (num_local_dofs, 2d layout rank).\n"
      "       Expected sizes: (" + std::to_string(m_num_local_dofs) + "," + std::to_string(get_2d_scalar_layout().rank()) + ")\n"
      "       Input map sizes: (" + std::to_string(lid_to_idx.extent(0)) + "," + std::to_string(lid_to_idx.extent(1)) + ")\n");

#ifndef NDEBUG
  EKAT_REQUIRE_MSG(this->valid_lid_to_idx_map(lid_to_idx), "Error! Invalid lid->idx map.\n");
#endif

  m_lid_to_idx = lid_to_idx;
  m_lid_to_idx_set = true;
}

const AbstractGrid::lid_to_idx_map_type&
AbstractGrid::get_lid_to_idx_map () const {
  // Sanity check
  EKAT_REQUIRE_MSG (m_dofs_gids.size()>0, "Error! You must call 'set_dofs' first.\n");

  return m_lid_to_idx;
}

void AbstractGrid::
set_geometry_data (const std::string& name, const geo_view_type& data) {
  m_geo_views[name] = data;
  m_geo_views_host[name] = Kokkos::create_mirror_view(data);
  Kokkos::deep_copy(m_geo_views_host[name],data);
}

const AbstractGrid::geo_view_type&
AbstractGrid::get_geometry_data (const std::string& name) const {
  EKAT_REQUIRE_MSG (m_geo_views.find(name)!=m_geo_views.end(),
                    "Error! Grid '" + m_name + "' does not store geometric data '" + name + "'.\n");
  return m_geo_views.at(name);
}

const AbstractGrid::geo_view_h_type&
AbstractGrid::get_geometry_data_host (const std::string& name) const {
  EKAT_REQUIRE_MSG (m_geo_views_host.find(name)!=m_geo_views_host.end(),
                    "Error! Grid '" + m_name + "' does not store geometric data '" + name + "'.\n");
  return m_geo_views_host.at(name);
}

void
AbstractGrid::set_global_min_dof_gid ()
{
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! You need to set dofs gids before you can compute the global min dof.\n");
  gid_type local_min = get_num_global_dofs(); // the local min should be <= than the size of the grid.
  gid_type global_min;
  if (get_num_local_dofs()>0) {
    auto dofs = get_dofs_gids();
    Kokkos::parallel_reduce(Kokkos::RangePolicy<>(0,get_num_local_dofs()),
        KOKKOS_LAMBDA (const int& i, gid_type& lmin) {
          if (dofs(i) < lmin) {
            lmin = dofs(i);
          }
        },Kokkos::Min<gid_type>(local_min));
    Kokkos::fence();
  }

  m_comm.all_reduce(&local_min,&global_min,1,MPI_MIN);

  m_global_min_dof_gid = global_min;
}

void
AbstractGrid::set_global_max_dof_gid ()
{
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! You need to set dofs gids before you can compute the global max dof.\n");
  gid_type local_max = -1;
  gid_type global_max;
  if (get_num_local_dofs()>0) {
    auto dofs = get_dofs_gids();
    Kokkos::parallel_reduce(Kokkos::RangePolicy<>(0,get_num_local_dofs()),
        KOKKOS_LAMBDA (const int& i, gid_type& lmax) {
          if (dofs(i) > lmax) {
            lmax = dofs(i);
          }
        },Kokkos::Max<gid_type>(local_max));
    Kokkos::fence();
  }

  m_comm.all_reduce(&local_max,&global_max,1,MPI_MAX);

  m_global_max_dof_gid = global_max;
}

std::list<std::string>
AbstractGrid::get_geometry_data_names () const
{
  std::list<std::string> names;
  for (const auto& it : m_geo_views) {
    names.push_back(it.first);
  }
  return names;
}

void AbstractGrid::reset_num_vertical_lev (const int num_vertical_lev) {
  m_num_vert_levs = num_vertical_lev;

  // TODO: when the PR storing geo data as Field goes in, you should
  //       invalidate all geo data whose FieldLayout contains LEV/ILEV
}

auto
AbstractGrid::get_unique_gids () const
 -> dofs_list_type
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
  MPI_Allgatherv (m_dofs_gids_host.data(),m_num_local_dofs,mpi_gid_t,
                  all_gids.data(),ngids.data(),offsets.data(),
                  mpi_gid_t,m_comm.mpi_comm());

  // Figure out unique dofs
  std::vector<gid_type> unique_dofs;
  const auto all_gids_beg = all_gids.begin();
  const auto all_gids_end = all_gids.begin() + offsets[m_comm.rank()];
  const auto my_gids_beg = m_dofs_gids_host.data();
  const auto my_gids_end = m_dofs_gids_host.data() + m_num_local_dofs;
  for (auto it=my_gids_beg; it!=my_gids_end; ++it) {
    if (std::find(all_gids_beg,all_gids_end,*it)==all_gids_end) {
      unique_dofs.push_back(*it);
    }
  }

  dofs_list_type unique_gids_d("",unique_dofs.size());
  auto unique_gids_h = Kokkos::create_mirror_view(unique_gids_d);
  std::memcpy(unique_gids_h.data(),unique_dofs.data(),sizeof(gid_type)*unique_dofs.size());
  Kokkos::deep_copy(unique_gids_d,unique_gids_h);
  return unique_gids_d;
}

auto AbstractGrid::
get_owners (const hview_1d<const gid_type>& gids) const
 -> hview_1d<int>
{
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! Cannot retrieve gids owners until dofs gids have been set.\n");

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
  auto my_gids_h = m_dofs_gids_host;
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
  EKAT_REQUIRE_MSG (num_found==owners.size(),
      "Error! Could not locate the owner of one of the input GIDs.\n"
      "  - rank: " + std::to_string(comm.rank()) + "\n"
      "  - num found: " + std::to_string(num_found) + "\n"
      "  - num gids in: " + std::to_string(num_gids_in) + "\n");


  // Now create and fill output view
  hview_1d<int> result("",num_gids_in);
  for (int i=0; i<num_gids_in; ++i) {
    result[i] = owners.at(gids[i]);
  }

  return result;
}

void AbstractGrid::copy_views (const AbstractGrid& src, const bool shallow)
{
  if (src.m_dofs_set) {
    m_dofs_set = false;
    if (shallow) {
      set_dofs (src.m_dofs_gids);
      // set_dof created a new host mirror, so assing that too
      m_dofs_gids_host = src.m_dofs_gids_host;
    } else {
      decltype (src.m_dofs_gids) dofs ("",m_num_local_dofs);
      Kokkos::deep_copy (dofs,src.get_dofs_gids());
      set_dofs (dofs);
    }
  }

  if (src.m_lid_to_idx_set) {
    m_lid_to_idx_set = false;
    if (shallow) {
      set_lid_to_idx_map (src.get_lid_to_idx_map());
    } else {
      decltype (src.m_lid_to_idx) lid2idx ("",m_num_local_dofs,get_2d_scalar_layout().rank());
      Kokkos::deep_copy (lid2idx,src.m_lid_to_idx);
      set_lid_to_idx_map (lid2idx);
    }
  }

  for (const auto& it : src.m_geo_views) {
    const auto& name = it.first;
    if (shallow) {
      m_geo_views[name] = it.second;
    } else {
      m_geo_views[name] = geo_view_type("",it.second.size());
      Kokkos::deep_copy (m_geo_views[name],it.second);
    }
  }

  for (const auto& it : src.m_geo_views_host) {
    const auto& name = it.first;
    if (shallow) {
      m_geo_views_host[name] = it.second;
    } else {
      m_geo_views_host[name] = geo_view_h_type("",it.second.size());
      Kokkos::deep_copy (m_geo_views_host[name],it.second);
    }
  }
}

} // namespace scream
