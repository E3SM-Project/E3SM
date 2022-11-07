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

AbstractGrid::gid_type
AbstractGrid::get_global_min_dof_gid () const
{
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! You need to set dofs gids before you can compute the global min dof.\n");
  // TODO: we could probably cache these into mutable variables.
  //       But unless we call this method *many* times, it won't matter
  gid_type local_min, global_min;
  auto dofs = get_dofs_gids();
  Kokkos::parallel_reduce(Kokkos::RangePolicy<>(0,get_num_local_dofs()),
      KOKKOS_LAMBDA (const int& i, gid_type& lmin) {
        if (dofs(i) < lmin) {
          lmin = dofs(i);
        }
      },Kokkos::Min<gid_type>(local_min));
  Kokkos::fence();

  m_comm.all_reduce(&local_min,&global_min,1,MPI_MIN);

  return global_min;
}

AbstractGrid::gid_type
AbstractGrid::get_global_max_dof_gid () const
{
  EKAT_REQUIRE_MSG (m_dofs_set,
      "Error! You need to set dofs gids before you can compute the global max dof.\n");
  // TODO: we could probably cache these into mutable variables.
  //       But unless we call this method *many* times, it won't matter
  gid_type local_max, global_max;
  auto dofs = get_dofs_gids();
  Kokkos::parallel_reduce(Kokkos::RangePolicy<>(0,get_num_local_dofs()),
      KOKKOS_LAMBDA (const int& i, gid_type& lmax) {
        if (dofs(i) > lmax) {
          lmax = dofs(i);
        }
      },Kokkos::Max<gid_type>(local_max));
  Kokkos::fence();

  m_comm.all_reduce(&local_max,&global_max,1,MPI_MAX);

  return global_max;
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
  // In order to ship information around across ranks, it is easier to use
  // an auxiliary grid, where dofs are partitioned across ranks linearly
  // NOTE: we actually don't need the grid itself. We only need to know
  //       what the local number of dofs would be on this rank.
  const int ngdofs = get_num_global_dofs();
  const auto& comm = get_comm();
  int nldofs_linear = ngdofs / comm.size();
  if (comm.rank()<(ngdofs % comm.size())) {
    ++ nldofs_linear;
  }

  // For each pid, compute offsets in the global gids array.
  std::vector<int> offsets(comm.size()+1,0);
  const int ndofs_per_rank = ngdofs / comm.size();
  const int remainder = ngdofs % comm.size();
  for (int pid=1; pid<=comm.size(); ++pid) {
    offsets[pid] = offsets[pid-1] + ndofs_per_rank;
    if ( (pid-1)<remainder ){
      ++offsets[pid];
    }
  }
  EKAT_REQUIRE_MSG (offsets.back()==ngdofs,
      "Error! Something went wrong while calling get_gids_owners.\n"
      "  - grid name: " + this->name() + "\n");

  // Utility lambda: given a GID, retrieve the PID that would own it in a
  // linearly distributed grid, as well as the corresponding LID it would
  // have in that grid on that rank. This is doable without any communication
  // since the gids are partitioned linearly
  auto pid_and_lid = [&] (const gid_type gid) -> std::pair<int,int> {
    auto it = std::upper_bound (offsets.begin(),offsets.end(),gid);
    int pid = std::distance(offsets.begin(),it) - 1;
    int lid = gid - offsets[pid];
    EKAT_REQUIRE_MSG (pid>=0 && pid<comm.size(),
        "Error! Failed to retrieve owner of GID in the linear grid.\n");
    return std::make_pair(pid,lid);
  };

  // The idea is to create a "global" array (partitioned across ranks) of the
  // gids in the linear map, use it to store the owner of each dofs (in the original grid),
  // and finally read that global array for all the input gids.
  struct PidLid {
    int pid;
    int lid;
  };
  std::map<int,std::vector<PidLid>> rma_data;
  std::map<int,std::vector<int>> rma_offsets;
  std::map<int,MPI_Datatype> rma_dtypes;

  auto clear_dtypes = [&] () {
    for (auto& it : rma_dtypes) {
      MPI_Type_free(&it.second);
    }
    rma_dtypes.clear();
  };

  MPI_Win win;
  PidLid* pids_lids_linear;
  MPI_Win_allocate (nldofs_linear*sizeof(MPI_2INT),sizeof(MPI_2INT),
                    MPI_INFO_NULL,comm.mpi_comm(),&pids_lids_linear,&win);
  MPI_Win_fence(0,win);

  // Step 1: each rank loops through its grid gids, and sets owners_linear=rank
  // for all its gids in the grid.
  MPI_Win_fence(0,win);

  //  - 1.a Figure where each local dof specs will be written in the linearly
  //        distributed global array
  for (int i=0; i<get_num_local_dofs(); ++i) {
    const auto gid = m_dofs_gids_host[i];
    const auto pidlid = pid_and_lid (gid);
    const auto pid = pidlid.first;
    const auto lid = pidlid.second;
    rma_data[pid].push_back({comm.rank(),i});
    rma_offsets[pid].push_back(lid);
  }

  //  - 1.b: build data types for writing all the data on each tgt rank at once
  for (const auto& it : rma_offsets) {
    auto& dtype = rma_dtypes[it.first];
    std::vector<int> ones(it.second.size(),1);
    MPI_Type_indexed (it.second.size(),ones.data(),it.second.data(),MPI_2INT,&dtype);
    MPI_Type_commit (&dtype);
  }

  //  - 1.c: write on the window
  for (const auto& it : rma_data) {
    const auto pid = it.first;
    const auto dtype = rma_dtypes.at(pid);
    // Note: the dtype already encodes offsets in the remote window, so tgt_disp=0
    MPI_Put (it.second.data(),it.second.size(),MPI_2INT,pid,
             0,1,dtype,win);
  }
  clear_dtypes();
  rma_data.clear();
  rma_offsets.clear();
  MPI_Win_fence(0,win);

  // Step 2: each rank loops over its input gids, and retrieves the owner from the window

  //  - 1.a: Figure out what needs to be read from each rank
  for (int i=0; i<gids.extent_int(0); ++i) {
    const auto gid = gids[i];
    const auto pidlid = pid_and_lid (gid);
    const auto pid = pidlid.first;
    const auto lid = pidlid.second;
    rma_offsets[pid].push_back(lid);
    rma_data[pid].push_back({-1,-1});
  }

  //  - 1.b: build data types for reading all the data from each tgt rank at once
  for (const auto& it : rma_offsets) {
    auto& dtype = rma_dtypes[it.first];
    std::vector<int> ones(it.second.size(),1);
    MPI_Type_indexed (it.second.size(),ones.data(),it.second.data(),MPI_2INT,&dtype);
    MPI_Type_commit (&dtype);
  }

  //  - 1.c: read from the window
  for (auto& it : rma_data) {
    const auto pid = it.first;
    const auto dtype = rma_dtypes.at(pid);
    // Note: the dtype already encodes offsets in the remote window, so tgt_disp=0
    MPI_Get (it.second.data(),it.second.size(),MPI_2INT,pid,
             0,1,dtype,win);
  }
  clear_dtypes();
  rma_offsets.clear();
  MPI_Win_fence(0,win);

  // Step 3: copy data in rma types into output view, making sure we keep correct order
  hview_1d<int> owners("",gids.size());
  std::map<int,int> curr_data_index;
  for (int i=0; i<gids.extent_int(0); ++i) {
    const auto gid = gids[i];
    const auto pidlid = pid_and_lid (gid);
    const auto pid = pidlid.first;
    auto it_bool = curr_data_index.emplace(pid,0);
    auto& idx = it_bool.first->second;
    owners(i) = rma_data[pid][idx].pid;
    ++idx;
  }
  rma_data.clear();

  // Clean up
  MPI_Win_free(&win);

  return owners;
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
