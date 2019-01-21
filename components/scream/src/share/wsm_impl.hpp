#ifndef WSM_IMPL_HPP
#define WSM_IMPL_HPP

#include "util.hpp"

namespace util {

/*
 * An implementation header for wsm.hpp, this helps keep wsm.hpp clean by avoiding
 * mixing interface and implementation. Clients should NOT even include this file.
 */

template <typename T, typename D>
WorkspaceManager<T, D>::WorkspaceManager(int size, int max_used, TeamPolicy policy) :
  m_tu(policy),
  m_concurrent_teams(m_tu.get_num_concurrent_teams()),
  m_reserve( (sizeof(T) > 2*sizeof(int)) ? 1 :
             (2*sizeof(int) + sizeof(T) - 1)/sizeof(T) ),
  m_size(size),
  m_total(m_size + m_reserve),
  m_max_used(max_used),
#ifndef NDEBUG
  m_num_used("Workspace.m_num_used", m_concurrent_teams),
  m_high_water("Workspace.m_high_water", m_concurrent_teams),
  m_active("Workspace.m_active", m_concurrent_teams, m_max_used),
  m_curr_names("Workspace.m_curr_names", m_concurrent_teams, m_max_used, m_max_name_len),
  m_all_names("Workspace.m_all_names", m_concurrent_teams, m_max_names, m_max_name_len),
  // A name's index in m_all_names is used to index into m_counts
  m_counts("Workspace.m_counts", m_concurrent_teams, m_max_names, 2),
#endif
  m_next_slot("Workspace.m_next_slot", m_pad_factor*m_concurrent_teams),
  m_data(Kokkos::ViewAllocateWithoutInitializing("Workspace.m_data"),
         m_concurrent_teams, m_total * m_max_used)
{
  init(*this, m_data, m_concurrent_teams, m_max_used, m_total);
}

template <typename T, typename D>
void WorkspaceManager<T, D>::report() const
{
#ifndef NDEBUG
  auto host_num_used   = Kokkos::create_mirror_view(m_num_used);
  auto host_high_water = Kokkos::create_mirror_view(m_high_water);
  auto host_all_names  = Kokkos::create_mirror_view(m_all_names);
  auto host_counts     = Kokkos::create_mirror_view(m_counts);

  std::cout << "\nWS usage (capped at " << m_max_used << "): " << std::endl;
  for (int t = 0; t < m_concurrent_teams; ++t) {
    std::cout << "WS " << t << " currently using " << host_num_used(t) << std::endl;
    std::cout << "WS " << t << " high-water " << host_high_water(t) << std::endl;
  }

  std::cout << "\nWS deep analysis" << std::endl;
  std::map<std::string, std::tuple<int, int, int> > ws_usage_map;
  for (int t = 0; t < m_concurrent_teams; ++t) {
    std::cout << "  For wsidx " << t << std::endl;
    for (int n = 0; n < m_max_names; ++n) {
      const char* name = &(host_all_names(t, n, 0));
      if (util::strcmp(name, "") == 0) {
        break;
      }
      else {
        const int takes    = host_counts(t, n, 0);
        const int releases = host_counts(t, n, 1);
        std::cout << "    workspace '" << name << "' was taken " << takes
                  << " times and released " << releases << " times" << std::endl;
        if (takes != releases) {
          std::cout << "      POSSIBLE LEAK" << std::endl;
        }
        std::string sname(name);
        if (ws_usage_map.find(sname) == ws_usage_map.end()) {
          ws_usage_map[sname] = std::make_tuple(1, takes, releases);
        }
        else {
          std::get<0>(ws_usage_map[sname]) += 1;
          std::get<1>(ws_usage_map[sname]) += takes;
          std::get<2>(ws_usage_map[sname]) += releases;
        }
      }
    }
  }

  std::cout << "\nWS workspace summary" << std::endl;
  for (auto& kv : ws_usage_map) {
    auto data = kv.second;
    std::cout << "Workspace '" << kv.first << "' was used by " << std::get<0>(data) << " wsindices with "
              << std::get<1>(data) << " takes and " << std::get<2>(data) << " releases." << std::endl;
  }
#endif
}

template <typename T, typename D>
KOKKOS_INLINE_FUNCTION
typename WorkspaceManager<T, D>::Workspace
WorkspaceManager<T, D>::get_workspace(const MemberType& team) const
{ return Workspace(*this, m_tu.get_workspace_idx(team), team); }

template <typename T, typename D>
void WorkspaceManager<T, D>::init(const WorkspaceManager<T, D>& wm, const view_2d<T>& data,
                                  const int concurrent_teams, const int max_used, const int total)
{
  Kokkos::parallel_for(
    "WorkspaceManager ctor",
    util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(concurrent_teams, max_used),
    KOKKOS_LAMBDA(const MemberType& team) {
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, max_used), [&] (int i) {
          wm.init_metadata(team.league_rank(), i);
        });
    });
}

template <typename T, typename D>
template <typename S>
KOKKOS_FORCEINLINE_FUNCTION
int WorkspaceManager<T, D>::set_next_and_get_index(const Unmanaged<view_1d<S> >& space, int next) const
{
  const auto metadata = reinterpret_cast<int*>(reinterpret_cast<T*>(space.data()) - m_reserve);
  metadata[1] = next;
  return metadata[0];
}

template <typename T, typename D>
template <typename S>
KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<typename WorkspaceManager<T, D>::template view_1d<S> >
WorkspaceManager<T, D>::get_space_in_slot(const int team_idx, const int slot) const
{
  return Unmanaged<view_1d<S> >(
    reinterpret_cast<S*>(&m_data(team_idx, slot*m_total) + m_reserve),
    sizeof(T) == sizeof(S) ?
    m_size :
    (m_size*sizeof(T))/sizeof(S));
}

template <typename T, typename D>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::init_metadata(const int ws_idx, const int slot) const
{
  int* const metadata = reinterpret_cast<int*>(&m_data(ws_idx, slot*m_total));
  metadata[0] = slot;     // idx
  metadata[1] = slot + 1; // next
}

template <typename T, typename D>
KOKKOS_INLINE_FUNCTION
WorkspaceManager<T, D>::Workspace::Workspace(
  const WorkspaceManager& parent, int ws_idx, const MemberType& team) :
  m_parent(parent), m_team(team), m_ws_idx(ws_idx),
  m_next_slot(parent.m_next_slot(m_pad_factor*ws_idx))
{}

template <typename T, typename D>
template <typename S>
KOKKOS_INLINE_FUNCTION
Unmanaged<typename WorkspaceManager<T, D>::template view_1d<S> > WorkspaceManager<T, D>::Workspace::take(
  const char* name) const
{
#ifndef NDEBUG
  change_num_used(1);
#endif

  const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, m_next_slot);

  // We need a barrier here so get_space_in_slot returns consistent results
  // w/in the team.
  m_team.team_barrier();
  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
    m_next_slot = m_parent.get_next<S>(space);
#ifndef NDEBUG
    change_indv_meta<S>(space, name);
#endif
  });
  // We need a barrier here so that a subsequent call to take or release
  // starts with the metadata in the correct state.
  m_team.team_barrier();

  return space;
}

template <typename T, typename D>
template <size_t N, typename S>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::Workspace::take_many_contiguous_unsafe(
  const Kokkos::Array<const char*, N>& names,
  const view_1d_ptr_array<S, N>& ptrs) const
{
#ifndef NDEBUG
  change_num_used(N);
  // Verify contiguous
  for (int n = 0; n < static_cast<int>(N); ++n) {
    const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, m_next_slot + n);
    micro_kassert(m_parent.get_next<S>(space) == m_next_slot + n + 1);
  }
#endif

  for (int n = 0; n < static_cast<int>(N); ++n) {
    const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, m_next_slot+n);
    *ptrs[n] = space;
  }

  // We need a barrier here so get_space_in_slot above returns consistent results
  // w/in the team.
  m_team.team_barrier();
  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
    m_next_slot += N;
#ifndef NDEBUG
    for (int n = 0; n < static_cast<int>(N); ++n) {
      change_indv_meta<S>(*ptrs[n], names[n]);
    }
#endif
  });
  // We need a barrier here so that a subsequent call to take or release
  // starts with the metadata in the correct state.
  m_team.team_barrier();
}

template <typename T, typename D>
template <size_t N, typename S>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::Workspace::take_many(
  const Kokkos::Array<const char*, N>& names,
  const view_1d_ptr_array<S, N>& ptrs) const
{
#ifndef NDEBUG
  change_num_used(N);
#endif

  int next_slot = m_next_slot;
  for (int n = 0; n < static_cast<int>(N); ++n) {
    auto& space = *ptrs[n];
    space = m_parent.get_space_in_slot<S>(m_ws_idx, next_slot);
    next_slot = m_parent.get_next<S>(space);
  }

  // We need a barrier here so get_space_in_slot above returns consistent results
  // w/in the team.
  m_team.team_barrier();
  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
    m_next_slot = next_slot;
#ifndef NDEBUG
    for (int n = 0; n < static_cast<int>(N); ++n) {
      change_indv_meta<S>(*ptrs[n], names[n]);
    }
#endif
  });
  // We need a barrier here so that a subsequent call to take or release
  // starts with the metadata in the correct state.
  m_team.team_barrier();
}

template <typename T, typename D>
template <size_t N, typename S>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::Workspace::take_many_and_reset(
  const Kokkos::Array<const char*, N>& names,
  const view_1d_ptr_array<S, N>& ptrs) const
{
#ifndef NDEBUG
  change_num_used(N - m_parent.m_num_used(m_ws_idx));
#endif

  for (int n = 0; n < static_cast<int>(N); ++n) {
    const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, n);
    *ptrs[n] = space;
  }

  // We only need to reset the metadata for spaces that are being left free
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(m_team, m_parent.m_max_used - N), [&] (int i) {
      m_parent.init_metadata(m_ws_idx, i+N);
    });

  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
    m_next_slot = N;
#ifndef NDEBUG
    // Mark all old spaces as released
    for (int a = 0; a < m_parent.m_max_used; ++a) {
      if (m_parent.m_active(m_ws_idx, a)) {
        change_indv_meta<S>(m_parent.get_space_in_slot<S>(m_ws_idx, a), "", true);
      }
    }

    // Mark all new spaces as taken
    for (int n = 0; n < static_cast<int>(N); ++n) {
      change_indv_meta<S>(*ptrs[n], names[n]);
    }
#endif
    });
  // We need a barrier here so that a subsequent call to take or release
  // starts with the metadata in the correct state.
  m_team.team_barrier();
}

template <typename T, typename D>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::Workspace::reset() const
{
  m_team.team_barrier();
#ifndef NDEBUG
  change_num_used(-m_parent.m_num_used(m_ws_idx));
#endif
  m_next_slot = 0;
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(m_team, m_parent.m_max_used), [&] (int i) {
      m_parent.init_metadata(m_ws_idx, i);
    });

#ifndef NDEBUG
  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
    // Mark all old spaces as released
    for (int a = 0; a < m_parent.m_max_used; ++a) {
      if (m_parent.m_active(m_ws_idx, a)) {
        change_indv_meta<T>(m_parent.get_space_in_slot<T>(m_ws_idx, a), "", true);
      }
    }
  });
#endif

  m_team.team_barrier();
}

// Print the linked list. Obviously not a device function.
template <typename T, typename D>
void WorkspaceManager<T, D>::Workspace::print() const
{
  m_team.team_barrier();
  Kokkos::single(
    Kokkos::PerTeam(m_team), [&] () {
      std::stringstream ss;
      ss << m_ws_idx << ":";
      auto space = m_parent.get_space_in_slot<T>(m_ws_idx, m_next_slot);
      for (int cnt = 0, nmax = m_parent.m_max_used;
           cnt < nmax;
           ++cnt) {
        ss << " (" << m_parent.get_index<T>(space) << ", "
           << m_parent.get_next<T>(space) << ")";
        space = m_parent.get_space_in_slot<T>(m_ws_idx, m_parent.get_next<T>(space));
      }
      ss << "\n";
      std::cout << ss.str();
    });
}

#ifndef NDEBUG
template <typename T, typename D>
template <typename S>
KOKKOS_INLINE_FUNCTION
const char* WorkspaceManager<T, D>::Workspace::get_name_impl(const Unmanaged<view_1d<S> >& space) const
{
  const int slot = m_parent.get_index<S>(space);
  return &(m_parent.m_curr_names(m_ws_idx, slot, 0));
}

template <typename T, typename D>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::Workspace::change_num_used(int change_by) const
{
  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
    int curr_used = m_parent.m_num_used(m_ws_idx) += change_by;
    micro_kassert(curr_used <= m_parent.m_max_used);
    micro_kassert(curr_used >= 0);
    if (curr_used > m_parent.m_high_water(m_ws_idx)) {
      m_parent.m_high_water(m_ws_idx) = curr_used;
    }
  });
}

template <typename T, typename D>
template <typename S>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::Workspace::change_indv_meta(
  const Unmanaged<view_1d<S> >& space, const char* name, bool release) const
{
  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
    const int slot = m_parent.get_index<S>(space);
    if (!release) {
      micro_kassert(util::strlen(name) < m_max_name_len); // leave one char for null terminator
      micro_kassert(util::strlen(name) > 0);
      micro_kassert(!m_parent.m_active(m_ws_idx, slot));
      char* val = &(m_parent.m_curr_names(m_ws_idx, slot, 0));
      util::strcpy(val, name);
    }
    else {
      micro_kassert(m_parent.m_active(m_ws_idx, slot));
      name = get_name(space);
    }
    const int name_idx = get_name_idx(name, !release);
    const int count_idx = release ? 1 : 0;
    m_parent.m_counts(m_ws_idx, name_idx, count_idx) += 1;
    m_parent.m_active(m_ws_idx, slot) = !release;
  });
}

template <typename T, typename D>
KOKKOS_INLINE_FUNCTION
int WorkspaceManager<T, D>::Workspace::get_name_idx(const char* name, bool add) const
{
  int name_idx = -1;
  for (int n = 0; n < m_max_names; ++n) {
    char* old_name = &(m_parent.m_all_names(m_ws_idx, n, 0));
    if (util::strcmp(old_name, name) == 0) {
      name_idx = n;
      break;
    }
    else if (add && util::strcmp(old_name, "") == 0) {
      util::strcpy(old_name, name);
      name_idx = n;
      break;
    }
  }
  micro_kassert(name_idx != -1);
  return name_idx;
}
#endif

template <typename T, typename D>
template <typename S>
KOKKOS_INLINE_FUNCTION
void WorkspaceManager<T, D>::Workspace::release_impl(const Unmanaged<view_1d<S> >& space) const
{
#ifndef NDEBUG
  change_num_used(-1);
  change_indv_meta<S>(space, "", true);
#endif

  // We don't need a barrier before this block b/c it's OK for metadata to
  // change while some threads in the team are still using the bulk data.
  Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
      m_next_slot = m_parent.set_next_and_get_index<S>(space, m_next_slot);
  });
  m_team.team_barrier();
}

} // namespace util

#endif
