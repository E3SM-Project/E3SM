#ifndef WSM_HPP
#define WSM_HPP

#include "types.hpp"
#include "scream_arch.hpp"
#include "kokkos_util.hpp"

namespace unit_test {
struct UnitWrap;
}

namespace util {

/*
 * WorkspaceManager is a utility for requesting workspaces
 * (temporary memory blocks) from within a Kokkos kernel. Workspaces
 * are shared within thread teams.
 *
 * The WorkspaceManager should be initialized before Kokkos kernels
 * are run.  Users need to specify the size of each individual
 * sub-block (in terms of number of T's), the maximum number of
 * sub-blocks (AKA "spaces") they intend to use, and the team policy
 * that they intend to use for their Kokkos kernels.
 *
 * Once inside a Kokkos kernel, the user will call get_workspace on
 * their WorkspaceManager. This will return a Workspace object which
 * is the user's iterface to the entire memory block, which mostly
 * involves taking/releasing memory sub-blocks (represented by kokkos
 * views). The API for taking/releasing sub-blocks has a number of
 * variants in order to offer both simplicity/clarity and high
 * performance. The simplest way to use your Workspace is to manage
 * your sub-blocks individually with the "take" and "release"
 * methods. You'll get better performance with the methods that deal
 * with the sub-blocks in bulk.  We've seen the best performance by
 * making a single take_many_and_reset call at the beginning of the
 * Kokkos kernel, though this may not be practical for large
 * kernels. High-granularity take/release calls can allow you to use a
 * smaller Workspace (by allowing you to use a smaller max_used).
 *
 * You can fine-tune the max_used by making a large overestimate for
 * this value, running your kernel, and then calling report, which
 * will tell you the actual maximum number of sub-blocks that you
 * used. Note that all sub-blocks have a name.
 */

template <typename T, typename DeviceT=DefaultDevice>
class WorkspaceManager
{
 public:

  //
  // ---------- Types ---------
  //

  using Device = DeviceT;

  using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;
  using MemberType = typename KokkosTypes<Device>::MemberType;
  using ExeSpace   = typename KokkosTypes<Device>::ExeSpace;

  template <typename S>
  using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
  template <typename S>
  using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
  template <typename S>
  using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KokkosTypes<Device>::template view_1d_ptr_array<S, N>;

  //
  // ------- public API ---------
  //

  // Constructor, call from host
  //   size: The number of T's per sub-block
  //   max_used: The maximum number of active sub-blocks
  //   policy: The team policy for Kokkos kernels using this WorkspaceManager
  WorkspaceManager(int size, int max_used, TeamPolicy policy);

  // call from host.
  //
  // Will report usage statistics for your workspaces. These statistics will
  // have much more detail for debug builds.
  void report() const;

  class Workspace;

  // call from device
  //
  // Returns a Workspace object which provides access to sub-blocks.
  KOKKOS_INLINE_FUNCTION
  Workspace get_workspace(const MemberType& team) const;

  class Workspace {
   public:

    // Take an individual sub-block
    template <typename S=T>
    KOKKOS_INLINE_FUNCTION
    Unmanaged<view_1d<S> > take(const char* name) const;

    // Take several sub-blocks. The user gets pointers to their sub-blocks
    // via the ptrs argument.
    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many(const Kokkos::Array<const char*, N>& names,
                   const view_1d_ptr_array<S, N>& ptrs) const;

    // Similar to take_many except assumes that there is enough contiguous
    // memory avaiable in the Workspace for N sub-blocks. This method is higher-performing
    // than take_many.
    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many_contiguous_unsafe(const Kokkos::Array<const char*, N>& names,
                                     const view_1d_ptr_array<S, N>& ptrs) const;

    // Combines reset and take_many_contiguous_unsafe. This is the most-performant
    // option for a kernel to use N sub-blocks that are needed for the duration of the
    // kernel.
    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many_and_reset(const Kokkos::Array<const char*, N>& names,
                             const view_1d_ptr_array<S, N>& ptrs) const;

    // Release an individual sub-block.
    template <typename View>
    KOKKOS_FORCEINLINE_FUNCTION
    void release(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { release_impl<typename View::value_type>(space); }

#ifndef NDEBUG
    // Get the name of a sub-block
    template <typename View>
    KOKKOS_INLINE_FUNCTION
    const char* get_name(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { return get_name_impl<typename View::value_type>(space); }
#endif

    // Reset back to initial state. All sub-blocks will be considered inactive.
    KOKKOS_INLINE_FUNCTION
    void reset() const;

    // Print the linked list. Obviously not a device function.
    void print() const;

    //
    // ---------- Private --------------
    //

#ifndef KOKKOS_ENABLE_CUDA
   private:
#endif

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void release_impl(const Unmanaged<view_1d<S> >& space) const;

#ifndef NDEBUG
    template <typename S>
    KOKKOS_INLINE_FUNCTION
    const char* get_name_impl(const Unmanaged<view_1d<S> >& space) const;

    KOKKOS_INLINE_FUNCTION
    void change_num_used(int change_by) const;

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void change_indv_meta(const Unmanaged<view_1d<S> >& space, const char* name, bool release=false) const;

    KOKKOS_INLINE_FUNCTION
    int get_name_idx(const char* name, bool add) const;

    KOKKOS_INLINE_FUNCTION
    int get_alloc_count(const char* name) const
    { return m_parent.m_counts(m_ws_idx, get_name_idx(name), 0); }

    KOKKOS_INLINE_FUNCTION
    int get_release_count(const char* name) const
    { return m_parent.m_counts(m_ws_idx, get_name_idx(name), 1); }

    KOKKOS_INLINE_FUNCTION
    int get_num_used() const
    { return m_parent.m_num_used(m_ws_idx); }

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    bool is_active(const Unmanaged<view_1d<S> >& space) const
    { return m_parent.m_active(m_ws_idx, m_parent.template get_index<S>(space));}
#endif

    KOKKOS_INLINE_FUNCTION
    Workspace(const WorkspaceManager& parent, int ws_idx, const MemberType& team);

    friend struct unit_test::UnitWrap;
    friend class WorkspaceManager;

    const WorkspaceManager& m_parent;
    const MemberType& m_team;
    const int m_ws_idx;
    int& m_next_slot;
  }; // class Workspace

#ifndef KOKKOS_ENABLE_CUDA
 private:
#endif

  friend struct unit_test::UnitWrap;

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_index(const Unmanaged<view_1d<S> >& space) const
  { return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[0]; }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_next(const Unmanaged<view_1d<S> >& space) const
  { return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[1]; }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int set_next_and_get_index(const Unmanaged<view_1d<S> >& space, int next) const;

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  Unmanaged<view_1d<S> > get_space_in_slot(const int team_idx, const int slot) const;

  KOKKOS_INLINE_FUNCTION
  void init_metadata(const int ws_idx, const int slot) const;

  static void init(const WorkspaceManager& wm, const view_2d<T>& data,
                   const int concurrent_teams, const int max_used, const int total);

  //
  // data
  //

  enum { m_pad_factor   = OnGpu<ExeSpace>::value ? 1 : 32,
         m_max_name_len = 128,
         m_max_names    = 256
  };

  util::TeamUtils<ExeSpace> m_tu;
  int m_concurrent_teams, m_reserve, m_size, m_total, m_max_used;
#ifndef NDEBUG
  view_1d<int> m_num_used;
  view_1d<int> m_high_water;
  view_2d<bool> m_active;
  view_3d<char> m_curr_names;
  view_3d<char> m_all_names;
  view_3d<int> m_counts;
#endif
  view_1d<int> m_next_slot;
  view_2d<T> m_data;
}; // class WorkspaceManager

} // namespace util

#include "wsm_impl.hpp"

#endif
