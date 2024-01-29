#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "share/field/field_header.hpp"
#include "share/util/scream_combine_ops.hpp"
#include "share/scream_types.hpp"

#include "ekat/std_meta/ekat_std_type_traits.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include <memory>   // For std::shared_ptr
#include <string>

namespace scream
{

// Enum used when quering Field for a view on a specific mem space
enum HostOrDevice {
  Device = 0,
  Host
};

// ======================== FIELD ======================== //

// A field is composed of metadata info (the header) and a pointer to a view.
// Fields are always stored as 1D arrays of char-valued data. The associated
// view can be reshaped as needed to match a desired layout/datatype for a given client.
class Field {
public:

  // The syntax of std::enable_if is way too long...
  template<bool c, typename T, typename F>
  using cond_t = typename std::conditional<c,T,F>::type;

  template<bool c, typename T = void>
  using if_t = typename std::enable_if<c,T>::type;

  // Various kokkos-related types
  using device_t      = DefaultDevice;
  using host_device_t = HostDevice;
  using kt_dev        = KokkosTypes<device_t>;
  using kt_host       = KokkosTypes<host_device_t>;

  template<HostOrDevice HD>
  using get_device = cond_t<HD==Device, device_t,host_device_t>;

  // The data type of an N-dimensional array of T, with all dynamic extents
  template<typename T, int N>
  using data_nd_t = typename ekat::DataND<T,N>::type;

  // Types of device and host views given data type and memory traits
  template<typename DT, typename MT = Kokkos::MemoryManaged>
  using view_dev_t = typename kt_dev::template view<DT,MT>;
  template<typename DT, typename MT = Kokkos::MemoryManaged>
  using view_host_t = typename kt_host::template view<DT,MT>;

  // Analogue of the above, but with LayoutStride
  template<typename DT, typename MT = Kokkos::MemoryManaged>
  using strided_view_dev_t = typename kt_dev::template sview<DT,MT>;
  template<typename DT, typename MT = Kokkos::MemoryManaged>
  using strided_view_host_t = typename kt_host::template sview<DT,MT>;

private:
  // A bare DualView-like struct. This is an impl detail, so don't expose it.
  // NOTE: we could use DualView, but all we need is a container-like struct.
  template<typename DT, typename MT = Kokkos::MemoryManaged>
  struct dual_view_t {
    view_dev_t<DT,MT>   d_view;
    view_host_t<DT,MT>  h_view;

    template<HostOrDevice HD>
    const if_t<HD==Device,view_dev_t<DT,MT>>& get_view() const {
      return d_view;
    }
    template<HostOrDevice HD>
    const if_t<HD==Host,view_host_t<DT,MT>>& get_view() const {
      return h_view;
    }
  };
public:

  // Type of a view given data type, HostOrDevice enum, and memory traits
  template<typename DT, HostOrDevice HD, typename MT = Kokkos::MemoryManaged>
  using get_view_type = cond_t<HD==Device,view_dev_t<DT,MT>,view_host_t<DT,MT>>;

  template<typename DT, HostOrDevice HD, typename MT = Kokkos::MemoryManaged>
  using get_strided_view_type = cond_t<HD==Device,strided_view_dev_t<DT,MT>,strided_view_host_t<DT,MT>>;

  // Field stack classes types
  using header_type          = FieldHeader;
  using identifier_type      = FieldIdentifier;

  static constexpr int MaxRank = 6;

  // Constructor(s)
  Field () = default;
  explicit Field (const identifier_type& id);
  template<typename ViewT,
           typename = typename std::enable_if<Kokkos::is_view<ViewT>::value>::type>
  Field (const identifier_type& id,
         const ViewT& view_d);
  Field (const Field& src) = default;
  ~Field () = default;

  Field& operator= (const Field& src) = default;

  // ---- Getters and const methods---- //
  const header_type& get_header () const { return *m_header; }
        header_type& get_header ()       { return *m_header; }
  const std::shared_ptr<header_type>& get_header_ptr () const { return m_header; }

  // Returns a Field copy of this field, which cannot be modified
  Field get_const () const;

  bool is_read_only () const { return m_is_read_only; }

  // Creates a deep copy version of this field.
  // It is created with a pristine header (no providers/customers)
  Field clone () const;
  Field clone (const std::string& name) const;
  Field alias (const std::string& name) const;

  // Allows to get the underlying view, reshaped for a different data type.
  // The class will check that the requested data type is compatible with the
  // allocation. This allows each field to be stored as a 1d array, but then
  // be reshaped to the desired layout before being used.
  template<typename DT, HostOrDevice HD = Device>
  get_view_type<DT,HD>
  get_view () const;

  // Like the method above, but only for rank-1 fields, returning a view with LayoutStride.
  // This is safer to use for fields that could be a subfield of another one, since a
  // rank-1 view that is the subview of a 2d one along the 2nd index cannot have
  // LayoutRight, and must have LayoutStride instead.
  template<typename DT, HostOrDevice HD = Device>
  get_strided_view_type<DT,HD>
  get_strided_view () const;

  // These two getters are convenience function for commonly accessed metadata.
  // The same info can be extracted from the metadata stored in the FieldHeader
  DataType data_type () const { return get_header().get_identifier().data_type(); }
  const std::string& name () const { return get_header().get_identifier().name(); }
  int rank () const { return get_header().get_identifier().get_layout().rank(); }

  // WARNING: this is a power-user method. Its implementation, including assumptions
  //          on pre/post conditions, may change in the future. Use at your own risk!
  //          Read carefully the instructions below.
  // Allows to get a raw pointer (host or device) from the view stored in the field.
  // The user must provide the pointed type for the returned pointer. Such type must
  // either be char or the correct data type of this field.
  // Notice that the view stored may contain more data than just the
  // data of the current field. This can happen in two cases (possibly simultaneously).
  //   - The field was allocated in a way that allows packing. In this case,
  //     there may be padding along the last *physical* dimension of the field.
  //     In the stored 1d view, the padding may appear interleaved with actual
  //     data (due to the view layout being LayoutRight).
  //   - The field is a subfield of another field. In this case, this class
  //     actually stores the "parent" field view. So when calling this method,
  //     you will actually get the raw pointer of the parent field view.
  template<typename ST, HostOrDevice HD = Device>
  ST* get_internal_view_data () const {
    // Check that the scalar type is correct
    using nonconst_ST = typename std::remove_const<ST>::type;
    EKAT_REQUIRE_MSG ((field_valid_data_types().at<nonconst_ST>()==m_header->get_identifier().data_type()
                       or std::is_same<nonconst_ST,char>::value),
        "Error! Attempt to access raw field pointere with the wrong scalar type.\n");
    EKAT_REQUIRE_MSG (not m_is_read_only || std::is_const<ST>::value,
        "Error! Cannot get a non-const raw pointer to the field data if the field is read-only.\n");

    return reinterpret_cast<ST*>(get_view_impl<HD>().data());
  }

  // WARNING: this is a power-user method. Its implementation, including assumptions
  //          on pre/post conditions, may change in the future. Use at your own risk!
  //          Read carefully the instructions below.
  // Same as above, but does allows the field to be read-only. This is unsafe as one
  // could alter data that should not be changed. An example of where this is needed
  // is in SurfaceCoupling where exports need to access read-only fields via their
  // view.data ptr.
  template<typename ST, HostOrDevice HD = Device>
  ST* get_internal_view_data_unsafe () const {
    // Check that the scalar type is correct
    using nonconst_ST = typename std::remove_const<ST>::type;
    EKAT_REQUIRE_MSG ((field_valid_data_types().at<nonconst_ST>()==m_header->get_identifier().data_type()
                       or std::is_same<nonconst_ST,char>::value),
          "Error! Attempt to access raw field pointere with the wrong scalar type.\n");

    return reinterpret_cast<ST*>(get_view_impl<HD>().data());
  }

  // If someone needs the host view, some sync routines might be needed.
  // Note: this class takes no responsibility in keeping track of whether
  //       a sync is required in either direction. Mainly because we expect
  //       host views to be seldom used, and even less frequently modified.
  void sync_to_host () const;
  void sync_to_dev () const;

  // Set the field to a constant value (on host or device)
  template<typename T, HostOrDevice HD = Device>
  void deep_copy (const T value);

  // Copy the data from one field to this field
  template<HostOrDevice HD = Device>
  void deep_copy (const Field& src);

  // Updates this field y as y=alpha*x+beta*y
  // NOTE: ST=void is just so we can give a default to HD,
  //       but ST will *always* be deduced from input arguments.
  // NOTE: the type ST  must be such that no narrowing happens when
  //       casting the values to whatever the data type of this field is.
  //       E.g., if data_type()=IntType, you can't pass double's.
  template<HostOrDevice HD = Device, typename ST = void>
  void update (const Field& x, const ST alpha, const ST beta);

  // Special case of update with alpha=0
  template<HostOrDevice HD = Device, typename ST = void>
  void scale (const ST beta);

  // Scale a field y as y=y*x where x is also a field
  template<HostOrDevice HD = Device>
  void scale (const Field& x);

  // Returns a subview of this field, slicing at entry k along dimension idim
  // NOTES:
  //   - the output field stores *the same* 1d view as this field. In order
  //     to get the N-1 dimensional view, call get_view<DT>(), using
  //     the correct N-1 dimensional data type DT.
  //   - when calling get_view<DT>() on the N-1 dimensional subfield,
  //     we first get an N-dimensional view, then subview it at index k along
  //     dimension idim.
  //   - idim must be either 0 or 1. This is b/c we cannot subview an N-dim
  //     view along idim=2+ while keeping LayoutRight. Kokkos would force the
  //     resulting view to have layout stride, which would conflict with the
  //     return type of get_view<DT>().
  //   - If the field rank is 2, then idim cannot be 1. This is b/c Kokkos
  //     specializes view's traits for LayoutRight of rank 1, not allowing
  //     to store a stride for the slowest dimension.
  //   - If dynamic = true, it is possible to "reset" the slice index (k) at runtime.
  Field subfield (const std::string& sf_name, const ekat::units::Units& sf_units,
                       const int idim, const int index, const bool dynamic = false) const;
  Field subfield (const std::string& sf_name, const int idim,
                       const int index, const bool dynamic = false) const;
  Field subfield (const int idim, const int k, const bool dynamic = false) const;

  // If this field is a vector field, get a subfield for the ith component.
  // If dynamic = true, it is possible to "reset" the component index at runtime.
  // Note: throws if this is not a vector field.
  Field get_component (const int i, const bool dynamic = false);

  // Checks whether the underlying view has been already allocated.
  bool is_allocated () const { return m_data.d_view.data()!=nullptr; }

  // Whether this field is equivalent to the rhs. To be equivalent is
  // less strict than to have all the members equal. In particular,
  // this method returns true if and only if:
  //  - this==&rhs OR all the following apply
  //    - both views are allocated (if not, allocating one won't be reflected on the other)
  //    - both fields have the same header
  //    - both fields have the same views
  // We need to SFINAE on RhsRT, cause this==&rhs only works if the
  // two are the same. And we do want to check this==&rhs for the
  // same type, since if we didn't, f.equivalent(f) would return false
  // if f is not allocated...
  bool equivalent (const Field& rhs) const;

  // ---- Setters and non-const methods ---- //

  // Allocate the actual view
  void allocate_view ();

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  template<HostOrDevice HD, typename ST>
  void deep_copy_impl (const ST value);

  template<HostOrDevice HD, typename ST>
  void deep_copy_impl (const Field& src);

  template<CombineMode CM, HostOrDevice HD, typename ST>
  void update_impl (const Field& x, const ST alpha, const ST beta, const ST fill_val);

protected:

  template<HostOrDevice HD>
  const get_view_type<char*,HD>&
  get_view_impl () const {
    EKAT_REQUIRE_MSG (is_allocated (), "Error! View was not yet allocated.\n");
    return m_data.get_view<HD>();
  }

  // These SFINAE impl of get_subview are needed since subview_1 does not
  // exist for rank2 (or less) views.
  template<HostOrDevice HD, typename T, int N>
  if_t<(N>2),
       get_view_type<data_nd_t<T,N-1>,HD>>
  get_subview_1 (const get_view_type<data_nd_t<T,N>,HD>& v, const int k) const {
    return ekat::subview_1(v,k);
  }

  template<HostOrDevice HD, typename T, int N>
  if_t<(N<=2),
       get_view_type<data_nd_t<T,N-1>,HD>>
  get_subview_1 (const get_view_type<data_nd_t<T,N>,HD>&, const int) const {
    EKAT_ERROR_MSG ("Error! Cannot subview a rank2 view along the second dimension without losing LayoutRight.\n");
  }

  template<HostOrDevice HD,typename T,int N>
  auto get_ND_view () const
    -> if_t<N==MaxRank, get_view_type<data_nd_t<T,N>,HD>>;

  template<HostOrDevice HD,typename T,int N>
  auto get_ND_view () const
    -> if_t<(N<MaxRank), get_view_type<data_nd_t<T,N>,HD>>;

  // Metadata (name, rank, dims, customere/providers, time stamp, ...)
  std::shared_ptr<header_type>            m_header;

  // Actual data.
  dual_view_t<char*>    m_data;

  // Whether this field is read-only
  bool                  m_is_read_only = false;
};

// We use this to find a Field in a std container.
// We do NOT allow two entries with same field identifier in such containers.
inline bool operator== (const Field& lhs, const Field& rhs) {
  return lhs.get_header().get_identifier() == rhs.get_header().get_identifier();
}

} // namespace scream

// Include template methods implementation
#include "share/field/field_impl.hpp"

#endif // SCREAM_FIELD_HPP
