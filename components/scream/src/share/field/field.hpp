#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "field_header.hpp"
#include <share/scream_types.hpp>
#include <share/util/scream_std_meta.hpp>
#include <share/util/scream_kokkos_utils.hpp>

#include <memory>   // For std::shared_ptr

namespace scream
{

// Tiny structure holding a few properties of the field allocation
template<typename ValueType>
struct FieldAllocProp {
  using value_type = ValueType;
  static constexpr int value_size = sizeof(value_type);

  FieldAllocProp () = default;

  template<typename SrcVT>
  FieldAllocProp (const FieldAllocProp<SrcVT>& src) {
    allocated  = src.allocated;
    alloc_size = src.alloc_size;
    if (std::is_same<SrcVT,value_type>::value) {
      pack_size  = src.pack_size;
    } else {
      int ratio = value_size>src.value_size ? value_size/src.value_size : src.value_size/value_size;
      pack_size = value_size>src.value_size ? ratio*src.pack_size : src.pack_size/ratio;
    }
  }
  FieldAllocProp& operator= (const FieldAllocProp&) = delete;

  int pack_size;
  int alloc_size;
  bool allocated;
};

// ======================== FIELD ======================== //

// A field should be composed of metadata info (the header) and a pointer to the view

template<typename DataType, typename MemSpace, typename MemManagement>
class Field {
public:
  using data_type           = DataType;
  using memory_space        = MemSpace;
  using memory_management   = MemManagement;
  using view_type           = ViewType<data_type,MemSpace,MemManagement>;
  using const_data_type     = typename view_type::traits::const_data_type;
  using non_const_data_type = typename view_type::traits::non_const_data_type;
  using value_type          = typename view_type::traits::value_type;
  using const_field_type    = Field<const_data_type,MemSpace,MemManagement>;
  using header_type         = FieldHeader;
  using alloc_prop_type     = FieldAllocProp<typename std::remove_const<value_type>::type>;
  using identifier_type     = header_type::identifier_type;

  // Constructor(s)
  Field () = delete;
  Field (const identifier_type& id);
  Field (const header_type& header);

  // WARNING: this is an expert user constructor. We do perform some checks, to make
  //          sure things are consistent, but it is not to be used frequently.
  Field (const std::shared_ptr<header_type>& fh, const alloc_prop_type& alloc_prop, const view_type& view);

  // This constructor allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field (const Field<SrcDT,MemSpace,MemManagement>& src);

  // Assignment: allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field& operator= (const Field<SrcDT,MemSpace,MemManagement>& src);

  // ---- Getters ---- //
  const header_type& get_header () const { return *m_header; }
        header_type& get_header ()       { return *m_header; }
  const view_type&   get_view   () const { return  m_view;   }

  const std::shared_ptr<header_type>& get_header_ptr () const { return m_header; }

  const alloc_prop_type& get_alloc_properties () const { return m_alloc_prop; }

  // ---- Setters ---- //

  // Set alloc properties
  void request_pack_size (const int pack_size);

  //Allocate the actual view
  void allocate_view ();

protected:

  // Metadata (name, rank, dims, customere/providers, time stamp, ...)
  std::shared_ptr<header_type>    m_header;

  // Actual data.
  view_type                       m_view;

  // Allocation properties
  alloc_prop_type                 m_alloc_prop;
};

// ================================= IMPLEMENTATION ================================== //

template<typename DataType, typename MemSpace, typename MemManagement>
Field<DataType,MemSpace,MemManagement>::
Field (const identifier_type& id)
 : m_header (new header_type(id))
{
  m_alloc_prop.allocated  = false;
  m_alloc_prop.alloc_size = 0;
  m_alloc_prop.pack_size  = 1;
}

template<typename DataType, typename MemSpace, typename MemManagement>
template<typename SrcDT>
Field<DataType,MemSpace,MemManagement>::
Field (const Field<SrcDT,MemSpace,MemManagement>& src)
 : m_header     (src.get_header_ptr())
 , m_view       (src.get_view())
 , m_alloc_prop (src.get_alloc_properties())
{

  // Check that underlying array layout is the same
  static_assert(std::is_same<non_const_data_type,
                             typename util::remove_all_consts<SrcDT>::type
                            >::value,
                "Error! Cannot use copy constructor if the underlying array layout is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<data_type,const_data_type>::value ||
                std::is_same<SrcDT,non_const_data_type>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
}

template<typename DataType, typename MemSpace, typename MemManagement>
Field<DataType,MemSpace,MemManagement>::
Field (const std::shared_ptr<header_type>& header, const alloc_prop_type& alloc_prop, const view_type& view)
 : m_header     (header)
 , m_view       (view)
 , m_alloc_prop (alloc_prop)
{
  // Perform some checks
  error::runtime_check(static_cast<bool>(header),"Error! Input header is null.\n");
  error::runtime_check(view.size()>0, "Error! Input view has 0 size.\n");

  error::runtime_check(static_cast<int>(view.size()*sizeof(value_type))==alloc_prop.alloc_size,
                       "Error! Input view seems to have the wrong allocation size:" +
                       std::to_string(view.size()*sizeof(value_type)) + " rather than " + std::to_string(alloc_prop.alloc_size) + ".\n");
  error::runtime_check(alloc_prop.allocated,
                       "Error! Alloc properties seem to store wrong information.\n");
  // if (view_type::Rank>1) {
    // if (util::is_pack<DataType>::value) {
      // error::runtime_check(view.extent_int(header->get_identifier().rank()-1) ==
                           // alloc_prop.pack_size * header->get_identifier().dims().back(),
                           // "Error! Input view seems to have been incorrectly allocated.\n");
    // } else {
    // }
  // }
}

template<typename DataType, typename MemSpace, typename MemManagement>
template<typename SrcDT>
Field<DataType,MemSpace,MemManagement>&
Field<DataType,MemSpace,MemManagement>::
operator= (const Field<SrcDT,MemSpace,MemManagement>& src) {
  // Check that underlying array layout is the same
  static_assert(std::is_same<non_const_data_type,
                             typename util::remove_all_consts<SrcDT>::type
                            >::value,
                "Error! Cannot use copy constructor if the underlying array layout is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<data_type,const_data_type>::value ||
                std::is_same<SrcDT,non_const_data_type>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
  if (&src!=*this) {
    m_header     = src.get_header_ptr();
    m_view       = src.get_view();
    m_alloc_prop = src.get_alloc_properties();
  }
}

template<typename DataType, typename MemSpace, typename MemManagement>
void Field<DataType,MemSpace,MemManagement>::request_pack_size (const int pack_size_request) {
  error::runtime_check(pack_size_request>0, "Error! Pack size must be positive.");
  // This check looks confusing, but it works (if input is >0). Write down the bits if you don't believe it.
  error::runtime_check((pack_size_request & (pack_size_request-1))==0, "Error! Pack size must be a power of two.\n");
  error::runtime_check(!m_alloc_prop.allocated, "Error! Cannot request a pack size after the view has been allocated.\n");

  // Set the pack size if the request is bigger than the one already stored
  m_alloc_prop.pack_size = std::max(pack_size_request,m_alloc_prop.pack_size);
}

template<typename DataType, typename MemSpace, typename MemManagement>
void Field<DataType,MemSpace,MemManagement>::allocate_view ()
{
  // Short name
  const auto& id = m_header->get_identifier();

  // Check the identifier has all the dimensions set
  error::runtime_check(id.are_dimensions_set(), "Error! Cannot create a field until all the field's dimensions are set.\n");

  // Determine allocation size,using identifier's dimensions and alloc_prop's pack_size
  const int last_dim = id.dims().back();
  const int ps = m_alloc_prop.pack_size;
  int last_dim_num_packs = last_dim;
  if (ps>1 && (last_dim % ps!=0)) {
    // Determine the size along the last dimension
    last_dim_num_packs = (last_dim + ps - 1) / ps;
  }

  // Create the layout. At the last position, we set the number of packs:
  // if DataType is Pack<T,N>, this will be smaller than last_dim, while if
  // DataType is not a pack, then last_dim_num_packs = last_dim.
  typename view_type::traits::array_layout layout;
  if (view_type::Rank==1) {
    // The field is stored as 1d
    layout.dimension[0] = (id.size() / last_dim) * last_dim_num_packs;
  } else {
    for (int idim=0; idim<id.rank()-1; ++idim) {
      layout.dimension[idim] = id.dim(idim);
    }
    layout.dimension[id.rank()-1] = last_dim_num_packs;
  }

  m_view = view_type(id.name(),layout);

  m_alloc_prop.alloc_size = m_view.size() * sizeof(value_type);

  m_alloc_prop.allocated = true;
}

} // namespace scream

#endif // SCREAM_FIELD_HPP
