#ifndef SCREAM_SC_DATA_MANAGER_HPP
#define SCREAM_SC_DATA_MANAGER_HPP

#include "share/eamxx_types.hpp"
#include "ekat/ekat_assert.hpp"

namespace scream {

// This struct provides a intermediary between the AD and the SurfaceCouplingImport/Export AtmophereProcess
// classes, allowing SCREAM to access CPL data pointers and info about the fields.
struct SCDataManager {

  template<typename DevT, typename DataT>
  using view_1d = typename KokkosTypes<DevT>::template view_1d<DataT>;
  template<typename DevT, typename DataT>
  using view_2d = typename KokkosTypes<DevT>::template view_2d<DataT>;

  using name_t = char[32];

  SCDataManager() = default;

  ~SCDataManager() = default;

  void setup_internals (const int num_cpl_fields, const int num_scream_fields, const int field_size,
                        Real* field_data_ptr,
#ifdef HAVE_MOAB
                        Real* field_data_moab_ptr,
#endif
                        char* field_names, int* field_cpl_indices_ptr,
                        int* field_vector_components_ptr, Real* field_constant_multiple_ptr,
                        bool* transfer_during_init_ptr)
  {
    m_num_cpl_fields    = num_cpl_fields;
    m_num_scream_fields = num_scream_fields;
    m_field_size        = field_size;

    // Data ptr is of size ncols,num_cpl_fields. All other views are of size num_scream_fields.
    EKAT_ASSERT_MSG(field_data_ptr              !=nullptr, "Error! Ptr for field data is null.");
    EKAT_ASSERT_MSG(field_names                 !=nullptr, "Error! Ptr for field names is null.");
    EKAT_ASSERT_MSG(field_cpl_indices_ptr       !=nullptr, "Error! Ptr for cpl indices is null.");
    EKAT_ASSERT_MSG(field_vector_components_ptr !=nullptr, "Error! Ptr for field vector components is null.");
    EKAT_ASSERT_MSG(field_constant_multiple_ptr !=nullptr, "Error! Ptr for constant multiple is null.");
    EKAT_ASSERT_MSG(transfer_during_init_ptr    !=nullptr, "Error! Ptr for initial transfer boolean is null.");
    m_field_data                 = decltype(m_field_data)                (field_data_ptr,              m_field_size, m_num_cpl_fields);
#ifdef HAVE_MOAB
    m_field_data_moab            = decltype(m_field_data_moab)           (field_data_moab_ptr,         m_num_cpl_fields, m_field_size);
#endif
    m_field_cpl_indices          = decltype(m_field_cpl_indices)         (field_cpl_indices_ptr,       m_num_scream_fields);
    m_field_vector_components    = decltype(m_field_vector_components)   (field_vector_components_ptr, m_num_scream_fields);
    m_field_constant_multiple    = decltype(m_field_constant_multiple)   (field_constant_multiple_ptr, m_num_scream_fields);
    m_field_transfer_during_init = decltype(m_field_transfer_during_init)(transfer_during_init_ptr,    m_num_scream_fields);

    // Fortran gives a 1d array of 32char strings. So let's memcpy the input char
    // strings into 2d char arrays. Each string is null-terminated (atm_mct_mod
    // makes sure of that).
    m_field_names = new name_t[num_scream_fields];
    std::memcpy(m_field_names, field_names, num_scream_fields*32*sizeof(char));
  }

  int get_field_size () const {
    return m_field_size;
  }

  int get_num_cpl_fields () const {
    return m_num_cpl_fields;
  }

  int get_num_scream_fields () const {
    return m_num_scream_fields;
  }

  Real* get_field_data_ptr () const {
    return m_field_data.data();
  }
#ifdef HAVE_MOAB
  Real* get_field_data_moab_ptr () const {
    return m_field_data_moab.data();
  }
#endif
  Real get_field_data_view_entry(const int i, const int f) {
    return m_field_data(i, f);
  }

  int* get_field_cpl_indices_ptr() const {
    return m_field_cpl_indices.data();
  }

  int get_cpl_index(const int f) {
    return m_field_cpl_indices(f);
  }

  char* get_field_name_ptr () const {
    return m_field_names[0];
  }

  std::string get_field_name(const int f) {
    return m_field_names[f];
  }

  int* get_field_vector_components_ptr () const {
    return m_field_vector_components.data();
  }

  Real* get_field_constant_multiple_ptr () const {
    return m_field_constant_multiple.data();
  }

  bool* get_field_transfer_during_init_ptr () const {
    return m_field_transfer_during_init.data();
  }

protected:

  int m_field_size;
  int m_num_cpl_fields;
  int m_num_scream_fields;

  view_2d<HostDevice, Real> m_field_data;
#ifdef HAVE_MOAB
  view_2d<HostDevice, Real> m_field_data_moab;
#endif
  name_t*                   m_field_names;
  view_1d<HostDevice, int>  m_field_cpl_indices;
  view_1d<HostDevice, int>  m_field_vector_components;
  view_1d<HostDevice, Real> m_field_constant_multiple;
  view_1d<HostDevice, bool> m_field_transfer_during_init;
};

} // scream

#endif // SCREAM_SC_DATA_MANAGER_HPP
