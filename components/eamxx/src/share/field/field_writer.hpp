#ifndef EAMXX_FIELD_WRITER_HPP
#define EAMXX_FIELD_WRITER_HPP

#include "share/field/field.hpp"

#include <ekat_comm.hpp>

namespace scream
{

class FieldWriter
{
public:
  FieldWriter () = default;

  FieldWriter (const FieldWriter&) = delete;
  ~FieldWriter ();

  FieldWriter& operator= (const FieldWriter&) = delete;

  void set_fields (const std::vector<Field>& fields);
  void set_file_specs(const std::string& filename, const std::string& io_type = "default");
  void set_file_specs(const std::string& filename, const std::map<FieldTag,std::string>& tag2name, const std::string& io_type = "default");
  void set_dim_decomp (const Field& gids, const ekat::Comm& comm);
  void set_time_dependent (const bool time_dep,
                           const std::string& time_units = "",
                           const std::string& time_name = "time");

  void init_scorpio_structures ();

  void write ();
  void write (const double time);

  void clean_up ();

protected:
  bool                m_inited = false;
  bool                m_time_dep = false;

  std::string         m_time_units;
  std::string         m_time_name = "time";
  std::string         m_filename;

  std::vector<Field>  m_fields;
  std::vector<Field>  m_io_fields;

  std::map<FieldTag, std::string> m_tag2name;

  Field               m_decomp_dim_gids;
  ekat::Comm          m_comm;

  void write_impl ();
};

void write_fields (const std::string& filename,
                   const std::vector<Field>& fields);
void write_fields (const std::string& filename,
                   const std::initializer_list<Field>& fields);
void write_fields (const std::string& filename,
                   const std::map<std::string,Field>& fields);

void write_fields (const std::string& filename,
                   const std::vector<Field>& fields,
                   const Field& decomp_gids,
                   const ekat::Comm& comm);
void write_fields (const std::string& filename,
                   const std::initializer_list<Field>& fields,
                   const Field& decomp_gids,
                   const ekat::Comm& comm);
void write_fields (const std::string& filename,
                   const std::map<std::string,Field>& fields,
                   const Field& decomp_gids,
                   const ekat::Comm& comm);

} //namespace scream

#endif // EAMXX_FIELD_WRITER_HPP
