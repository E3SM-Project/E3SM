#ifndef EAMXX_FIELD_READER_HPP
#define EAMXX_FIELD_READER_HPP

#include "share/field/field.hpp"

#include <ekat_comm.hpp>

namespace scream
{

class FieldReader
{
public:
  // --- Constructor(s) & Destructor --- //
  FieldReader () = default;

  // Due to resource acquisition (in scorpio), avoid copies
  FieldReader (const FieldReader&) = delete;
  ~FieldReader ();

  // Due to resource acquisition (in scorpio), avoid copies
  FieldReader& operator= (const FieldReader&) = delete;

  // --- Methods --- //

  // Expose the ability to set/reset fields or filename for cases like data interpolation,
  // where we swap pointers but all the scorpio data structures are unchanged.
  void set_fields (const std::vector<Field>& fields);
  void set_file_specs(const std::string& filename, const std::string& io_type = "default");
  void set_file_specs(const std::string& filename, const std::map<FieldTag,std::string>& tag2name, const std::string& io_type = "default");
  void set_dim_decomp (const Field& gids, const ekat::Comm& comm);

  // To be called AFTER all of the above, to set up the data structures in scorpio
  void init_scorpio_structures ();

  // Read fields that were required via parameter list.
  void read (const int time_index = -1);

  // Cleans up the reader and closes scorpio stuff
  // NOTE: mostly useful for tests, when scorpio::finalize_session is in the same scope
  //       as the reader, and we MUST close all files before finalizing scorpio
  void clean_up ();

protected:

  // Once set to true, re-setting fields or filename will automatically
  // trigger a call to init_scorpio_structures();
  bool                m_inited = false;

  std::string         m_filename;

  // Entries in m_io_fields may alias entries in m_fields if there is no padding, and they are not subfields.
  std::vector<Field>  m_fields;
  std::vector<Field>  m_io_fields;

  // If the input file has non-standard names, we store them here, so we can alias the fields
  std::map<FieldTag, std::string> m_tag2name;

  Field               m_decomp_dim_gids;
  ekat::Comm          m_comm;
};

// Shorthand functions if you don't need to keep the reader around and the default io type is ok
void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const int time_index = -1);
void read_fields (const std::string& filename,
                  const std::initializer_list<Field>& fields,
                  const int time_index = -1);
void read_fields (const std::string& filename,
                  const std::map<std::string,Field>& fields,
                  const int time_index = -1);

void read_fields (const std::string& filename,
                  const std::vector<Field>& fields,
                  const Field& decomp_gids,
                  const ekat::Comm& comm,
                  const int time_index = -1);
void read_fields (const std::string& filename,
                  const std::initializer_list<Field>& fields,
                  const Field& decomp_gids,
                  const ekat::Comm& comm,
                  const int time_index = -1);
void read_fields (const std::string& filename,
                  const std::map<std::string,Field>& fields,
                  const Field& decomp_gids,
                  const ekat::Comm& comm,
                  const int time_index = -1);

} //namespace scream

#endif // EAMXX_FIELD_READER_HPP
