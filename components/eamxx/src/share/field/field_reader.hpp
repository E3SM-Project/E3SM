#ifndef EAMXX_FIELD_READER_HPP
#define EAMXX_FIELD_READER_HPP

#include "share/field/field.hpp"

#include <ekat_comm.hpp>

#include <string>
#include <map>

namespace scream
{

class FieldReader
{
public:
  // --- Constructor(s) & Destructor --- //
  FieldReader () = default;
  ~FieldReader ();

  // Due to resource acquisition (in scorpio), avoid copies
  FieldReader (const FieldReader&) = delete;
  FieldReader& operator= (const FieldReader&) = delete;

  // --- Methods --- //

  // Expose the ability to set/reset fields or filename for cases like data interpolation,
  // where we swap pointers but all the scorpio data structures are unchanged.
  void set_fields (const std::vector<Field>& fields);
  void set_file_specs(const std::string& filename, const std::string& io_type = "default");
  void set_file_specs(const std::string& filename, const std::map<std::string,std::string>& tag_rename, const std::string& io_type = "default");
  void set_dim_decomp (const Field& gids, const ekat::Comm& comm);

  // Read fields that were required via parameter list.
  void read (const int time_index = -1);

  // Cleans up the reader and closes scorpio stuff
  // NOTE: mostly useful for tests, when scorpio::finalize_session is in the same scope
  //       as the reader, and we MUST close all files before finalizing scorpio
  void clean_up ();

protected:

  // Flags to determine if anything needs to be re-inited before we read
  static constexpr int CLEAN      = 0;
  static constexpr int NEW_FILE   = 1;
  static constexpr int NEW_FIELDS = 2;
  static constexpr int NEW_DECOMP = 4;

  // Called lazily at the beginning of read if m_reader_state!=CLEAN
  void setup_internals ();

  std::string         m_filename;

  // Entries in m_io_fields may alias entries in m_fields if there is no padding, they are not subfields,
  // and the data type matches what's on the file.
  std::vector<Field>  m_fields;
  std::vector<Field>  m_io_fields;

  // Store info about the decomposed dim (if any), so if the file changes
  // we can replay it in the new file
  std::string          m_dim_decomp_name;
  std::vector<int64_t> m_dim_decomp_offsets;

  std::map<std::string,Field> m_layout_to_io_field;

  // If the input file has non-standard names, we store them here, so we can alias the fields
  std::map<std::string, std::string> m_tag_rename;

  // When read is called, if this is not READER_CLEAN, we must proceed to setup some things
  int m_reader_state = CLEAN;
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
