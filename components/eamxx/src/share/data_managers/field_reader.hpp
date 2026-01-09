#ifndef SCREAM_FIELD_READER_HPP
#define SCREAM_FIELD_READER_HPP

#include "share/data_managers/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat_parameter_list.hpp>
#include <ekat_logger.hpp>

//  The FieldReader class handles reading Field from file

namespace scream
{

class FieldReader 
{
public:
  // --- Constructor(s) & Destructor --- //
  FieldReader () = default;
  FieldReader (const std::string& filename,
               const std::shared_ptr<const AbstractGrid>& grid,
               const std::vector<Field>& fields,
               const std::string& iotype = "default");

  // This constructor defers initialization
  // to when set_fields and reset_filename are called
  FieldReader (const std::shared_ptr<const AbstractGrid>& grid);

  // Due to resource acquisition (in scorpio), avoid copies
  FieldReader (const FieldReader&) = delete;
  ~FieldReader ();

  // Due to resource acquisition (in scorpio), avoid copies
  FieldReader& operator= (const FieldReader&) = delete;

  // --- Methods --- //
  void init (const ekat::ParameterList& params,
             const std::vector<const Field>& fields);

  // Read fields that were required via parameter list.
  void read_variables (const int time_index = -1);

  // Cleans up the class
  void finalize();

  // Getters
  std::string get_filename() { return m_filename; } // Simple getter to query the filename for this stream.

  // Expose the ability to set/reset fields for cases like data interpolation,
  // where we swap pointers but all the scorpio data structures are unchanged.
  void set_fields (const std::vector<Field>& fields);
  void reset_filename (const std::string& filename,
                       const std::string& iotype = "default");

  // Option to set a logger
  void set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger);
protected:

  void set_grid (const std::shared_ptr<const AbstractGrid>& grid);
  void init_scorpio_structures (const std::string& iotype);

  void set_decompositions();

  std::vector<std::string> get_vec_of_dims (const FieldLayout& layout);

  // Internal variables
  std::shared_ptr<FieldManager>         m_fields_from_user;
  std::shared_ptr<FieldManager>         m_fields_for_scorpio;
  std::shared_ptr<const AbstractGrid>   m_io_grid;

  std::vector<std::string>  m_fields_names;
  std::string               m_filename;

  bool m_fields_inited  = false;
  bool m_scorpio_inited = false;

  std::shared_ptr<ekat::logger::LoggerBase> m_atm_logger = console_logger(ekat::logger::LogLevel::warn);
}; // Class FieldReader

} //namespace scream

#endif // SCREAM_FIELD_READER_HPP
