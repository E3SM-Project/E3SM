#ifndef SCREAM_SCORPIO_TYPES_HPP
#define SCREAM_SCORPIO_TYPES_HPP

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <cstdint>

namespace scream {
namespace scorpio {

// Enum denoting how we open a file
enum FileMode {
  Unset = 0,
  Read = 1,
  Write = 2,
  Append = Read | Write
};
std::string e2str (const FileMode mode);

// I/O types supported
enum IOType {
  // Default I/O type is used to let the code choose I/O type as needed (via CIME)
  DefaultIOType = 0,
  NetCDF,
  PnetCDF,
  Adios,
  Adiosc,
  Hdf5,
  Invalid
};

IOType str2iotype(const std::string &str);
std::string iotype2str(const IOType iotype);

// The type used by PIOc for offsets
using offset_t = std::int64_t;

/*
 * The following PIOxyz types each represent an entity that
 * is associated in scorpio with an id. For each of them, we add some
 * additional metadata, to avoid having to call PIOc_inq_xyz every time
 * we need such information.
 * While these types are defined in this public header, the customers
 * of scream_io will never be able to get any object of these types out
 * of the internal database. In particular, all data is stored in a
 * ScorpioSession singleton class, whose declaration is hidden inside
 * eamxx_scorpio_interface.cpp.
 */

// The basic common data of any PIO entity
struct PIOEntity {
  int ncid = -1;      // PIO access all data via their id
  std::string name;   // In EAMxx, we prefer to use a string name than a number id
};

// An entity that is linked to a specific file
// This allows an entity to retrieve its file
struct PIOFileEntity : public PIOEntity {
  int fid = -1;       // Allow retrieving file from the entity
};

// A dimension
struct PIODim : public PIOFileEntity {
  int length       = -1;
  bool unlimited   = false;

  // In case we decompose the dimension, this will store
  // the owned offsets on this rank
  // NOTE: use a pointer, so we can detect if a decomposition already
  //       existed or not when we set one.
  std::shared_ptr<std::vector<offset_t>> offsets;
};

// A decomposition
// NOTE: the offsets of a PIODecomp are not the same as the offsets of
//       stored in its dim. The latter are the offsets *along that dim*.
//       A PIODecomp is associated with a Nd layout, which includes decomp_dim
//       among its dimensions. PIODecomp::offsets are the offsets of the full
//       array layout owned by this rank. Hence, there can be many PIODecomp
//       all storing the same dim
struct PIODecomp : public PIOEntity {
  std::vector<offset_t>           offsets;  // Owned offsets
  std::shared_ptr<const PIODim>   dim; 
};

// A variable
struct PIOVar : public PIOFileEntity {
  // Note: if time_dep=true, we will add it to the list of dims passed
  // to scorpio, but the time dim will not appear in this list.
  std::vector<std::shared_ptr<const PIODim>>          dims;

  std::vector<std::string> dim_names () const {
    std::vector<std::string> n;
    for (auto d : dims) {
      n.push_back(d->name);
    }
    return n;
  }

  std::string dtype;
  std::string nc_dtype;
  std::string units;

  bool time_dep = false;
  
  // Extra safety measure: for time_dep vars, use this to check that
  // we are not writing more slices than the current time dim length.
  int num_records = 0;

  std::shared_ptr<const PIODecomp> decomp;

  // Used only if a) var is not decomposed, and b) dtype!=nc_dtype
  int size = -1; // Product of all dims
  std::vector<char> buf;
};

// A file, which is basically a container for dims and vars
struct PIOFile : public PIOEntity {
  std::map<std::string,std::shared_ptr<PIODim>>   dims;
  std::map<std::string,std::shared_ptr<PIOVar>>   vars;

  std::shared_ptr<PIODim> time_dim;
  FileMode mode;
  IOType   iotype;
  bool enddef = false;

  // We keep track of how many places are currently using this file, so that we
  // can close it only when they are all done.
  int num_customers = 0;
};

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
