#include "eamxx_scorpio_types.hpp"

#include <ekat/ekat_assert.hpp>

namespace scream {
namespace scorpio {

std::string e2str (const FileMode mode)
{
  auto mode_int = static_cast<typename std::underlying_type<FileMode>::type>(mode);
  std::string s;
  switch (mode) {
    case Unset:  s = "UNSET";  break;
    case Read:   s = "READ";   break;
    case Write:  s = "WRITE";  break;
    case Append: s = "APPEND"; break;
    default:
      EKAT_ERROR_MSG (
          "Error! Unsupported/unrecognized FileMode value.\n"
          " - value: " + std::to_string(mode_int) + "\n");
  }
  return s;
}

IOType str2iotype(const std::string &str)
{
  if(str == "default"){
    return IOType::DefaultIOType;
  } else if(str == "netcdf") {
    return IOType::NetCDF;
  } else if(str == "pnetcdf") {
    return IOType::PnetCDF;
  } else if(str == "adios") {
    return IOType::Adios;
  } else if (str == "adiosc") {
    return IOType::Adiosc;
  } else if(str == "hdf5") {
    return IOType::Hdf5;
  } else {
    return IOType::Invalid;
  }
}

std::string iotype2str(const IOType iotype)
{
  std::string s;
  switch(iotype){
    case IOType::DefaultIOType: s = "default";  break;
    case IOType::NetCDF:        s = "netcdf";   break;
    case IOType::PnetCDF:       s = "pnetcdf";  break;
    case IOType::Adios:         s = "adios";    break;
    case IOType::Adiosc:        s = "adiosc";   break;
    case IOType::Hdf5:          s = "hdf5";     break;
    case IOType::Invalid:       s = "invalid";  break;
    default:
      EKAT_ERROR_MSG ("Unrecognized iotype.\n");
  }
  return s;
}

} // namespace scorpio
} // namespace scream
