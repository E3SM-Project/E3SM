#ifndef PYGRID_HPP
#define PYGRID_HPP

#include "share/grid/mesh_free_grids_manager.hpp"

#include "pyeamxx_ext.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <mpi.h>

namespace nb = nanobind;

namespace scream {

inline void create_grids_manager (int ncols, int nlevs, const std::string& latlon_nc_file)
{
  EKAT_REQUIRE_MSG (PySession::get().inited,
      "Error! You did not initialize pyeamxx, or you already finalized it!\n");
  auto& comm = PySession::get().comm;
  ekat::ParameterList gm_params;
  std::vector<std::string> grids_names = {"Physics"};
  auto& pl = gm_params.sublist("Physics");
  pl.set("type",std::string("point_grid"));
  pl.set("number_of_global_columns",ncols);
  pl.set("number_of_vertical_levels",nlevs);
  gm_params.set("grids_names",grids_names);
  gm_params.set("geo_data_source",std::string("CREATE_EMPTY_DATA"));

  if (latlon_nc_file!="") {
    gm_params.set("geo_data_source",std::string("IC_FILE"));
    gm_params.set("ic_filename",latlon_nc_file);
  }

  PySession::get().gm = create_mesh_free_grids_manager (comm, gm_params);
  PySession::get().gm->build_grids();
}
inline void create_grids_manager (int ncols, int nlevs)
{
  create_grids_manager(ncols,nlevs,"");
}

inline void nb_pygrid (nb::module_& m) {
  m.def("create_grids_manager",nb::overload_cast<int,int>(&create_grids_manager));
  m.def("create_grids_manager",nb::overload_cast<int,int,const std::string&>(&create_grids_manager));
}

} // namespace scream

#endif // PYGRID_HPP
