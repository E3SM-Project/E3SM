#ifndef PYGRID_HPP
#define PYGRID_HPP

#include "share/grid/mesh_free_grids_manager.hpp"
#include "pyunits.hpp"
#include "pyfield.hpp"

#include <pybind11/pybind11.h>

namespace scream {

struct PyGrid {
  std::shared_ptr<const AbstractGrid> grid;

  // 2d scalar, managed/unmanaged
  PyField scalar2d (const std::string& name,
                    const PyUnits& u)
  {
    FieldIdentifier fid(name,grid->get_2d_scalar_layout(),u.units,grid->name());
    return PyField(fid);
  }
  PyField scalar2d (const std::string& name,
                    const PyUnits& u,
                    pybind11::array_t<double> arr)
  {
    FieldIdentifier fid(name,grid->get_2d_scalar_layout(),u.units,grid->name());
    return PyField(fid,arr);
  }

  // 2d vector, managed/unmanaged
  PyField vector2d (const std::string& name,
                    const PyUnits& u,
                    int vdim)
  {
    FieldIdentifier fid(name,grid->get_2d_vector_layout(vdim),u.units,grid->name());
    return PyField(fid);
  }
  PyField vector2d (const std::string& name,
                    const PyUnits& u,
                    int vdim,
                    pybind11::array_t<double> arr)
  {
    FieldIdentifier fid(name,grid->get_2d_vector_layout(vdim),u.units,grid->name());
    return PyField(fid,arr);
  }

  // 3d midpoints scalar, managed/unmanaged
  PyField scalar3d_mid (const std::string& name,
                        const PyUnits& u)
  {
    FieldIdentifier fid(name,grid->get_3d_scalar_layout(true),u.units,grid->name());
    return PyField(fid);
  }
  PyField scalar3d_mid (const std::string& name,
                        const PyUnits& u,
                        pybind11::array_t<double> arr)
  {
    FieldIdentifier fid(name,grid->get_3d_scalar_layout(true),u.units,grid->name());
    return PyField(fid,arr);
  }

  // 3d interfaces scalar, managed/unmanaged
  PyField scalar3d_int (const std::string& name,
                        const PyUnits& u)
  {
    FieldIdentifier fid(name,grid->get_3d_scalar_layout(false),u.units,grid->name());
    return PyField(fid);
  }
  PyField scalar3d_int (const std::string& name,
                        const PyUnits& u,
                        pybind11::array_t<double> arr)
  {
    FieldIdentifier fid(name,grid->get_3d_scalar_layout(false),u.units,grid->name());
    return PyField(fid,arr);
  }
private:
  PyField create (const FieldIdentifier& fid,
                  pybind11::array_t<double> arr)
  {
    const auto rank = fid.get_layout().rank();
    PyField pyf;
    EKAT_REQUIRE_MSG (rank==arr.ndim(),
        "Error! Rank mismatch between input FieldIdentifier and pybind11::array_t.\n"
        "  - field name: " + fid.name() + "\n"
        "  - identifier rank: " + std::to_string(rank) + "\n"
        "  - array_t rank   : " + std::to_string(arr.ndim()) + "\n");
    switch (rank) {
      case 1:
      {
        Field::view_dev_t<double*> v(arr.mutable_data(0),arr.shape(0));
        pyf.f = Field(fid,v);
        break;
      }
      case 2:
      {
        Field::view_dev_t<double**> v(arr.mutable_data(0,0),arr.shape(0),arr.shape(1));
        pyf.f = Field(fid,v);
        break;
      }
      default:
        EKAT_ERROR_MSG ("AAARGH!\n");
    }
    return pyf;
  }
};

struct PyGridsManager {
  std::shared_ptr<MeshFreeGridsManager> gm;
  PyGrid grid;

  // Create field and allocate memory
  PyGridsManager(const std::string& name, int ncols, int nlevs)
  {
    ekat::ParameterList gm_pl;
    auto& pg_pl = gm_pl.sublist(name);
    pg_pl.set("type",std::string("point_grid"));
    pg_pl.set("number_of_global_columns",ncols);
    pg_pl.set("number_of_vertical_levels",nlevs);
    pg_pl.set<std::vector<std::string>>("aliases",{"Physics"});

    gm_pl.set("grids_names",std::vector<std::string>{name});

    ekat::Comm comm(MPI_COMM_WORLD);
    gm = std::make_shared<MeshFreeGridsManager>(comm,gm_pl);
    gm->build_grids();

    grid.grid = gm->get_grid(name);
  }

  const PyGrid& get_grid () const {
    return grid;
  }

  void cleanup () {
    gm = nullptr;
    grid = {};
  }
};

} // namespace scream

#endif // PYGRID_HPP
