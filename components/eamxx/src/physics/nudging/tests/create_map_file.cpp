#include <catch2/catch.hpp>

#include "share/io/scream_scorpio_interface.hpp"

TEST_CASE("create_map_file")
{
  using namespace scream;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::eam_init_pio_subsystem(comm);

  // Add a dof in the middle of two coarse dofs
  const int ngdofs_src = 12;
  const int ngdofs_tgt = 2*ngdofs_src-1;

  std::string filename = "map_ncol" + std::to_string(ngdofs_src)
                       + "_to_"     + std::to_string(ngdofs_src) + ".nc";

  // Existing dofs are "copied", added dofs are averaged from neighbors
  const int nnz = ngdofs_src + 2*(ngdofs_src-1);

  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::register_dimension(filename, "n_a", "n_a", ngdofs_src, false);
  scorpio::register_dimension(filename, "n_b", "n_b", ngdofs_tgt, false);
  scorpio::register_dimension(filename, "n_s", "n_s", nnz,        false);

  scorpio::register_variable(filename, "col", "col", "1", {"n_s"}, "int",    "int",    "");
  scorpio::register_variable(filename, "row", "row", "1", {"n_s"}, "int",    "int",    "");
  scorpio::register_variable(filename, "S",   "S",   "1", {"n_s"}, "double", "double", "");

  std::vector<scorpio::offset_t> dofs(nnz);
  std::iota(dofs.begin(),dofs.end(),0);
  scorpio::set_dof(filename,"col",dofs.size(),dofs.data());
  scorpio::set_dof(filename,"row",dofs.size(),dofs.data());
  scorpio::set_dof(filename,"S",  dofs.size(),dofs.data());

  scorpio::eam_pio_enddef(filename);

  std::vector<int> col(nnz), row(nnz);
  std::vector<double> S(nnz);
  for (int i=0; i<ngdofs_src; ++i) {
    col[i] = i;
    row[i] = i;
      S[i] = 1.0;
  }
  for (int i=0; i<ngdofs_src-1; ++i) {
    col[ngdofs_src+2*i] = i;
    row[ngdofs_src+2*i] = ngdofs_src+i;
      S[ngdofs_src+2*i] = 0.5;

    col[ngdofs_src+2*i+1] = i+1;
    row[ngdofs_src+2*i+1] = ngdofs_src+i;
      S[ngdofs_src+2*i+1] = 0.5;
  }

  scorpio::grid_write_data_array(filename,"row",row.data(),nnz);
  scorpio::grid_write_data_array(filename,"col",col.data(),nnz);
  scorpio::grid_write_data_array(filename,"S",  S.data(),  nnz);

  scorpio::eam_pio_closefile(filename);
  scorpio::eam_pio_finalize();
}
