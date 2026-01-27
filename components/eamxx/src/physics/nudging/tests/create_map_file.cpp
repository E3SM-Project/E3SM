#include <catch2/catch.hpp>

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

TEST_CASE("create_map_file")
{
  using namespace scream;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  // Add a dof in the middle of two coarse dofs
  const int ngdofs_src = 12;
  const int ngdofs_tgt = 2*ngdofs_src-1;

  std::string filename = "map_ncol" + std::to_string(ngdofs_src)
                       + "_to_"     + std::to_string(ngdofs_tgt) + ".nc";

  // Existing dofs are "copied", added dofs are averaged from neighbors
  const int nnz = ngdofs_src + 2*(ngdofs_src-1);

  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::define_dim(filename, "n_a", ngdofs_src);
  scorpio::define_dim(filename, "n_b", ngdofs_tgt);
  scorpio::define_dim(filename, "n_s", nnz);

  scorpio::define_var(filename, "col", {"n_s"}, "int");
  scorpio::define_var(filename, "row", {"n_s"}, "int");
  scorpio::define_var(filename, "S",   {"n_s"}, "double");

  scorpio::enddef(filename);

  std::vector<int> col(nnz), row(nnz);
  std::vector<double> S(nnz);
  const int gid_base = 1;
  for (int i=0; i<ngdofs_src; ++i) {
    col[i] = gid_base + i;
    row[i] = gid_base + i;
      S[i] = 1.0;
  }
  for (int i=0; i<ngdofs_src-1; ++i) {
    col[ngdofs_src+2*i] = gid_base + i;
    row[ngdofs_src+2*i] = gid_base + ngdofs_src+i;
      S[ngdofs_src+2*i] = 0.5;

    col[ngdofs_src+2*i+1] = gid_base + i+1;
    row[ngdofs_src+2*i+1] = gid_base + ngdofs_src+i;
      S[ngdofs_src+2*i+1] = 0.5;
  }

  scorpio::write_var(filename,"row",row.data());
  scorpio::write_var(filename,"col",col.data());
  scorpio::write_var(filename,"S",  S.data());

  scorpio::release_file(filename);
  scorpio::finalize_subsystem();
}
