#include <catch2/catch.hpp>

#include "scream_config.h"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/io.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"

#include "ekat/ekat_pack.hpp"

namespace {
/* ======= Internal Functions needed for test ======*/

void set_spatial_vectors(Real* x, Real* y, Real* z, std::vector<int> dims3d);
void update_index_output(Kokkos::View<Int*> index_1d, Kokkos::View<Int**> index_2d, Kokkos::View<Int***> index_3d, std::vector<Int> dims3d, Int tt);
void update_data_output(Kokkos::View<Real*> data_1d, Kokkos::View<Real**> data_2d, Kokkos::View<Real***> data_3d, Real* x, Real* y, Real* z, const std::vector<int> dims3d, const Real t);
 
/* ================================================================================================================ */
TEST_CASE("scorpio_yaml", "") {

  using namespace scream;
  using namespace ekat::units;
  using Device = DefaultDevice;

/* The first step is to establish a Field Manager Repo to work with.  This example is fashioned from the 'field_repo'
 * test from /share/tests/field_tests.hpp                                                                          */ 
  std::vector<FieldTag> tag1d = {FieldTag::Column};
  std::vector<FieldTag> tag2d = {FieldTag::Column, FieldTag::VerticalLevel};
  std::vector<FieldTag> tag3d = {FieldTag::Column, FieldTag::VerticalLevel, FieldTag::Component};

  FieldIdentifier fid_x("x", tag1d, m);
  FieldIdentifier fid_y("y", tag1d, m);
  FieldIdentifier fid_z("z", tag1d, m);
  FieldIdentifier fid_index_1d("index_1d", tag1d, kg/s);
  FieldIdentifier fid_index_2d("index_2d", tag2d, kg/s);
  FieldIdentifier fid_index_3d("index_3d", tag3d, kg/s);
  FieldIdentifier fid_data_1d("data_1d", tag1d, m/s);
  FieldIdentifier fid_data_2d("data_2d", tag2d, m/s);
  FieldIdentifier fid_data_3d("data_3d", tag3d, m/s);

  std::vector<int> dimsx = {10};
  std::vector<int> dimsy = {5};
  std::vector<int> dimsz = {2};
  std::vector<int> dims1d = {dimsx[0]};
  std::vector<int> dims2d = {dimsx[0], dimsy[0]};
  std::vector<int> dims3d = {dimsx[0], dimsy[0], dimsz[0]};

  fid_x.set_dimensions(dimsx);
  fid_y.set_dimensions(dimsy);
  fid_z.set_dimensions(dimsz);
  fid_index_1d.set_dimensions(dims1d);
  fid_index_2d.set_dimensions(dims2d);
  fid_index_3d.set_dimensions(dims3d);
  fid_data_1d.set_dimensions(dims1d);
  fid_data_2d.set_dimensions(dims2d);
  fid_data_3d.set_dimensions(dims3d);

  // Create fields
  Field<Real,Device> scalar_x (fid_x);
  Field<Real,Device> scalar_y (fid_y);
  Field<Real,Device> scalar_z (fid_z);

  Field<Int,Device> scalar_index_1d (fid_index_1d);
  Field<Int,Device> scalar_index_2d (fid_index_2d);
  Field<Int,Device> scalar_index_3d (fid_index_3d);

  Field<Real,Device> scalar_data_1d (fid_data_1d);
  Field<Real,Device> scalar_data_2d (fid_data_2d);
  Field<Real,Device> scalar_data_3d (fid_data_3d);

  scalar_x.allocate_view();
  scalar_y.allocate_view();
  scalar_z.allocate_view();

  scalar_index_1d.allocate_view();
  scalar_index_2d.allocate_view();
  scalar_index_3d.allocate_view();

  scalar_data_1d.allocate_view();
  scalar_data_2d.allocate_view();
  scalar_data_3d.allocate_view();

  // Initialize spatial vectors
  auto xd = scalar_x.get_view(); 
  auto yd = scalar_y.get_view();
  auto zd = scalar_z.get_view();
  auto xh = Kokkos::create_mirror_view( xd );
  auto yh = Kokkos::create_mirror_view( yd );
  auto zh = Kokkos::create_mirror_view( zd );
  set_spatial_vectors(xh.data(), yh.data(), zh.data(), dims3d);
  Kokkos::deep_copy(xd,xh);
  Kokkos::deep_copy(yd,yh);
  Kokkos::deep_copy(zd,zh);

  // Field Repo
  FieldRepository<Real,DefaultDevice>  repo;
  repo.registration_begins();
  repo.register_field(fid_x,"group_1");
  repo.register_field(fid_y,"group_1");
  repo.register_field(fid_z,"group_1");
  repo.register_field(fid_index_1d,"group_2");
  repo.register_field(fid_index_2d,"group_2");
  repo.register_field(fid_index_3d,"group_2");
  repo.register_field(fid_data_1d,"group_3");
  repo.register_field(fid_data_2d,"group_3");
  repo.register_field(fid_data_3d,"group_3");
  repo.registration_ends();

/* The next step is to initialize the set of IO class objects to handle output,
 * and sync it with the field manager we just created. */ 

  // Update index variables
  auto index_1d_dev = scalar_index_1d.get_view(); 
  auto index_2d_dev = scalar_index_2d.get_reshaped_view<Int**>();
  auto index_3d_dev = scalar_index_3d.get_reshaped_view<Int***>();
  auto index_1d_h = Kokkos::create_mirror_view( index_1d_dev );
  auto index_2d_h = Kokkos::create_mirror_view( index_2d_dev );
  auto index_3d_h = Kokkos::create_mirror_view( index_3d_dev );
  update_index_output(index_1d_dev, index_2d_dev, index_3d_dev, dims3d, 0);
  Kokkos::deep_copy(index_1d_dev,index_1d_h);
  Kokkos::deep_copy(index_2d_dev,index_2d_h);
  Kokkos::deep_copy(index_3d_dev,index_3d_h);
  // Update data variables
  auto data_1d_dev = scalar_data_1d.get_view(); 
  auto data_2d_dev = scalar_data_2d.get_reshaped_view<Real**>();
  auto data_3d_dev = scalar_data_3d.get_reshaped_view<Real***>();
  auto data_1d_h = Kokkos::create_mirror_view( data_1d_dev );
  auto data_2d_h = Kokkos::create_mirror_view( data_2d_dev );
  auto data_3d_h = Kokkos::create_mirror_view( data_3d_dev );
  update_data_output(data_1d_dev, data_2d_dev, data_3d_dev, xh.data(), yh.data(), zh.data(), dims3d, 0);
  Kokkos::deep_copy(data_1d_dev,data_1d_h);
  Kokkos::deep_copy(data_2d_dev,data_2d_h);
  Kokkos::deep_copy(data_3d_dev,data_3d_h);
} // TEST_CASE scorpio_yaml
/* ================================================================================================================ */

void set_spatial_vectors(Real* x, Real* y, Real* z, std::vector<int> dims3d) 
{
  Real pi = 2*acos(0.0);

  for (int ii=0;ii<dims3d[0];ii++)
  {
    x[ii] = 2.0*pi/dims3d[0]*(ii+1);
  }
  for (int jj=0;jj<dims3d[1];jj++) 
  {
    y[jj] = 4.0*pi/dims3d[1]*(jj+1);
  }
  for (int kk=0;kk<dims3d[2];kk++)
  {
    z[kk] = 100*(kk+1);
  }
}

void update_index_output(Kokkos::View<Int*> index_1d, Kokkos::View<Int**> index_2d, Kokkos::View<Int***> index_3d, const std::vector<int> dims3d, const Int tt)
{
  for (int ii=0;ii<dims3d[0];ii++)
  {
    index_1d(ii) = ii + 10000*tt;
    for (int jj=0;jj<dims3d[1];jj++) 
    {
      index_2d(ii,jj) = ii + jj*100 + 10000*tt;
      for (int kk=0;kk<dims3d[2];kk++)
      {
        index_3d(ii,jj,kk) = ii + jj*100 + kk*1000 + 10000*tt;
      }
    }
  }
}

void update_data_output(Kokkos::View<Real*> data_1d, Kokkos::View<Real**> data_2d, Kokkos::View<Real***> data_3d, Real* x, Real* y, Real* z, const std::vector<int> dims3d, const Real t)
{
  for (int ii=0;ii<dims3d[0];ii++)
  {
    data_1d(ii) = 0.1 * cos(x[ii]+t);
    for (int jj=0;jj<dims3d[1];jj++) 
    {
      data_2d(ii,jj) = sin(y[jj]+t) * data_1d(ii);
      for (int kk=0;kk<dims3d[2];kk++)
      {
        data_3d(ii,jj,kk) = z[kk] + data_2d(ii,jj);
      }
    }
  }
}

} // namespace
