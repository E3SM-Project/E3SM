#include "catch2/catch.hpp"

#include "diagnostics/register_diagnostics.hpp"

#include "share/io/eamxx_io_utils.hpp"
#include "share/grid/point_grid.hpp"

namespace scream {

TEST_CASE("create_diag")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  register_diagnostics();

  // Create a grid
  const int ncols = 3*comm.size();
  const int nlevs = 10;
  auto grid = create_point_grid("Physics",ncols,nlevs,comm);

  SECTION ("field_at") {
    // FieldAtLevel
    auto d1 = create_diagnostic("BlaH_123_at_model_top",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(d1)!=nullptr);
    auto d2 = create_diagnostic("BlaH_123_at_model_bot",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(d2)!=nullptr);
    auto d3 = create_diagnostic("BlaH_123_at_lev_10",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(d3)!=nullptr);

    REQUIRE_THROWS(create_diagnostic("BlaH_123_at_modeltop",grid)); // misspelled

    // FieldAtPressureLevel
    auto d4 = create_diagnostic("BlaH_123_at_10mb",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtPressureLevel>(d4)!=nullptr);
    auto d5 = create_diagnostic("BlaH_123_at_10hPa",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtPressureLevel>(d5)!=nullptr);
    auto d6 = create_diagnostic("BlaH_123_at_10Pa",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtPressureLevel>(d6)!=nullptr);

    REQUIRE_THROWS(create_diagnostic("BlaH_123_at_400KPa",grid)); // invalid units

    // FieldAtHeight
    auto d7 = create_diagnostic("BlaH_123_at_10m_above_sealevel",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtHeight>(d7)!=nullptr);
    auto d8 = create_diagnostic("BlaH_123_at_10m_above_surface",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtHeight>(d8)!=nullptr);

    REQUIRE_THROWS(create_diagnostic("BlaH_123_at_10.5m",grid)); // missing _above_X
    REQUIRE_THROWS(create_diagnostic("BlaH_123_at_1km_above_sealevel",grid)); // invalid units
    REQUIRE_THROWS(create_diagnostic("BlaH_123_at_1m_above_the_surface",grid)); // invalid reference
  }

  SECTION ("precip_mass_flux") {
    auto d1 = create_diagnostic("precip_liq_surf_mass_flux",grid);
    REQUIRE (std::dynamic_pointer_cast<PrecipSurfMassFlux>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("precip_type")=="liq");

    auto d2 = create_diagnostic("precip_ice_surf_mass_flux",grid);
    REQUIRE (std::dynamic_pointer_cast<PrecipSurfMassFlux>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("precip_type")=="ice");

    auto d3 = create_diagnostic("precip_total_surf_mass_flux",grid);
    REQUIRE (std::dynamic_pointer_cast<PrecipSurfMassFlux>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("precip_type")=="total");
  }

  SECTION ("water_and_number_path") {
    auto d1 = create_diagnostic("LiqWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPathDiagnostic>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("Water Kind")=="Liq");

    auto d2 = create_diagnostic("IceWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPathDiagnostic>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("Water Kind")=="Ice");

    auto d3 = create_diagnostic("RainWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPathDiagnostic>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("Water Kind")=="Rain");

    auto d4 = create_diagnostic("RimeWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPathDiagnostic>(d4)!=nullptr);
    REQUIRE (d4->get_params().get<std::string>("Water Kind")=="Rime");

    auto d5 = create_diagnostic("VapWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPathDiagnostic>(d5)!=nullptr);
    REQUIRE (d5->get_params().get<std::string>("Water Kind")=="Vap");

    auto d6 = create_diagnostic("LiqNumberPath",grid);
    REQUIRE (std::dynamic_pointer_cast<NumberPathDiagnostic>(d6)!=nullptr);
    REQUIRE (d6->get_params().get<std::string>("Number Kind")=="Liq");

    auto d7 = create_diagnostic("IceNumberPath",grid);
    REQUIRE (std::dynamic_pointer_cast<NumberPathDiagnostic>(d7)!=nullptr);
    REQUIRE (d7->get_params().get<std::string>("Number Kind")=="Ice");

    auto d8 = create_diagnostic("RainNumberPath",grid);
    REQUIRE (std::dynamic_pointer_cast<NumberPathDiagnostic>(d8)!=nullptr);
    REQUIRE (d8->get_params().get<std::string>("Number Kind")=="Rain");
  }

  SECTION ("aerocom_cld") {
    auto d1 = create_diagnostic("AeroComCldTop",grid);
    REQUIRE (std::dynamic_pointer_cast<AeroComCld>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("AeroComCld Kind")=="Top");

    auto d2 = create_diagnostic("AeroComCldBot",grid);
    REQUIRE (std::dynamic_pointer_cast<AeroComCld>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("AeroComCld Kind")=="Bot");
  }

  SECTION ("vapor_flux") {
    auto d1 = create_diagnostic("MeridionalVapFlux",grid);
    REQUIRE (std::dynamic_pointer_cast<VaporFluxDiagnostic>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("Wind Component")=="Meridional");

    auto d2 = create_diagnostic("ZonalVapFlux",grid);
    REQUIRE (std::dynamic_pointer_cast<VaporFluxDiagnostic>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("Wind Component")=="Zonal");
  }

  SECTION ("atm_tend") {
    auto d1 = create_diagnostic("BlaH_123_atm_backtend",grid);
    REQUIRE (std::dynamic_pointer_cast<AtmBackTendDiag>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("Tendency Name")=="BlaH_123");
  }

  SECTION ("pot_temp") {
    auto d1 = create_diagnostic("LiqPotentialTemperature",grid);
    REQUIRE (std::dynamic_pointer_cast<PotentialTemperatureDiagnostic>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("Temperature Kind")=="Liq");

    auto d2 = create_diagnostic("PotentialTemperature",grid);
    REQUIRE (std::dynamic_pointer_cast<PotentialTemperatureDiagnostic>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("Temperature Kind")=="Tot");
  }

  SECTION ("vert_layer") {
    auto d1 = create_diagnostic("z_mid",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayerDiagnostic>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("vert_location")=="mid");
    auto d2 = create_diagnostic("z_int",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayerDiagnostic>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("vert_location")=="int");

    auto d3 = create_diagnostic("height_mid",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayerDiagnostic>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("vert_location")=="mid");
    auto d4 = create_diagnostic("height_int",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayerDiagnostic>(d4)!=nullptr);
    REQUIRE (d4->get_params().get<std::string>("vert_location")=="int");

    auto d5 = create_diagnostic("geopotential_mid",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayerDiagnostic>(d5)!=nullptr);
    REQUIRE (d5->get_params().get<std::string>("vert_location")=="mid");
    auto d6 = create_diagnostic("geopotential_int",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayerDiagnostic>(d6)!=nullptr);
    REQUIRE (d6->get_params().get<std::string>("vert_location")=="int");

    auto d7 = create_diagnostic("dz",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayerDiagnostic>(d7)!=nullptr);
    REQUIRE (d7->get_params().get<std::string>("vert_location")=="mid");
  }

}

} // namespace scream
