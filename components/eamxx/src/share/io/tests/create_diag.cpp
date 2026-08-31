#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"

#include "share/io/eamxx_dexpr_diags.hpp"
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
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // Some diags look up geometry data as they are constructed
  grid->create_geometry_data("area",grid->get_2d_scalar_layout()).deep_copy(1.0);
  grid->create_geometry_data("lat", grid->get_2d_scalar_layout()).deep_copy(0.0);

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
    REQUIRE (std::dynamic_pointer_cast<WaterPath>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("water_kind")=="Liq");

    auto d2 = create_diagnostic("IceWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPath>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("water_kind")=="Ice");

    auto d3 = create_diagnostic("RainWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPath>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("water_kind")=="Rain");

    auto d4 = create_diagnostic("RimeWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPath>(d4)!=nullptr);
    REQUIRE (d4->get_params().get<std::string>("water_kind")=="Rime");

    auto d5 = create_diagnostic("VapWaterPath",grid);
    REQUIRE (std::dynamic_pointer_cast<WaterPath>(d5)!=nullptr);
    REQUIRE (d5->get_params().get<std::string>("water_kind")=="Vap");

    auto d6 = create_diagnostic("LiqNumberPath",grid);
    REQUIRE (std::dynamic_pointer_cast<NumberPath>(d6)!=nullptr);
    REQUIRE (d6->get_params().get<std::string>("number_kind")=="Liq");

    auto d7 = create_diagnostic("IceNumberPath",grid);
    REQUIRE (std::dynamic_pointer_cast<NumberPath>(d7)!=nullptr);
    REQUIRE (d7->get_params().get<std::string>("number_kind")=="Ice");

    auto d8 = create_diagnostic("RainNumberPath",grid);
    REQUIRE (std::dynamic_pointer_cast<NumberPath>(d8)!=nullptr);
    REQUIRE (d8->get_params().get<std::string>("number_kind")=="Rain");
  }

  SECTION ("vapor_flux") {
    auto d1 = create_diagnostic("MeridionalVapFlux",grid);
    REQUIRE (std::dynamic_pointer_cast<VaporFlux>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("wind_component")=="Meridional");

    auto d2 = create_diagnostic("ZonalVapFlux",grid);
    REQUIRE (std::dynamic_pointer_cast<VaporFlux>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("wind_component")=="Zonal");
  }

  SECTION ("atm_backtend") {
    // _atm_backtend is a built-in alias: X_atm_backtend → X_minus_X_prev_over_dt
    // The returned diagnostic is FieldOverDt with field_name = "BlaH_123_minus_BlaH_123_prev"
    auto d1 = create_diagnostic("BlaH_123_atm_backtend",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldOverDt>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("field_name")=="BlaH_123_minus_BlaH_123_prev");
    // The diag would name its field after the expansion, but the customer
    // asked for the alias. Without this, an output stream throws.
    REQUIRE (d1->get_params().get<std::string>("output_field_name")=="BlaH_123_atm_backtend");
  }

  SECTION ("field_prev") {
    auto d1 = create_diagnostic("BlaH_123_prev",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldPrev>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("field_name")=="BlaH_123");
  }

  SECTION ("field_over_dt") {
    auto d1 = create_diagnostic("BlaH_123_over_dt",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldOverDt>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("field_name")=="BlaH_123");
  }

  SECTION ("pot_temp") {
    auto d1 = create_diagnostic("LiqPotentialTemperature",grid);
    REQUIRE (std::dynamic_pointer_cast<PotentialTemperature>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("temperature_kind")=="Liq");

    auto d2 = create_diagnostic("PotentialTemperature",grid);
    REQUIRE (std::dynamic_pointer_cast<PotentialTemperature>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("temperature_kind")=="Tot");
  }

  SECTION ("vert_layer") {
    auto d1 = create_diagnostic("z_mid",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayer>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("vert_location")=="mid");
    auto d2 = create_diagnostic("z_int",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayer>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("vert_location")=="int");

    auto d3 = create_diagnostic("height_mid",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayer>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("vert_location")=="mid");
    auto d4 = create_diagnostic("height_int",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayer>(d4)!=nullptr);
    REQUIRE (d4->get_params().get<std::string>("vert_location")=="int");

    auto d5 = create_diagnostic("geopotential_mid",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayer>(d5)!=nullptr);
    REQUIRE (d5->get_params().get<std::string>("vert_location")=="mid");
    auto d6 = create_diagnostic("geopotential_int",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayer>(d6)!=nullptr);
    REQUIRE (d6->get_params().get<std::string>("vert_location")=="int");

    auto d7 = create_diagnostic("dz",grid);
    REQUIRE (std::dynamic_pointer_cast<VerticalLayer>(d7)!=nullptr);
    REQUIRE (d7->get_params().get<std::string>("vert_location")=="mid");
  }

  // ==========================================================================
  //  Expression syntax
  // ==========================================================================
  // Names no pattern above claims go to the dexpr front end, which drives the
  // same factory from the parsed expression. Operands go down by canonical
  // expression string and come back through create_diagnostic one at a time.

  SECTION ("dexpr_operators") {
    auto d1 = create_diagnostic("(qc+qv)*p_mid",grid);
    REQUIRE (std::dynamic_pointer_cast<BinaryOp>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("arg1")=="(qc+qv)");
    REQUIRE (d1->get_params().get<std::string>("arg2")=="p_mid");
    REQUIRE (d1->get_params().get<std::string>("binary_op")=="times");
    // The field is named after the request, not after "(qc+qv)_times_p_mid"
    REQUIRE (d1->get_params().get<std::string>("output_field_name")=="(qc+qv)*p_mid");

    // ...and the operand resolves on its own when it comes back around.
    auto inner = create_diagnostic("(qc+qv)",grid);
    REQUIRE (std::dynamic_pointer_cast<BinaryOp>(inner)!=nullptr);
    REQUIRE (inner->get_params().get<std::string>("arg1")=="qc");
    REQUIRE (inner->get_params().get<std::string>("arg2")=="qv");
    REQUIRE (inner->get_params().get<std::string>("binary_op")=="plus");

    auto d2 = create_diagnostic("T_mid-p_mid",grid);
    REQUIRE (d2->get_params().get<std::string>("binary_op")=="minus");
    auto d3 = create_diagnostic("T_mid/p_mid",grid);
    REQUIRE (d3->get_params().get<std::string>("binary_op")=="over");
  }

  SECTION ("dexpr_operator_precedence") {
    // 'qc_plus_qv_times_p_mid' groups by regex greediness: (qc+qv)*p_mid. As
    // an expression, '*' binds tighter.
    auto d1 = create_diagnostic("qc+qv*p_mid",grid);
    REQUIRE (d1->get_params().get<std::string>("binary_op")=="plus");
    REQUIRE (d1->get_params().get<std::string>("arg1")=="qc");
    REQUIRE (d1->get_params().get<std::string>("arg2")=="(qv*p_mid)");

    auto legacy = create_diagnostic("qc_plus_qv_times_p_mid",grid);
    REQUIRE (legacy->get_params().get<std::string>("binary_op")=="times");
    REQUIRE (legacy->get_params().get<std::string>("arg1")=="qc_plus_qv");
  }

  SECTION ("dexpr_numeric_literals") {
    // The underscore syntax cannot write a literal, so BinaryOp only ever saw
    // named constants. p_mid/100 must not look for a field called "100".
    auto d1 = create_diagnostic("p_mid/100",grid);
    REQUIRE (std::dynamic_pointer_cast<BinaryOp>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("arg2")=="100");
    // 100 is a scalar, so it is not requested as an input field
    REQUIRE (d1->get_input_fields_names().size()==1);
    REQUIRE (d1->get_input_fields_names().front()=="p_mid");

    auto d2 = create_diagnostic("T_mid*0.5",grid);
    REQUIRE (d2->get_params().get<std::string>("arg2")=="0.5");
    REQUIRE (d2->get_input_fields_names().size()==1);
  }

  SECTION ("dexpr_composition_is_unambiguous") {
    // 'qc.where(p_mid>100) + qv' has no unambiguous underscore spelling:
    // 'qc_where_p_mid_gt_100_plus_qv' hits the conditional-sampling pattern
    // first with a greedy field name, silently giving
    // ConditionalSampling(qc, p_mid, gt, 100_plus_qv). Canonical expression
    // strings cannot be mis-bound: no legacy pattern matches parentheses.
    auto d1 = create_diagnostic("qc.where(p_mid>100) + qv",grid);
    REQUIRE (std::dynamic_pointer_cast<BinaryOp>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("arg1")=="qc.where((p_mid>100))");
    REQUIRE (d1->get_params().get<std::string>("arg2")=="qv");

    auto lhs = create_diagnostic("qc.where((p_mid>100))",grid);
    REQUIRE (std::dynamic_pointer_cast<ConditionalSampling>(lhs)!=nullptr);
    REQUIRE (lhs->get_params().get<std::string>("field_name")=="qc");
    REQUIRE (lhs->get_params().get<std::string>("condition_lhs")=="p_mid");
    REQUIRE (lhs->get_params().get<std::string>("condition_cmp")=="gt");
    REQUIRE (lhs->get_params().get<std::string>("condition_rhs")=="100");

    // Two levels of nesting
    auto d2 = create_diagnostic("(qc.where(p_mid>100) + qv) * T_mid",grid);
    REQUIRE (d2->get_params().get<std::string>("arg1")=="(qc.where((p_mid>100))+qv)");
    REQUIRE (d2->get_params().get<std::string>("binary_op")=="times");
  }

  SECTION ("dexpr_isel") {
    auto d1 = create_diagnostic("T_mid.isel(lev=3)",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("field_name")=="T_mid");
    REQUIRE (d1->get_params().get<std::string>("vertical_location")=="lev_3");

    // xarray negative indexing; -1 is the one we can resolve without the
    // level count.
    auto d2 = create_diagnostic("T_mid.isel(lev=-1)",grid);
    REQUIRE (d2->get_params().get<std::string>("vertical_location")=="model_bot");
    REQUIRE_THROWS (create_diagnostic("T_mid.isel(lev=-2)",grid));

    // lev=0 is the model top; there is no separate label for it
    auto d3 = create_diagnostic("T_mid.isel(lev=0)",grid);
    REQUIRE (d3->get_params().get<std::string>("vertical_location")=="lev_0");
    REQUIRE_THROWS (create_diagnostic("T_mid.isel(lev='model_top')",grid));
  }

  SECTION ("dexpr_interp") {
    auto d1 = create_diagnostic("T_mid.interp(p_mid=500,units='hPa')",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtPressureLevel>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("pressure_value")=="500");
    REQUIRE (d1->get_params().get<std::string>("pressure_units")=="hPa");

    // Pa is the default, matching the units everything else in EAMxx uses
    auto d2 = create_diagnostic("T_mid.interp(p_mid=50000)",grid);
    REQUIRE (d2->get_params().get<std::string>("pressure_units")=="Pa");

    auto d3 = create_diagnostic("T_mid.interp(z_mid=10,reference='surface')",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldAtHeight>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("height_value")=="10");
    REQUIRE (d3->get_params().get<std::string>("height_units")=="m");
    REQUIRE (d3->get_params().get<std::string>("surface_reference")=="surface");

    REQUIRE_THROWS (create_diagnostic("T_mid.interp(p_mid=500,z_mid=10)",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.interp()",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.interp(p_mid=500,reference='surface')",grid));
  }

  SECTION ("dexpr_reductions") {
    auto d1 = create_diagnostic("T_mid.mean('col')",grid);
    REQUIRE (std::dynamic_pointer_cast<HorizAvg>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("field_name")=="T_mid");

    auto d2 = create_diagnostic("T_mid.mean('lev')",grid);
    REQUIRE (std::dynamic_pointer_cast<VertContract>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("contract_method")=="avg");

    auto d3 = create_diagnostic("T_mid.sum('lev')",grid);
    REQUIRE (std::dynamic_pointer_cast<VertContract>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("contract_method")=="sum");

    auto d4 = create_diagnostic("T_mid.mean('lev',weights='dp')",grid);
    REQUIRE (d4->get_params().get<std::string>("weighting_method")=="dp");
    // ...same as the legacy spelling
    auto legacy = create_diagnostic("T_mid_vert_avg_dp_weighted",grid);
    REQUIRE (legacy->get_params().get<std::string>("weighting_method")=="dp");

    REQUIRE_THROWS (create_diagnostic("T_mid.mean('time')",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.sum('col')",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.mean('col',weights='dp')",grid));
  }

  SECTION ("dexpr_where") {
    auto d1 = create_diagnostic("T_mid.where(qv>0.01)",grid);
    REQUIRE (std::dynamic_pointer_cast<ConditionalSampling>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("condition_rhs")=="0.01");
    REQUIRE (d1->get_params().get<std::string>("condition_cmp")=="gt");

    // The special operand names carry over from the legacy spelling
    auto d2 = create_diagnostic("mask.where(lev>5)",grid);
    REQUIRE (d2->get_params().get<std::string>("field_name")=="mask");
    REQUIRE (d2->get_params().get<std::string>("condition_lhs")=="lev");

    auto d3 = create_diagnostic("T_mid.where(qv<=0.01)",grid);
    REQUIRE (d3->get_params().get<std::string>("condition_cmp")=="le");
    auto d4 = create_diagnostic("T_mid.where(qv!=0)",grid);
    REQUIRE (d4->get_params().get<std::string>("condition_cmp")=="ne");

    // ConditionalSampling takes exactly one comparison
    REQUIRE_THROWS (create_diagnostic("T_mid.where(qv>0 and qc<1)",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.where(qv)",grid));
    // ...and a comparison is meaningless anywhere else
    REQUIRE_THROWS (create_diagnostic("T_mid>0",grid));
  }

  SECTION ("dexpr_other_calls") {
    auto d1 = create_diagnostic("T_mid.differentiate('p_mid')",grid);
    REQUIRE (std::dynamic_pointer_cast<VertDerivative>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("derivative_method")=="p");
    REQUIRE_THROWS (create_diagnostic("T_mid.differentiate('t')",grid));
    // the bare shorthands are gone; name the coordinate
    REQUIRE_THROWS (create_diagnostic("T_mid.differentiate('p')",grid));

    auto d2 = create_diagnostic("T_mid.histogram(bins=[0,1.5,2])",grid);
    REQUIRE (std::dynamic_pointer_cast<Histogram>(d2)!=nullptr);
    REQUIRE (d2->get_params().get<std::string>("bin_configuration")=="0_1.5_2");
    // Edges render in the shortest form that reads back, so an exponent in
    // the request need not survive into the configuration...
    auto d2b = create_diagnostic("T_mid.histogram(bins=[0,1e3])",grid);
    REQUIRE (d2b->get_params().get<std::string>("bin_configuration")=="0_1000.0");
    // ...but one that cannot be written without an exponent must be rejected,
    // since the diag re-splits the configuration on '_' and parses each piece.
    REQUIRE_THROWS (create_diagnostic("T_mid.histogram(bins=[0,1e30])",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.histogram(bins=[0])",grid));

    auto d3 = create_diagnostic("T_mid.zonal_mean(bins=20)",grid);
    REQUIRE (std::dynamic_pointer_cast<ZonalAvg>(d3)!=nullptr);
    REQUIRE (d3->get_params().get<std::string>("number_of_zonal_bins")=="20");
    REQUIRE_THROWS (create_diagnostic("T_mid.zonal_mean(bins=-1)",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.zonal_mean(bins='many')",grid));

    auto d4 = create_diagnostic("T_mid.shift(time=1)",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldPrev>(d4)!=nullptr);
    REQUIRE (d4->get_params().get<std::string>("field_name")=="T_mid");
    // A positive shift looks back, as in xarray. -1 would look ahead, which we
    // cannot do, and is the mistake isel(lev=-1) invites.
    REQUIRE_THROWS (create_diagnostic("T_mid.shift(time=-1)",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.shift(time=0)",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.shift(time=2)",grid));
    REQUIRE_THROWS (create_diagnostic("T_mid.prev()",grid));

    auto d5 = create_diagnostic("T_mid.over_dt()",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldOverDt>(d5)!=nullptr);
  }

  SECTION ("dexpr_tend") {
    // Shorthand for (X-X.shift(time=1)).over_dt(): the _atm_backtend tree.
    auto d1 = create_diagnostic("T_mid.tend()",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldOverDt>(d1)!=nullptr);
    REQUIRE (d1->get_params().get<std::string>("field_name")=="(T_mid-T_mid.shift((time=1)))");

    // The expansion is the canonical rendering, so writing it out by hand
    // resolves to the same diag.
    auto diff = create_diagnostic("(T_mid-T_mid.shift(time=1))",grid);
    REQUIRE (std::dynamic_pointer_cast<BinaryOp>(diff)!=nullptr);
    REQUIRE (diff->get_params().get<std::string>("binary_op")=="minus");
    REQUIRE (diff->get_params().get<std::string>("arg2")=="T_mid.shift((time=1))");

    auto prev = create_diagnostic("T_mid.shift(time=1)",grid);
    REQUIRE (std::dynamic_pointer_cast<FieldPrev>(prev)!=nullptr);
  }

  SECTION ("dexpr_rejects") {
    // Unknown function: the message names what IS available
    REQUIRE_THROWS (create_diagnostic("T_mid.iselect(lev=1)",grid));
    // Unknown keyword
    REQUIRE_THROWS (create_diagnostic("T_mid.isel(levl=1)",grid));
    // Free-standing spelling of a method
    REQUIRE_THROWS (create_diagnostic("mean(T_mid,'lev')",grid));
    // Operators with no diagnostic behind them
    REQUIRE_THROWS (create_diagnostic("T_mid**2",grid));
    REQUIRE_THROWS (create_diagnostic("-T_mid",grid));
    // Attribute access with no call
    REQUIRE_THROWS (create_diagnostic("T_mid.prev",grid));
    // An expression must produce a field
    REQUIRE_THROWS (create_diagnostic("[1,2]",grid));
    // A non-name on the lhs of '=' is not a keyword argument, so it counts
    // against the positional arity rather than dereferencing a null name
    REQUIRE_THROWS (create_diagnostic("T_mid.isel(1=2)",grid));
  }

  SECTION ("dexpr_leaves_plain_names_alone") {
    // A bare identifier is a diag class name or a typo; the factory decides.
    auto d1 = create_diagnostic("AtmosphereDensity",grid);
    REQUIRE (std::dynamic_pointer_cast<AtmDensity>(d1)!=nullptr);
    REQUIRE_FALSE (d1->get_params().isParameter("from_expression"));

    REQUIRE_THROWS (create_diagnostic("NotADiagnostic",grid));

    // Legacy names matched a pattern long before the expression front end.
    auto d2 = create_diagnostic("BlaH_123_where_qv_gt_0.01",grid);
    REQUIRE (std::dynamic_pointer_cast<ConditionalSampling>(d2)!=nullptr);
    REQUIRE_FALSE (d2->get_params().isParameter("from_expression"));
  }

  SECTION ("dexpr_every_registered_function_is_buildable") {
    // validate_registry() already proved each example matches its spec; this
    // proves the translator turns each one into a diagnostic. A function
    // registered with no case is otherwise caught only when someone writes it.
    const auto examples = dexpr_diagnostic_examples();
    REQUIRE (examples.size()>0);
    for (const auto& e : examples) {
      INFO ("example: " + e);
      REQUIRE (create_diagnostic(e,grid)!=nullptr);
    }
  }

  SECTION ("dexpr_marks_what_it_built") {
    // Customers must tell an expression from a plain name, since an
    // expression needs a ':=' output name.
    auto d1 = create_diagnostic("(qc+qv)*p_mid",grid);
    REQUIRE (d1->get_params().isParameter("from_expression"));
    REQUIRE (d1->get_params().get<bool>("from_expression"));
  }

}

} // namespace scream
