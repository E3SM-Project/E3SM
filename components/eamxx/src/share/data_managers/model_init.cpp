#include "share/data_manager/default_model_init.hpp"

namespace scream
{

ModelInit::
ModelInit (const ekat::ParameterList& params)
 : m_params (params)
{
  // Nothing to do
}

void ModelInit::
run (const std::shared_ptr<FieldsManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  // Gather all fields that need to be inited
  auto g2f = gather_fields(fm,t0,run_type);

  // Read topography fields
  read_topography_fields(g2f);

  // Initial runs can also set fields to a constant
  if (run_type==RunType::Initial) {
    set_constant_fields(g2f);
  }

  // If we set ALL fields to a constant, we're done
  if (g2f.size()==0)
    return;

  auto filename = params.get<std::string>("filename");
  for (auto& [gn, fields] : grid2fields) {
    auto grid = m_gm->get_grid(gn);
    auto gids = grid->get_partitioned_dim_gids();
    auto& fields = g2f.at(gn);
    read_fields(filename,fields,gids,m_comm);
  }
}

auto ModelInit::
gather_fields (const std::shared_ptr<FieldsManager>& fm,
               const util::TimeStamp& t0,
               const RunType run_type)
 -> strmap_t<strmap_t<Field>>
{
  auto group_name = run_type==RunType::Restart ? "RESTART" : "STARTUP";

  strmap_t<strmap_t<Field>> g2f;
  auto gm = fm->get_grids_manager();
  for (auto gn : gm->get_grids_names()) {
    auto& n2f = g2f[gn];
    if (not fm->has_group(group_name,gn))
      continue;
    auto group = fm->get_group(group_name,gn);
    for (auto& [n,fptr] : group.m_individual_fields) {
      n2f[n] = *fptr;
      // Immediately set the time stamp. If we can't init a field later,
      // we error out anyways...
      fptr->get_header().get_tracking().update_time_stamp(t0);
    }
  }
}

void ModelInit::
read_topography_fields (strmap_t<strmap_t<Field>>& g2f)
{
  for (auto gn : gm->get_grids_names()) {
    auto& n2f = g2f[gn];
    std::vector<Field> topo_fields;
    if (n2f.count("phis")) {
      const auto& phis = n2f.at("phis");
      if (gn=="physics_pg2")
        topo_fields.push_back(phis);
      else if (gn=="physics_gll" or gn=="point_grid")
        topo_fields.push_back(phis.alias("PHIS_d");
      else
        EKAT_ERROR_MSG ("Error! Requesting phis on an unknown grid: " + gn + ".\n");
    }
    if (n2f.count("sgh30")) {
      // The eamxx field "sgh30" is called "SGH30" in the
      // topography file and is only available on the PG2 grid.
      EKAT_REQUIRE_MSG(gn == "physics_pg2",
        "Error! Requesting sgh30 field on " + gn +
        " topo file only has sgh30 for physics_pg2.\n");
      const auto& shg30 = n2f.at("sgh30");
      topo_fields.push_back(sgh30.alias("SGH30");
    }
    if (fm->has_field("sgh",gn)) {
      // The eamxx field "sgh" is called "SGH" in the
      // topography file and is only available on the PG2 grid.
      EKAT_REQUIRE_MSG(gn == "physics_pg2",
        "Error! Requesting sgh30 field on " + gn +
        " topo file only has sgh for physics_pg2.\n");
      const auto& shg = n2f.at("sgh");
      topo_fields.push_back(sgh.alias("SGH");
    }

    if (topo_fields.size()==0)
      continue;

    EKAT_REQUIRE_MSG(m_params.isParameter("topography_filename"),
      "Error! Topography data was requested in the FM, but no "
      "topography_filename or entry matching the field name "
      "was given in IC parameters.\n");

    auto topo_file = ic_pl.get<std::string>("topography_filename");
    auto grid = gm->get_grid(gn);

    FieldReader reader;
    if (grid_name=="physics_gll")
      reader.set_file_specs(file_name,{ {"ncol","ncol_d"} });
    else
      reader.set_file_specs(file_name);
    reader.set_dim_decomp(grid->get_partitioned_dim_gids(),m_comm);
    reader.set_fields(topo_fields);
    reader.read();

    for (auto& f : topo_fields) {
      n2f.erase(f.name()); // remove from list of fields to read
    }
  }
}

void ModelInit::
set_constant_fields (strmap_t<strmap_t<Field>>& g2f)
{
  // First, parse the param list entry and build name->value map
  const auto& names_and_vals = m_params.get<strvec_t>("constant_fields");
  strmap_t<double> name2val;
  for (auto s : names_and_vals) {
    auto tokens = ekat::split(ekat::strip(s,' '),"=");
    EKAT_REQUIRE_MSG (tokens.size()==2,
        "[ModelInit::set_constant_fields] Error! Badly formatted entry.\n"
        " - entry: '" + s + "'\n"
        " - expected format: 'name = VALUE'\n");

    // Attempt to convert tokens[1] to a double
    double val;
    auto [ptr, ec] = std::from_chars(tokens[1].data(), tokens[1].data() + tokens[1].size(), val);

    if (ec != std::errc{} or ptr != tokens[1].data() + tokens[1].size()) {
      // If the conversion failed, we check if tokens[1] is the name of a known physics constant
      auto& pc_dict = physics::Constants<double>::dictionary();
      EKAT_REQUIRE_MSG (pc_dict.count(tokens[1]==0),
          "[ModelInit::set_constant_fields] Error! Badly formatted entry.\n"
          " Entry '" + tokens[1] + "' is neither a double nor a known physics constant.\n");
      val = pc_dict.at(tokens[1]).value;
    }
    name2val[tokens[0]] = val;
  };

  // Loop over the Fields map and set fields to the requested value.
  // Keep track of what was requested/found, so we can error out if a name is invalid
  strmap<bool> found;
  std::set<std::string> valid;
  for (auto it : g2f) {
    auto& n2f = g2f.second;
    for (auto it2=n2f.begin(); it2!=n2f.end(); ) {
      valid.insert(it2.first);
      auto it3 = name2val.find(it2.first);
      if (it3==name2val.end()) {
        ++it2;  // No value prescribed for this field
      } else {
        found[it3.first] = true;
        it2.second.deep_copy(it2.second);
        it2 = n2f.erase(it2);
      }
    }
  }

  EKAT_REQUIRE_MSG (found.size()==name2val.size(),
      "[ModelInit::set_constant_fields] Error! One (or more) of the specified constant fields were not found.\n"
      " - requested constant fields: " + ekat::join(name2val,[](auto it){return it.first},",") + "\n"
      " - atm input fields: " + ekat::join(valid,",") + "\n"); 
}

void ModelInit::
set_constant_fields (std::vector<Field>& fields)
{
  // Helper lambda that splits "foo=123" into a string and value pair
  auto key_val = [](std::string s) {
    auto tokens = ekat::split(ekat::strip(s,' '),"=");
    EKAT_REQUIRE_MSG (tokens.size()==2,
        "[ModelInit::set_constant_fields] Error! Badly formatted entry.\n"
        " - entry: '" + s + "'\n"
        " - expected format: 'name = VALUE'\n");

    // Attempt to convert tokens[1] to a double
    double val;
    auto [ptr, ec] = std::from_chars(tokens[1].data(), tokens[1].data() + tokens[1].size(), val);

    if (ec != std::errc{} or ptr != tokens[1].data() + tokens[1].size()) {
      // If the conversion failed, we check if tokens[1] is the name of a known physics constant
      auto& pc_dict = physics::Constants<double>::dictionary();
      EKAT_REQUIRE_MSG (pc_dict.count(tokens[1]==0),
          "[ModelInit::set_constant_fields] Error! Badly formatted entry.\n"
          " Entry '" + tokens[1] + "' is neither a double nor a known physics constant.\n");
      val = pc_dict.at(tokens[1]).value;
    }
    return std::make_pair(tokens[0],val);
  };

  const auto& names_and_vals = m_params.get<strvec_t>("constant_fields");
  strmap_t<double> name2val;
  strmap_t<bool> name2used;
  for (const auto& item : name_and_vals) {
    auto [name, val] = key_val(item);
    name2val[name] = val;
  }

  for (auto it=fields.begin(); it!=fields.end()) {
    auto it2 = name2val.find(it->name());
    if (it2==name2val.end()) {
      ++it;
    } else {
      it->deep_copy(it2->second);
      it = fields.erase(it);
      name2used[it2->first] = true;
    }
  }

  EKAT_REQUIRE_MSG (name2val.size()==name2used.size(),
      "[ModelInit::set_constant_fields] Error! Unexpected constant field.\n"
      " - requested constant fields: " + ekat::join(name2val,[](auto it){return it.first},",") + "\n"
      " - atm input fields: " + ekat::join(fields,[](auto it){return it.name()},",") + "\n"); 
}

} // namespace scream
