#include "share/data_managers/model_init.hpp"

#include "share/field/field_reader.hpp"
#include "share/physics/physics_constants.hpp"

#include <ekat_string_utils.hpp>

#include <ranges>
#include <charconv>

namespace scream
{

ModelInit::
ModelInit (const ekat::ParameterList& params)
 : m_params (params)
{
  // Build the map name->value for constant fields
  const auto& names_and_vals = m_params.get<strvec_t>("constant_fields");
  for (auto s : names_and_vals) {
    ekat::strip(s,' ');
    auto tokens = ekat::split(s,"=");
    EKAT_REQUIRE_MSG (tokens.size()==2,
        "[ModelInit::set_constant_fields] Error! Badly formatted entry.\n"
        " - entry: '" + s + "'\n"
        " - expected format: 'name = VALUE'\n");

    // Attempt to convert tokens[1] to a double
    double dval;
    int ival;
    auto beg = tokens[1].data();
    auto end = beg + tokens[1].size();

    if (auto [ptr, ec] = std::from_chars(beg,end,ival); ec == std::errc{} and ptr == end) {
      // The int conversion worked. Store the value
      m_constant_values[tokens[0]] = ival;
    } else if (auto [ptr, ec] = std::from_chars(beg,end,dval); ec == std::errc{} and ptr == end) {
      // The double conversion worked. Store the value
      m_constant_values[tokens[0]] = dval;
    } else {
      // The string is not a number. Check if it is the name of a known physics constant
      auto& pc_dict = physics::Constants<double>::dictionary();
      EKAT_REQUIRE_MSG (pc_dict.count(tokens[1])==0,
          "[ModelInit::set_constant_fields] Error! Badly formatted entry.\n"
          " Entry '" + tokens[1] + "' is neither a double nor a known physics constant.\n");
      m_constant_values[tokens[0]] = pc_dict.at(tokens[1]).value;
    }
  }
}

void ModelInit::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  auto model_group_name = run_type==RunType::Restart ? "RESTART" : "STARTUP";
  auto gm = fm->get_grids_manager();
  for (auto gn : gm->get_grid_names()) {
    auto grid = gm->get_grid(gn);
    if (fm->has_group("TOPOGRAPHY",gn)) {
      const auto& filename = m_params.get<std::string>("topography_filename");
      auto topo_fields = get_fields(fm,"TOPOGRAPHY",gn);
      set_constant_fields(topo_fields,t0);
      read_fields (filename,topo_fields,grid);
    }

    if (fm->has_group(model_group_name,gn)) {
      const auto& filename = m_params.get<std::string>("filename");
      auto model_fields = get_fields(fm,model_group_name,gn);
      set_constant_fields(model_fields,t0);
      read_fields (filename,model_fields,grid);
    }
  }
}

void ModelInit::
read_fields (const std::string& filename,
             std::vector<Field>& fields,
             const std::shared_ptr<const AbstractGrid>& grid)
{
  // Read the fields via common utils
  FieldReader reader;
  reader.set_fields(fields);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),grid->get_comm());
  reader.set_file_specs(filename,get_tag_rename(grid->name()));
  reader.read();
}

void ModelInit::
set_constant_fields (std::vector<Field>& fields,
                     const util::TimeStamp& t0)
{
  for (auto it_f = fields.begin(); it_f!=fields.end(); ) {
    const auto& name = it_f->name();
    auto it_v = m_constant_values.find(name);
    if (it_v!=m_constant_values.end()) {
      it_f->deep_copy(it_v->second);
      it_f->get_header().get_tracking().update_time_stamp(t0);
      it_f = fields.erase(it_f);
    } else {
      ++it_f;
    }
  }
}

std::vector<Field>
ModelInit::
get_fields (const std::shared_ptr<FieldManager>& fm,
            const std::string& group_name,
            const std::string& grid_name)
{
  auto group = fm->get_field_group(group_name,grid_name);
  std::vector<Field> fields;
  for (const auto& f : std::views::values(group.individual_fields())) {
    // If for some reason the field was already inited, skip it)
    if (f.get_header().get_tracking().get_time_stamp().is_valid())
      continue;
    fields.push_back(f);
  }
  return fields;
}

} // namespace scream
