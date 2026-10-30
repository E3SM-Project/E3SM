#include "share/data_managers/model_init.hpp"

#include "share/field/field_reader.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_assert.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <ranges>

namespace scream
{

ModelInit::
ModelInit (const ekat::ParameterList& params)
 : m_params (params)
{
  // The constants are in the 'constant_fields' sublist, which has two sublists
  // (if not present, they are created empty): 'overrides' and 'defaults'.
  m_params.sublist("constant_fields").sublist("overrides");
  m_params.sublist("constant_fields").sublist("defaults");
}

void ModelInit::
run (const std::shared_ptr<FieldManager>& fm,
     const util::TimeStamp& t0,
     const RunType run_type)
{
  const bool restart = run_type==RunType::Restart;

  if (not restart) {
    check_no_constants_for_composite_fields(fm);
  }

  const auto group_name = restart ? "RESTART" : "STARTUP";
  for (const auto& gn : fm->get_grids_manager()->get_grid_names()) {
    if (not fm->has_group(group_name,gn)) {
      continue;
    }
    if (restart) {
      init_restart_fields(fm,gn,t0);
    } else {
      init_startup_fields(fm,gn,t0);
    }
  }

  if (restart) {
    load_restart_extra_data();
  }
}

void ModelInit::
init_startup_fields (const std::shared_ptr<FieldManager>& fm,
                     const std::string& grid_name,
                     const util::TimeStamp& t0)
{
  const auto& constants_pl = m_params.sublist("constant_fields");
  const auto& overrides_pl = constants_pl.sublist("overrides");
  const auto& defaults_pl  = constants_pl.sublist("defaults");

  // Check if a field is in the IC file (if any). The file is opened on the fly (and released
  // right away), and only if we need to check for a field with a default constant.
  const bool has_ic_file = m_params.isParameter("filename");
  const std::string ic_file = has_ic_file ? m_params.get<std::string>("filename") : "";
  auto in_ic_file = [&](const std::string& fname) {
    return has_ic_file and scorpio::has_var(ic_file,fname);
  };

  // Process all fields in the STARTUP group. For each, either init to
  // a constant (if provided), or add it to the list of fields to read from file.
  std::vector<Field> inited;
  std::vector<Field> to_read;
  for (auto f : get_fields(fm,"STARTUP",grid_name)) {
    const auto& fname = f.name();
    if (overrides_pl.isParameter(fname)) {
      set_constant_field(f,overrides_pl);
      inited.push_back(f);
    } else if (defaults_pl.isParameter(fname) and not in_ic_file(fname)) {
      set_constant_field(f,defaults_pl);
      inited.push_back(f);
    } else {
      to_read.push_back(f);
    }
  }

  if (to_read.size()>0) {
    EKAT_REQUIRE_MSG (has_ic_file,
        "Error! Some fields need to be loaded from the IC file, but no 'filename' "
        "was found in the initial conditions parameters.\n"
        " - grid name: " + grid_name + "\n"
        " - fields: " + to_read.front().name() + (to_read.size()>1 ? ", ..." : "") + "\n");
    read_fields(ic_file,to_read,fm->get_grids_manager()->get_grid(grid_name));
    inited.insert(inited.end(),to_read.begin(),to_read.end());
  }

  update_time_stamps(inited,t0);
}

void ModelInit::
init_restart_fields (const std::shared_ptr<FieldManager>& fm,
                     const std::string& grid_name,
                     const util::TimeStamp& t0)
{
  auto fields = get_fields(fm,"RESTART",grid_name);
  if (fields.size()==0) {
    return;
  }

  EKAT_REQUIRE_MSG (m_params.isParameter("filename"),
      "Error! Some fields need to be loaded from the restart file, but no 'filename' "
      "was found in the initial conditions parameters.\n"
      " - grid name: " + grid_name + "\n");
  const auto& filename = m_params.get<std::string>("filename");

  try {
    read_fields(filename,fields,fm->get_grids_manager()->get_grid(grid_name));
  } catch (std::exception& e) {
    EKAT_ERROR_MSG (
        std::string(e.what()) + "\n"
        "NOTE: if the error above complains about a missing variable 'U' or 'V',\n"
        "      you may be trying to restart from a restart file created before\n"
        "      PR #8606 (https://github.com/E3SM-Project/E3SM/pull/8606), which\n"
        "      changed EAMxx to init/restart 'U'/'V' directly, rather than the\n"
        "      vector field 'horiz_winds'. See the 'Restart file compatibility'\n"
        "      section in components/eamxx/docs/user/model_configuration.md for\n"
        "      instructions on how to convert an old restart file.\n");
  }

  update_time_stamps(fields,t0);
}

void ModelInit::
load_restart_extra_data () const
{
  if (m_restart_extra_data.size()==0) {
    return;
  }

  EKAT_REQUIRE_MSG (m_params.isParameter("filename"),
      "Error! Restart extra data needs to be loaded from the restart file, but no 'filename' "
      "was found in the initial conditions parameters.\n");
  const auto& filename = m_params.get<std::string>("filename");

  for (const auto& it : m_restart_extra_data) {
    const auto& name = it.first;
          auto& any  = *it.second;

    if (any.type()==typeid(int)) {
      std::any_cast<int&>(any) = scorpio::get_attribute<int>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(std::int64_t)) {
      std::any_cast<std::int64_t&>(any) = scorpio::get_attribute<std::int64_t>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(float)) {
      std::any_cast<float&>(any) = scorpio::get_attribute<float>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(double)) {
      std::any_cast<double&>(any) = scorpio::get_attribute<double>(filename,"GLOBAL",name);
    } else if (any.type()==typeid(std::string)) {
      std::any_cast<std::string&>(any) = scorpio::get_attribute<std::string>(filename,"GLOBAL",name);
    } else {
      EKAT_ERROR_MSG (
          "Error! Unrecognized/unsupported concrete type for restart extra data.\n"
          " - extra data name  : " + name + "\n"
          " - extra data typeid: " + std::string(any.type().name()) + "\n");
    }
  }
}

std::vector<Field> ModelInit::
get_fields (const std::shared_ptr<FieldManager>& fm,
            const std::string& group_name,
            const std::string& grid_name) const
{
  std::vector<Field> fields;
  auto group = fm->get_field_group(group_name,grid_name);
  for (const auto& f : std::views::values(group.individual_fields())) {
    // Skip fields that were already inited
    if (not f.get_header().get_tracking().get_time_stamp().is_valid()) {
      fields.push_back(f);
    }
  }
  return fields;
}

void ModelInit::
read_fields (const std::string& filename,
             const std::vector<Field>& fields,
             const std::shared_ptr<const AbstractGrid>& grid) const
{
  scream::read_fields(filename,fields,grid->get_partitioned_dim_gids(),grid->get_comm());
}

void ModelInit::
update_time_stamps (const std::vector<Field>& fields,
                    const util::TimeStamp& t0) const
{
  for (auto f : fields) {
    f.get_header().get_tracking().update_time_stamp(t0);

    // Also init the timestamp of children (if any)
    for (auto& c : f.get_header().get_children()) {
      c.lock()->get_tracking().update_time_stamp(t0);
    }

    // If there is a parent, and all children were inited,
    // also init the parent's timestamp
    if (auto p = f.get_header().get_parent(); p!=nullptr) {
      bool p_inited = true;
      for (auto c : p->get_children()) {
        p_inited &= c.lock()->get_tracking().get_time_stamp().is_valid();
      }
      if (p_inited) {
        p->get_tracking().update_time_stamp(t0);
      }
    }
  }
}

void ModelInit::
check_no_constants_for_composite_fields (const std::shared_ptr<FieldManager>& fm)
{
  const auto& constants_pl = m_params.sublist("constant_fields");
  const auto& overrides_pl = constants_pl.sublist("overrides");
  const auto& defaults_pl  = constants_pl.sublist("defaults");

  for (const auto& gn : fm->get_grids_manager()->get_grid_names()) {
    for (const auto& it : fm->get_repo(gn)) {
      const auto& f = *it.second;
      const auto& children = f.get_header().get_children();
      if (children.size()>0 and
          (overrides_pl.isParameter(f.name()) or defaults_pl.isParameter(f.name()))) {
        std::string child_names;
        for (auto c : children) {
          child_names += c.lock()->get_identifier().name() + " ";
        }
        EKAT_ERROR_MSG (
            "Error! Cannot set an initial condition for field '" + f.name() + "' directly, "
            "since it has subfields.\n"
            "       Grid name:     " + gn + "\n"
            "       Field:         " + f.name() + "\n"
            "       Subfields:     " + child_names + "\n"
            "       Please, set the initial condition for each subfield individually.\n");
      }
    }
  }
}

void ModelInit::
set_constant_field (Field& f, const ekat::ParameterList& constants_pl) const
{
  const auto& fid = f.get_header().get_identifier();
  const auto& layout = fid.get_layout();
  const auto& name = fid.name();

  // The constant can be specified in one of these forms
  //   fname: val
  //   fname: [val]
  //   fname: [val1,...,valN]
  // where the values are reals (write 1.0, not 1, since the ekat yaml parser reads 1 as an int).
  // The XML defaults enforce that all constants are array(real), so that is what atmchange
  // and buildnml produce. The scalar form is only accepted to keep the input.yaml files of
  // standalone tests simple(r).
  // In the first two cases, all entries of the field are inited to val (for any field).
  // In the last case, the field must be a vector field with N components, and each
  // component is inited to the corresponding entry.
  std::vector<double> values;
  if (constants_pl.isType<std::vector<double>>(name)) {
    values = constants_pl.get<std::vector<double>>(name);
  } else if (constants_pl.isType<double>(name)) {
    values.push_back(constants_pl.get<double>(name));
  } else {
    EKAT_ERROR_MSG ("Error! Invalid initial condition for field '" + name + "': only scalar "
                    "or vector double arguments are allowed.\n"
                    "Note: If you used int values, convert to double (e.g., write 1.0 instead of 1).\n");
  }

  // Check all inputs before modifying the field
  const int vec_dim = layout.is_vector_layout() ? layout.get_vector_dim() : 1;
  EKAT_REQUIRE_MSG (values.size()==1 or values.size()==static_cast<size_t>(vec_dim),
      "Error! Initial condition values array for '" + name + "' has the wrong dimension.\n"
      "       Grid name:       " + fid.get_grid_name() + "\n"
      "       Field layout:    " + layout.to_string() + "\n"
      "       Array dimension: " + std::to_string(values.size()) + "\n"
      "       Valid dimensions: 1 (set all entries to the same value)" +
      (vec_dim>1 ? " or " + std::to_string(vec_dim) + " (one value per vector component)\n" : "\n"));

  // Integer fields can be set from (real) values, but the conversion must be exact
  // (e.g., 2.0 is ok, 2.5 throws), since we do not want to silently truncate. The value
  // must also fit in an int (this also rules out inf, and NaN fails the first check).
  const bool is_int = f.data_type()==DataType::IntType;
  if (is_int) {
    for (const auto v : values) {
      EKAT_REQUIRE_MSG (v==std::floor(v) and std::abs(v)<=std::numeric_limits<int>::max(),
          "Error! Cannot initialize the integer field '" + name + "' with a non-integral value, "
          "or a value outside the range of int.\n"
          "       Value: " + std::to_string(v) + "\n");
    }
  }

  auto set_value = [&](const Field& x, const double v) {
    if (is_int) {
      x.deep_copy(static_cast<int>(v));
    } else {
      x.deep_copy(v);
    }
  };

  if (values.size()==1) {
    set_value(f,values[0]);
  } else {
    // Extract a subfield for each component. This is not "too" expensive, especially
    // considering that this code is executed during initialization only.
    for (int comp=0; comp<vec_dim; ++comp) {
      set_value(f.get_component(comp),values[comp]);
    }
  }
}

} // namespace scream
