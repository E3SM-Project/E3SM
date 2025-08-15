#ifndef SCREAM_FIELD_REQUEST_HPP
#define SCREAM_FIELD_REQUEST_HPP

#include "share/field/field_identifier.hpp"
#include "share/field/field_alloc_prop.hpp"
#include "share/util/eamxx_utils.hpp"

namespace scream {

enum RequestType {
  Required,
  Computed,
  Updated   // For convenience, triggers Required+Computed
};

// Whether a tracer should be advected by both Dynamics
// and Turbulance, or only by Dynamics
enum TracerAdvection {
  NoPreference, // Default. In the case that no process gives a preference,
                // the tracer is advected by dynamics and turbulence
  DynamicsAndTurbulence,
  DynamicsOnly,
};

// Whether the field group should be allocated as a monolithic field
enum class MonolithicAlloc : int {
  Required,
  NotRequired
};

/*
 * A struct used to request a group of fields.
 *
 * Groups are simply a labels attached to a Field object (see field_tracking.hpp).
 * They can be useful when a class need to access a certain group of fields,
 * in a way that is agnostic to how many fields are in said group.
 * A GroupRequest is a lightweight struct that an AP can expose if it needs a group
 * of fields, without caring how many there are, or how they are called.
 * A typical example is an AP that needs to advect tracers (like Dynamics does):
 * it treats tracers agnostically, and does not really care how many there are.
 * So the AP exposes this need as a GroupRequest. Later, it will be provided
 * with a FieldGroup, which allows to access all the fields in the group
 * individually, and, if the allocation permits it, as a single N+1 dimensional
 * field. For more details about the FieldGroup struct, see field_group.hpp.
 */
struct GroupRequest {
  // Main constructor method. Here is what the args are:
  //  - name: the name of the group
  //  - grid: the grid where the group is requested
  //  - ps: the pack size that the allocation of the fields in the group
  //        (and the monolithic field, if any) should accommodate (see field_alloc_prop.hpp)
  //  - monolithic_alloc: whether the group should be allocated as a monolithic group (see field_group.hpp)
  GroupRequest (const std::string& name_, const std::string& grid_, const int ps,
                const MonolithicAlloc monolithic_alloc_ = MonolithicAlloc::NotRequired)
   : name(name_), grid(grid_), pack_size(ps), monolithic_alloc(monolithic_alloc_)
  {
    EKAT_REQUIRE_MSG(pack_size>=1, "Error! Invalid pack size request.\n");
  }

  GroupRequest (const std::string& name_, const std::string& grid_,
                const MonolithicAlloc monolithic_alloc_ = MonolithicAlloc::NotRequired)
   : GroupRequest(name_,grid_,1,monolithic_alloc_)
  { /* Nothing to do here */ }

  // Default copy ctor is perfectly fine
  GroupRequest (const GroupRequest&) = default;

  // Main parts of a group request
  std::string name;                  // Group name
  std::string grid;                  // Grid name
  int pack_size;                     // Request an allocation that can accomodate Pack<Real,pack_size>
  MonolithicAlloc monolithic_alloc;  // Whether the group should be allocated as a single n+1 dimensional field
};

// In order to use GroupRequest in std sorted containers (like std::set),
// we need to provide an overload of op< (or a specialization of  std::less<T>).
inline bool operator< (const GroupRequest& lhs,
                       const GroupRequest& rhs)
{
  // Order by group name first
  if (lhs.name<rhs.name) {
    return true;
  } else if (lhs.name>rhs.name) {
    return false;
  }

  // Same group name, order by grid
  if (lhs.grid<rhs.grid) {
    return true;
  } else if (lhs.grid>rhs.grid) {
    return false;
  }

  // Same grid name, order by pack size
  if (lhs.pack_size < rhs.pack_size) {
    return true;
  } else if (lhs.pack_size>rhs.pack_size) {
    return false;
  }

  // Same pack size, order by monolithic allocation
  return etoi(lhs.monolithic_alloc)<etoi(rhs.monolithic_alloc);
}

/*
 * A struct used to request a field.
 *
 * The request contains at least a FieldIdentifier, but can also contain
 * a pack size, and a list of names of groups that the field should belong to.
 */

struct FieldRequest {
  using Units = ekat::units::Units;
  using FID   = FieldIdentifier;

  // Main constructor. Here is what the args are:
  //  - fid: the FieldIdentifier for the requested field (see field_identifier.hpp)
  //  - groups_: a list of groups that this field should be added to (see field_group.hpp)
  //  - ps: the pack size that the allocation of the field should accommodate (see field_alloc_prop.hpp)
  FieldRequest (const FID& fid_, const std::list<std::string>& groups_, const int ps)
   : fid(fid_), pack_size(ps), groups(groups_)
  {
    // Sanity checks
    EKAT_REQUIRE_MSG (ps>=1,
        "Error! Pack sizes must be >= 1.\n");
    // This seems funky, but write down a pow of 2 in binary format, and you'll see why it works
    EKAT_REQUIRE_MSG ( (ps & (ps-1))==0,
        "Error! We only support pack sizes that are (positive) powers of 2.\n");
  }

  // Convenience constructors.
  // They first three allow defaulting ps/groups_ and/or allow simpler sintax for single-group requests.
  // The last three do the same, but also allow to build a FieldIdentifier on the fly.
  FieldRequest (const FID& fid, const int ps = 1)
   : FieldRequest(fid,std::list<std::string>{},ps)
  { /* Nothing to do here */ }
  FieldRequest (const FID& fid, const std::list<std::string>& groups_)
   : FieldRequest(fid,groups_,1)
  { /* Nothing to do here */ }
  FieldRequest (const FID& fid, const std::string& group, const int ps = 1)
   : FieldRequest(fid,std::list<std::string>{group},ps)
  { /* Nothing to do here */ }
  FieldRequest (const std::string& name, const FieldLayout& layout, const Units& u, const std::string& grid,
                const std::list<std::string>& groups_, const int ps = 1)
   : FieldRequest(FID(name,layout,u,grid),groups_,ps)
  { /* Nothing to do here */ }
  FieldRequest (const std::string& name, const FieldLayout& layout, const Units& u, const std::string& grid,
                const int ps = 1)
   : FieldRequest(FID(name,layout,u,grid),std::list<std::string>{},ps)
  { /* Nothing to do here */ }
  FieldRequest (const std::string& name, const FieldLayout& layout, const Units& u, const std::string& grid,
                const std::string& group, const int ps = 1)
   : FieldRequest(FID(name,layout,u,grid),std::list<std::string>{group},ps)
  { /* Nothing to do here */ }

  FieldRequest (const FID& fid, const FieldRequest& parent, int idim, int k, bool dynamic)
   : FieldRequest (fid)
  {
    subview_info.dim_idx = idim;
    subview_info.slice_idx = k;
    subview_info.dim_extent = parent.fid.get_layout().dim(idim);
    subview_info.dynamic = dynamic;

    parent_name = parent.fid.name();
  }

  FieldRequest (const std::string& field_name, const std::string& grid_name,
                const std::list<std::string>& groups = {}, const int ps = 1)
   : FieldRequest (incomplete_fid(field_name,grid_name),groups,ps)
  {
    incomplete = true;
  }

  static FieldIdentifier
  incomplete_fid (const std::string& field_name, const std::string& grid_name)
  {
    return FieldIdentifier(field_name,FieldLayout::invalid(),Units::invalid(),grid_name,DataType::Invalid);
  }

  // Data
  FieldIdentifier           fid;
  int                       pack_size;
  std::list<std::string>    groups;
  SubviewInfo               subview_info;
  std::string               parent_name;
  bool                      incomplete = false;
  std::string               calling_process = "UNKNOWN";
};

// In order to use FieldRequest in std sorted containers (like std::set),
// we need to provide an overload of op< (or a specialization of  std::less<T>).
inline bool operator< (const FieldRequest& lhs,
                       const FieldRequest& rhs)
{
  if (lhs.fid<rhs.fid) {
    return true;
  } else if (lhs.fid==rhs.fid) {
    if (lhs.pack_size<rhs.pack_size) {
      return true;
    } else if (lhs.pack_size==rhs.pack_size) {
      if (lhs.groups < rhs.groups) {
        return true;
      } else if (lhs.groups==rhs.groups) {
        return lhs.calling_process < rhs.calling_process;
      }
    }
  }
  return false;
}

} // namespace scream

#endif // SCREAM_FIELD_REQUEST_HPP
