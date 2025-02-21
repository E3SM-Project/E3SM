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

// Whether the bundling of a field group (see below) is needed, optional, or not needed.
enum class Bundling : int {
  Required,
  Preferred,
  NotNeeded
};

// What's the relation of group A and group B. In particular, if group A is 'derived'
// from group B, this enum explains how it is derived.
//  - None: not a "derived" request
//  - Import: group B is on a different grid than group A. This derivation type makes
//            sure that there is a replica of B on the grid associated with A.
//  - Copy: a hard copy of the group B is created, with only the "bundled" field being
//          allocated, to avoid creating two copies of the samae field.
//          This derivation type requires group A and group B to be on the *same* grid.
//  - Superset: all fields in A also appear in B (the *same* fields)
//  - Subset: all fields in B also appear in A (the *same* fields)
enum class DerivationType : int {
  None,
  Import,
  Copy,
  Superset,
  Subset
};

inline std::string e2str (const DerivationType rt) {
  switch (rt) {
    case DerivationType::None:      return "None";
    case DerivationType::Import:    return "Import";
    case DerivationType::Copy:      return "Copy";
    case DerivationType::Subset:    return "Subset";
    case DerivationType::Superset:  return "Superset";
  }
  return "INVALID";
}

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
  //        (and the bundled field, if any) should accommodate (see field_alloc_prop.hpp)
  //  - bundling: whether the group should be bundled (see field_group.hpp)
  //  - r: allows to specify specs of this request in terms of another.
  //  - t: the relationshipt of group in request r compared to this one.
  //       E.g.: t=Superset means r->name is a superset of group this->name.
  //  - exclude: if t=Superset/Subset, this group contains all field in r, *except/plus* those in this list.
  GroupRequest (const std::string& name_, const std::string& grid_, const int ps, const Bundling b,
                const DerivationType rt, const std::string& src_name_, const std::string& src_grid_,
                const std::list<std::string>& exclude_ = {})
   : name(name_), grid(grid_), pack_size(ps), bundling(b), derived_type(rt)
  {
    EKAT_REQUIRE_MSG(pack_size>=1, "Error! Invalid pack size request.\n");
    if (rt!=DerivationType::None) {
      src_name = src_name_;
      src_grid = src_grid_;

      EKAT_REQUIRE_MSG (src_grid==grid || rt==DerivationType::Import,
          "Error! We only allow cross-grid group derivation if the derivation type is 'Import'.");

      EKAT_REQUIRE_MSG (exclude_.size()==0 || rt==DerivationType::Subset,
          "Error! You can only exclude fields when deriving a group from another if this group\n"
          "       is a subset of the source group.\n");
      exclude = exclude_;
    }
  }

  // Convenience ctors when some features are not needed
  GroupRequest (const std::string& name_, const std::string& grid_,
                const int ps, const Bundling b = Bundling::NotNeeded)
   : GroupRequest(name_,grid_,ps,b,DerivationType::None,"","",{})
  { /* Nothing to do here */ }

  GroupRequest (const std::string& name_, const std::string& grid_,
                const Bundling b = Bundling::NotNeeded)
   : GroupRequest(name_,grid_,1,b)
  { /* Nothing to do here */ }

  // Default copy ctor is perfectly fine
  GroupRequest (const GroupRequest&) = default;
  
  // Main parts of a group request
  std::string name;   // Group name
  std::string grid;   // Grid name
  int pack_size;      // Request an allocation that can accomodate Pack<Real,pack_size>
  Bundling bundling;  // Whether the group should be allocated as a single n+1 dimensional field

  // The following members allow to specify a request in terms of another group.
  // A possible use of this is when an atm proc wants to create G1 "excluding"
  // some fields from G2, and have the remaining ones still contiguous in mem (i.e.,
  // accessible with a bundled array). This will inform the FM to rearrange the fields
  // in G2 so that the subset of fields that are in G1 appear contiguously.
  // Another use is to create group G2 on grid B to contain the same fields of
  // group G1 on grid A, without knowing what's in G1 a priori.
  // See comments in FieldManager::registration_ends() for more details
  // Note: derived_type is what *this group* is for $src_name. So, if derived_type=Subset,
  //       then this group is a subset of the group $src_name.
  DerivationType derived_type;
  std::string src_name;
  std::string src_grid;
  std::list<std::string> exclude;  // Only for derived_type=Subset
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

  // Same pack size, order by bundling
  if (etoi(lhs.bundling)<etoi(rhs.bundling)) {
    return true;
  } else if (etoi(lhs.bundling)>etoi(rhs.bundling)) {
    return false;
  }

  // Same bundling, order by derivation type
  if (etoi(lhs.derived_type)<etoi(rhs.derived_type)) {
    return true;
  } else if (etoi(lhs.derived_type)<etoi(rhs.derived_type)) {
    return false;
  }

  // Same derivation type, order by source group name
  if (lhs.src_name<rhs.src_name) {
    return true;
  } else if (lhs.src_name>rhs.src_name) {
    return false;
  }

  // Same souce group name, order by source group grid
  if (lhs.src_grid<rhs.src_grid) {
    return true;
  } else if (lhs.src_grid>rhs.src_grid) {
    return false;
  }

  // Same source group grid, order by exclude fields
  return (lhs.exclude<rhs.exclude);
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
      return lhs.groups < rhs.groups;
    }
  }
  return false;
}

} // namespace scream

#endif // SCREAM_FIELD_REQUEST_HPP
