#ifndef SCREAM_FIELD_REQUEST_HPP
#define SCREAM_FIELD_REQUEST_HPP

#include "share/field/field_identifier.hpp"
#include "share/field/field_alloc_prop.hpp"
#include "share/util/eamxx_utils.hpp"

namespace scream {

enum RequestType {
  Invalid  = 0,  // To spot if a request was not set correctly
  Required = 1,
  Computed = 2,
  Updated  = 3,  // For convenience, triggers Required+Computed
};

// Allow bitwise op on request types
inline RequestType& operator&= (RequestType& lhs, RequestType rhs)
{
  lhs = static_cast<RequestType>(static_cast<int>(lhs) & static_cast<int>(rhs));
  return lhs;
}
inline RequestType& operator|= (RequestType& lhs, RequestType rhs)
{
  lhs = static_cast<RequestType>(static_cast<int>(lhs) | static_cast<int>(rhs));
  return lhs;
}
inline RequestType operator& (RequestType lhs, RequestType rhs)
{
  return static_cast<RequestType>(static_cast<int>(lhs) & static_cast<int>(rhs));
}
inline RequestType operator| (RequestType lhs, RequestType rhs)
{
  return static_cast<RequestType>(static_cast<int>(lhs) | static_cast<int>(rhs));
}

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

// Base class for the Field and Group requests
// NOTE: the CRTP is used to allow set_ps and add_usage to return refs to the derived type.
//       Otherwise, we'd need to add an impl in both derived classes
template<typename Derived>
struct Request{
  Request (const std::string& name_, const std::string& grid_)
   : name(name_), grid(grid_) {}

  Derived& set_ps (int ps) {
    EKAT_REQUIRE_MSG (ps>=1,
        "Error! Pack sizes must be >= 1.\n");
    // This seems funky, but write down a pow of 2 in binary format, and you'll see why it works
    EKAT_REQUIRE_MSG ( (ps & (ps-1))==0,
        "Error! We only support pack sizes that are (positive) powers of 2.\n");
    pack_size = ps;
    return static_cast<Derived&>(*this);
  }
  Derived& add_usage (RequestType usage_) {
    usage |= usage_;
    return static_cast<Derived&>(*this);
  }
  Derived& set_caller (const std::string& caller_) {
    caller = caller_;
    return static_cast<Derived&>(*this);
  }


  std::string  name;                // the name of the resource requested
  std::string  grid;                // the grid where the resource is needed
  RequestType  usage = Invalid;     // how the resource will be used
  int          pack_size = 1;       // the needed pack (SIMD) length for the resource
  std::string  caller = "UNKNOWN";  // for debug purposes
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
struct GroupRequest : Request<GroupRequest> {
  // Main constructor method. Here is what the args are:
  //  - name: the name of the group
  //  - grid: the grid where the group is requested
  //  - monolithic_alloc: whether the group should be allocated as a monolithic group (see field_group.hpp)
  GroupRequest (const std::string& name_, const std::string& grid_)
   : Request(name_, grid_) {}

  // Default copy ctor is perfectly fine
  GroupRequest (const GroupRequest&) = default;

  GroupRequest& set_alloc (MonolithicAlloc alloc) {
    monolithic_alloc = alloc;
    return *this;
  }

  // Main parts of a group request
  std::string name;                  // Group name
  MonolithicAlloc monolithic_alloc;  // Whether the group should be allocated as a single n+1 dimensional field
};

// // In order to use GroupRequest in std sorted containers (like std::set),
// // we need to provide an overload of op< (or a specialization of  std::less<T>).
// inline bool operator< (const GroupRequest& lhs,
//                        const GroupRequest& rhs)
// {
//   // Order by group name first
//   if (lhs.name<rhs.name) {
//     return true;
//   } else if (lhs.name>rhs.name) {
//     return false;
//   }

//   // Same group name, order by grid
//   if (lhs.grid<rhs.grid) {
//     return true;
//   } else if (lhs.grid>rhs.grid) {
//     return false;
//   }

//   // Same grid name, order by pack size
//   if (lhs.pack_size < rhs.pack_size) {
//     return true;
//   } else if (lhs.pack_size>rhs.pack_size) {
//     return false;
//   }

//   // Same pack size, order by monolithic allocation
//   return etoi(lhs.monolithic_alloc)<etoi(rhs.monolithic_alloc);
// }

/*
 * A struct used to request a field.
 *
 * The request contains at least a FieldIdentifier, but can also contain
 * a pack size, and a list of names of groups that the field should belong to.
 */

struct FieldRequest : Request<FieldRequest> {
  using FID   = FieldIdentifier;

  // Main constructor. Here is what the args are:
  //  - fid: the FieldIdentifier for the requested field (see field_identifier.hpp)
  //  - groups_: a list of groups that this field should be added to (see field_group.hpp)
  //  - ps: the pack size that the allocation of the field should accommodate (see field_alloc_prop.hpp)
  FieldRequest (const FID& fid_, const std::string& grid_)
   : Request(fid_.name(), grid_), fid(fid_)
  { /* Nothing to do here */ }

  FieldRequest& set_groups(const std::list<std::string>& groups_) {
    groups = groups_;
    return *this;
  }
  FieldRequest& set_groups(const std::string& group) {
    groups.push_back(group);
    return *this;
  }

  FieldRequest (const FID& fid, const FieldRequest& parent, int idim, int k, bool dynamic)
   : FieldRequest (fid,parent.grid)
  {
    subview_info.dim_idx = idim;
    subview_info.slice_idx = k;
    subview_info.dim_extent = parent.fid.get_layout().dim(idim);
    subview_info.dynamic = dynamic;

    parent_name = parent.fid.name();
  }

  FieldRequest (const std::string& field_name, const std::string& grid)
   : FieldRequest (incomplete_fid(field_name),grid)
  {
    incomplete = true;
  }

  static FieldIdentifier
  incomplete_fid (const std::string& field_name)
  {
    return FieldIdentifier(field_name,FieldLayout::invalid(),DataType::Invalid);
  }

  // Data
  FieldIdentifier           fid;
  std::list<std::string>    groups;
  SubviewInfo               subview_info;
  std::string               parent_name;
  bool                      incomplete = false;
};

// // In order to use FieldRequest in std sorted containers (like std::set),
// // we need to provide an overload of op< (or a specialization of  std::less<T>).
// inline bool operator< (const FieldRequest& lhs,
//                        const FieldRequest& rhs)
// {
//   if (lhs.fid<rhs.fid) {
//     return true;
//   } else if (lhs.fid==rhs.fid) {
//     if (lhs.pack_size<rhs.pack_size) {
//       return true;
//     } else if (lhs.pack_size==rhs.pack_size) {
//       if (lhs.groups < rhs.groups) {
//         return true;
//       } else if (lhs.groups==rhs.groups) {
//         return lhs.calling_process < rhs.calling_process;
//       }
//     }
//   }
//   return false;
// }

} // namespace scream

#endif // SCREAM_FIELD_REQUEST_HPP
