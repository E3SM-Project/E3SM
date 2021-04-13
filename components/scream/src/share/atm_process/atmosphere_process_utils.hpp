#ifndef SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP

#include "share/field//field_identifier.hpp"
#include "ekat/ekat_assert.hpp"

#include <string>
#include <list>
#include <memory>
#include <type_traits>

namespace scream {

enum class AtmosphereProcessType {
  Coupling,   // Process responsible of interfacing with the component coupler
  Dynamics,   // Process responsible of handling the dynamics
  Physics,    // Process handling a physics parametrization
  Group       // Process that groups a bunch of processes (so they look as a single process)
};

inline std::string e2str (const AtmosphereProcessType ap_type) {
  switch (ap_type) {
    case AtmosphereProcessType::Coupling:  return "Surface Coupling";
    case AtmosphereProcessType::Dynamics:  return "Atmosphere Dynamics";
    case AtmosphereProcessType::Physics:   return "Atmosphere Physics Parametrization";
    case AtmosphereProcessType::Group:     return "Atmosphere Process Group";
    default:
      ekat::error::runtime_abort("Error! Unrecognized atmosphere process type.\n");
  }
  return "INVALID";
}

// This enum is mostly used by AtmosphereProcessGroup to establish whether
// its atm procs are to be run concurrently or sequentially.
// We put the enum here so other files can easily access it.
enum class ScheduleType {
  Sequential,
  Parallel
};

// Whether a the bundling of a field group (see below) is needed, optional, or not needed.
enum class Bundling : int {
  Unspecified = 0,
  Required  = 1,
  Preferred = 2,
  NotNeeded = 3
};

enum class ChildGroupRequestType : int {
  Unspecified   = 0,
  SoftCopyAlias = 1,
  HardCopyAlias = 2,
  Subset        = 3
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
  // Main constructor method
  GroupRequest (const std::string& name_, const std::string& grid_, const int ps, const Bundling b,
                const GroupRequest* p, const ChildGroupRequestType t, const std::list<std::string>& excl)
   : name(name_), grid(grid_), pack_size(ps), bundling(b)
  {
    EKAT_REQUIRE_MSG(pack_size>=1, "Error! Invalid pack size request.\n");
    if (p!=nullptr) {
      parent = std::make_shared<GroupRequest>(*p);
      child_type = t;
      EKAT_REQUIRE_MSG (excl.size()==0 || t==ChildGroupRequestType::Subset,
          "Error! You can only exclude fields from a parent group if creating a Subset child group.\n");
      exclude = excl;

      // TODO: should we relax this? Not allowing multiple levels of nested
      //       parenting makes it easier for the AD and FM to correctly allocated
      //       fields...
      EKAT_REQUIRE_MSG (p->parent==nullptr,
          "Error! We cannot handle multiple levels of nested groups.\n");
    }
  }

  // Convenience overloads of the ctor
  GroupRequest (const std::string& name_, const std::string& grid_,
                const int ps, const Bundling b = Bundling::NotNeeded)
   : GroupRequest(name_,grid_,ps,b,nullptr,ChildGroupRequestType::Unspecified,{})
  { /* Nothing to do here */ }
  GroupRequest (const std::string& name_, const std::string& grid_,
                const Bundling b = Bundling::NotNeeded)
   : GroupRequest(name_,grid_,1,b)
  { /* Nothing to do here */ }

  // Default copy ctor is perfectly fine
  GroupRequest (const GroupRequest&) = default;

  // Group name
  std::string name;
  // Grid name
  std::string grid;
  // Request an allocation that can accomodate a value type like Pack<Real,pack_size>
  int       pack_size;

  // The following members allow to specify a request in terms of another group.
  // E.g., one can ask for group G1 to be an alias of G2, which will tell the FM
  // to create a field group G1 "equivalent" to G2. The alias can be soft or hard.
  // If soft, the group bundle of G1 will store the same view of the bundle in G2
  // (assuming G2 is bundled), and same for the fields. A hard alias, will create
  // G1 as a separate bundled field, and will *not* create the fields corresponding
  // to the individual group members (the user can still "subview" the bundled group
  // at particular entries, of course).
  // Another use of the parent group is when an atm proc wants to create G1 "excluding"
  // some fields from G2, and have the remaining ones still contiguous in mem (i.e.,
  // accessible with a bundled array). This will inform the FM to rearrange the fields
  // in G2 so that the subset of fields that are in G1 appear contiguously. Clearly,
  // the FM can only accommodate certain requests of this type, so if not possible,
  // the FM will throw.
  std::shared_ptr<GroupRequest>   parent;

  Bundling bundling;
  ChildGroupRequestType child_type;
  std::list<std::string> exclude;
};

template<typename EnumT>
constexpr typename
std::enable_if<std::is_enum<EnumT>::value,
                 typename std::underlying_type<EnumT>::type
              >::type
etoi (const EnumT e) {
  return static_cast<typename std::underlying_type<EnumT>::type>(e);
}

// In order to use GroupRequest in std sorted containers (like std::set),
// we need to provide an overload of op< or std::less.
inline bool operator< (const GroupRequest& lhs,
                       const GroupRequest& rhs)
{
  if (lhs.name<rhs.name) {
    return true;
  } else if (lhs.name==rhs.name) {
    if (lhs.grid<rhs.grid) {
      return true;
    } else if (lhs.grid==rhs.grid) {
      if (lhs.pack_size < rhs.pack_size) {
        return true;
      } else if (lhs.pack_size==rhs.pack_size) {
        if (lhs.parent==nullptr) {
          return rhs.parent!=nullptr;
        } else if (rhs.parent!=nullptr) {
          if ( (*lhs.parent) < (*rhs.parent) ) {
            return true;
          } else if ( ! ( (*rhs.parent)<(*lhs.parent)) ) {
            return etoi(lhs.bundling) < etoi(rhs.bundling) ||
                   (etoi(lhs.bundling) == etoi(rhs.bundling) &&
                    etoi(lhs.child_type) < etoi(rhs.child_type));
          }
        }
      }
    }
  }
  return false;
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

  // Main constructor
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

  // Convenience constructors
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

  // Data
  FieldIdentifier           fid;
  int                       pack_size;
  std::list<std::string>    groups;
};

// In order to use FieldRequest in std sorted containers (like std::set),
// we need to provide an overload of op< or std::less.
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

#endif // SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP
