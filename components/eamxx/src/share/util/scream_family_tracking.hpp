#ifndef SCREAM_FAMILY_TRACKING_CLASS
#define SCREAM_FAMILY_TRACKING_CLASS

#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/ekat_assert.hpp"

#include <list>
#include <memory>
#include <type_traits>

namespace scream
{

/*
 * Class that handles tracking of parent/children instances of the same class type
 *
 * This class is useful when we have a class T that needs to keep a pointer
 * to other instances of T, stored as weak pointer.  In particular, we allow
 * storing one "parent" and multiple "children". In order for the family
 * tracking mechanisms to kick in, the objects must be created in a way
 * that triggers the capabilities of the enable_shared_from_this class;
 * that is, calling 'weak_from_this' should always return a nonnull pointer
 * (though it can be expired, e.g., within this class' destructor).
 *
 * After creating a new instance of this class, you can create a parent-child
 * link by setting the parent of this class via 'create_parent_child_link'.
 * This method sets this object as a child in the input parent ptr, and sets
 * the input ptr as parent of this object. You *cannot* call this function
 * while this object is being constructed, since the shared pointer control
 * block is not yet built (hence, no weak_from_this available yet).
 *
 * When this object is destructed, the parent-child link will be removed.
 * Notice that if ~FamilyTracking is called, then the last shared ptr to
 * this object has been deleted/reset, so the ptr returned by weak_from_this
 * will be expired. Nevertheless, it will still store the same ptr, and
 * can still be used to do comparisons (though it can't be used to get
 * a shared_ptr via lock() anymore).
 */

template<typename DerivedType>
class FamilyTracking : public ekat::enable_shared_from_this<DerivedType>
{
public:
  using derived_type  = DerivedType;
  using tracking_type = FamilyTracking<derived_type>;

  FamilyTracking ();
  FamilyTracking (const tracking_type&) = default;
  FamilyTracking& operator= (const tracking_type&) = default;
  ~FamilyTracking ();

  void create_parent_child_link (const std::shared_ptr<derived_type>& parent);

  std::shared_ptr<derived_type> get_parent () const { return m_parent; }

  const std::list<std::weak_ptr<derived_type>>& get_children () const { return m_children; }
protected:

  // Check if a weak_ptr points to the same object as this class
  bool is_same (const std::weak_ptr<derived_type>& src) const;

  std::shared_ptr<derived_type>            m_parent;
  std::list<std::weak_ptr<derived_type>>   m_children;
};

// =================== IMPLEMENTATION ====================== //

template<typename DerivedType>
FamilyTracking<DerivedType>::FamilyTracking ()
{
  // Note: we cannot put the static assert in the class decl, cause DerivedType
  //       is still incomplete at that point.
  static_assert (std::is_base_of<tracking_type,derived_type>::value,
      "Error! Do not instantiate FamilyTracking<T> if T does not inherit from FamilyTracking.\n"
      "       This class exploits the Curiously Recurring Template Pattern (CRTP).\n");
}

template<typename DerivedType>
FamilyTracking<DerivedType>::~FamilyTracking ()
{
  // If we are the child/parent of another instance, we need to
  // remove ourself as their parent/child respectively.

  // Since derived_type inherits from tracking_type, the pointer
  // returned by lock(), can be assigned to a ptr to tracking_type.
  if (m_parent) {
    // Scan the children of my parent, then remove myself.
    // NOTE: since we're in the dtor, the counter of the shared ptr
    //       has reached 0, so we cannot call shared_from_this.
    //       HOWEVER, the weak ptr still stores the same ptr,
    //       and can be compared with ohter weak_ptr's, by using
    //       a symmetric version of owner_before. In particular,
    //       a==b if !a.owner_before(b) and !b.owner_before(a).
    auto me = this->weak_from_this();
    auto& siblings = m_parent->m_children;
    bool found = false;
    for (auto it=siblings.begin(); it!=siblings.end(); ++it) {
      if (is_same(*it)) {
        found = true;
        siblings.erase(it);
        break;
      }
    }

    // Note: Cannot throw in a destructor, so just print and call std::abort
    if (not found) {
        printf("Error! Could not find this object in the list of the parent's children.\n"
               "       Aborting...\n");
        std::abort();
    }
  }

  // Remove myself as the parent of all my children
  for (auto it : m_children) {
    // Since derived_type inherits from tracking_type, the pointer
    // returned by lock(), can be assigned to a ptr to tracking_type.
    std::shared_ptr<tracking_type> c = it.lock();
    if (not c) {
      // Note: Cannot throw in a destructor, so just print and call std::abort
      printf("Error! One of the weak_ptr's of my children has expired.\n");
      std::abort();
    }

    c->m_parent.reset();
  }
}

template<typename DerivedType>
void FamilyTracking<DerivedType>::
create_parent_child_link (const std::shared_ptr<derived_type>& parent)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (this->shared_from_this(),
      "Error! Failure to get a shared object from *this.\n");
  EKAT_REQUIRE_MSG (m_parent==nullptr,
      "Error! This object already stores a parent.\n");

  auto me = this->weak_from_this ();
  EKAT_REQUIRE_MSG (me.lock(),"Error! Unable to aquire a shared_ptr to this object.\n");

  // Set parent
  m_parent = parent;

  // Safety check. This should never happen, but just in case
  for (auto it : parent->get_children()) {
    EKAT_REQUIRE_MSG (not is_same(it),
        "Error! This object is already in the list of children of the input parent.\n");
  }

  // Add myself as child in my parent's list
  parent->m_children.push_back(me);
}

template<typename DerivedType>
bool FamilyTracking<DerivedType>::
is_same (const std::weak_ptr<derived_type>& src) const
{
  auto me = this->weak_from_this();
  return not src.owner_before(me) and not me.owner_before(src);
}

} // namespace scream

#endif // SCREAM_FAMILY_TRACKING_CLASS
