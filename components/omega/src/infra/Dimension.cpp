//===-- infra/Dimension.cpp - OMEGA dimension implementation -----*- C++
//-*-===//
//
// Implementation of the field dimension class used by the multi-dimensional
// Field class
//
// The Dimension class defines dimension information for array fields.
// Dimensions can be either distributed (most horizontal dimensions) or
// non-distributed (eg vertical). The dimension primarily defines lengths
// and offsets that are used for parallel IO and as dimension metadata in
// self-describing files.
//
//===----------------------------------------------------------------------===//

#include "Dimension.h"
#include "DataTypes.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include <map>
#include <memory>
#include <string>

namespace OMEGA {

// Create static class members that store instantiations of metadata
std::map<std::string, std::shared_ptr<Dimension>> Dimension::AllDims;

//------------------------------------------------------------------------------
// Creates a distributed dimension given a name, global (unpartitioned)
// length, size of local partition (including ghost cells and padding),
// and GlobalID/index for each local entry. Entries that should not be
// used in IO (eg ghost cells) should be given a GlobalID of -1.
// The created dimension will be added to the list of defined dimensions.
// If a dimension has already been defined for a dimension of the same
// name and global/local length, the existing dimension will be returned.
std::shared_ptr<Dimension> Dimension::create(
    const std::string &Name, // [in] name of dimension
    const I4 GlobalLength,   // [in] global (unpartitioned) size of dim
    const I4 LocalLength,    // [in] size of dim in local partition
    HostArray1DI4 Offset     // [in] glob indx offset for each local indx
) {

   auto Dim = std::make_shared<Dimension>(); // create empty dim

   if (exists(Name)) { // Dimension already exists

      // Retrieve the previously defined dim and use that unless...
      Dim = get(Name);

      // The already-defined dim has different properties, so the dimensions
      // are not the same and we have a conflict
      if ((Dim->GlobalLength != GlobalLength) or
          (Dim->LocalLength != LocalLength) or !Dim->Distributed) {
         LOG_ERROR("Attempt to create dimension {} but a dimension with"
                   " that name already exists with different properties",
                   Name);
         Dim = nullptr;
      }

   } else { // a new dimension, fill data members

      Dim->DimName      = Name;
      Dim->GlobalLength = GlobalLength;
      Dim->LocalLength  = LocalLength;
      Dim->Offset       = Offset;
      Dim->Distributed  = true;
      AllDims[Name]     = Dim; // add to list of dims
   }

   return Dim;

} // end create distributed

//------------------------------------------------------------------------------
// Creates a non-distributed dimension given a name and length
std::shared_ptr<Dimension>
Dimension::create(const std::string &Name, // [in] name of dimension
                  const I4 GlobalLength    // [in] length of dimension
) {
   auto Dim = std::make_shared<Dimension>(); // create empty dim

   if (exists(Name)) { // Dimension already exists

      // Retrieve the previously defined dim and use that unless...
      Dim = get(Name);

      // The already-defined dim has different properties, so the dimensions
      // are not the same and we have a conflict
      if ((Dim->GlobalLength != GlobalLength) or
          (Dim->LocalLength != GlobalLength) or Dim->Distributed) {
         LOG_ERROR("Attempt to create dimension {} but a dimension with"
                   " that name already exists with different properties",
                   Name);
         Dim = nullptr;
      }

   } else { // a new dimension, fill data members

      Dim->DimName           = Name;
      Dim->GlobalLength      = GlobalLength;
      Dim->LocalLength       = GlobalLength;
      Dim->Distributed       = false;
      std::string OffsetName = Name + "Offset";
      HostArray1DI4 NewOffset(OffsetName, GlobalLength);
      for (int I = 0; I < GlobalLength; ++I) {
         NewOffset(I) = I;
      }
      Dim->Offset   = NewOffset;
      AllDims[Name] = Dim; // add to list of dims
   }

   return Dim;

} // end create non-distributed

//------------------------------------------------------------------------------
// Checks to see if a dim of this name exists
bool Dimension::exists(
    const std::string &Name // [in] name of dimension to check
) {
   return (AllDims.find(Name) != AllDims.end());
} // end exists

//------------------------------------------------------------------------------
// Destroys a dimension
void Dimension::destroy(
    const std::string &Name // [in] name of dimension to destroy
) {
   if (exists(Name)) {
      AllDims.erase(Name);
   } else {
      LOG_ERROR("Attempt to destroy the dimension {} failed: "
                "dimension does not exist or has not been defined.",
                Name);
   }

} // end destroy

//------------------------------------------------------------------------------
// Destroys all defined dimensions
void Dimension::clear() { AllDims.clear(); }

//----------------------------------------------------------------------------//
// Retrieves full dimension instance by name
std::shared_ptr<Dimension>
Dimension::get(const std::string &Name // [in] Name of dimension
) {
   if (exists(Name)) {
      return AllDims[Name];
   } else {
      LOG_ERROR("Cannot retrieve dimension {}: dimension does not exist"
                " or has not net been defined",
                Name);
      return nullptr;
   }

} // end get full dimension instance

//------------------------------------------------------------------------------
// Get name of dimension from instance (iterator)
std::string Dimension::getName() const { return DimName; }

//------------------------------------------------------------------------------
// Check whether dimension is distributed
bool Dimension::isDistributed() const { return Distributed; }

//------------------------------------------------------------------------------
// Check whether dimension is distributed by dimension name
bool Dimension::isDistributedDim(
    const std::string &Name // [in] name of dimension
) {
   // Make sure dimension exists
   if (exists(Name)) {
      std::shared_ptr<Dimension> ThisDim = AllDims[Name];
      return ThisDim->Distributed;

   } else {
      LOG_ERROR("Cannot check distribution of dimension {}: "
                "dimension does not exist or has not been defined",
                Name);
      return false;
   }
}

//------------------------------------------------------------------------------
// Get global dimension from instance
I4 Dimension::getLengthGlobal() const { return GlobalLength; }

//------------------------------------------------------------------------------
// Get global dimension length by name
I4 Dimension::getDimLengthGlobal(
    const std::string &Name // [in] name of dimension
) {
   I4 Length;

   // Make sure dimension exists
   if (exists(Name)) {
      std::shared_ptr<Dimension> ThisDim = AllDims[Name];
      Length                             = ThisDim->GlobalLength;

   } else {
      LOG_ERROR("Cannot get global length of dimension {}: "
                "dimension does not exist or has not been defined",
                Name);
      Length = -1;
   }

   return Length;

} // end getDimLengthGlobal

//------------------------------------------------------------------------------
// Get length dimension in local partition from instance
I4 Dimension::getLengthLocal() const { return LocalLength; }

//------------------------------------------------------------------------------
// Get length dimension in local partition by name
I4 Dimension::getDimLengthLocal(
    const std::string &Name // [in] name of dimension
) {
   I4 Length;

   // Make sure dimension exists
   if (exists(Name)) {
      std::shared_ptr<Dimension> ThisDim = AllDims[Name];
      Length                             = ThisDim->LocalLength;

   } else {
      LOG_ERROR("Cannot get local length of dimension {}: "
                "dimension does not exist or has not been defined",
                Name);
      Length = -1;
   }

   return Length;

} // end get DimLengthLocal

//------------------------------------------------------------------------------
// Get global offset for each local address from dim instance
HostArray1DI4 Dimension::getOffset() const { return Offset; }

//------------------------------------------------------------------------------
// Get global offset for each local address given a dim name
HostArray1DI4
Dimension::getDimOffset(const std::string &Name // [in] name of dimension
) {

   // Make sure dimension exists
   if (exists(Name)) {
      std::shared_ptr<Dimension> ThisDim = AllDims[Name];
      return ThisDim->Offset;

   } else {
      LOG_ERROR("Cannot get offset array for dimension {}: "
                "dimension does not exist or has not been defined",
                Name);
   }

} // end getDimOffset

//----------------------------------------------------------------------------//
// Retrieves the number of currently defined dimensions
int Dimension::getNumDefinedDims() { return AllDims.size(); }

//------------------------------------------------------------------------------
// Returns an iterator to the first dimension stored in AllDims
Dimension::Iter Dimension::begin() { return Dimension::AllDims.begin(); }

//------------------------------------------------------------------------------
// Returns an iterator to the last dimension stored in AllDims
Dimension::Iter Dimension::end() { return Dimension::AllDims.end(); }

} // namespace OMEGA
