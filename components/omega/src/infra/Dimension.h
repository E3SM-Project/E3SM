#ifndef OMEGA_DIMENSION_H
#define OMEGA_DIMENSION_H
//===-- infra/Dimension.h - OMEGA dimension class ----------------*- C++
//-*-===//
//
/// \file
/// \brief Defines a dimension class used by the multi-dimensional Field class
///
/// The Dimension class defines dimension information for array fields.
/// Dimensions can be either distributed (most horizontal dimensions) or
/// non-distributed (eg vertical). The dimension primarily defines lengths
/// and offsets that are used for parallel IO and as dimension metadata in
/// self-describing files.
///
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Logging.h"
#include <map>
#include <memory>
#include <set>

namespace OMEGA {

//------------------------------------------------------------------------------
/// The Dimension class manages dimension information for fields that are needed
/// for parallel IO and metadata
class Dimension {
 private:
   /// Store and maintain all defined dimensions
   static std::map<std::string, std::shared_ptr<Dimension>> AllDims;

   /// Name of dimension
   std::string DimName;

   /// Flag to identify whether this is a distributed dimension or not
   bool Distributed;

   /// Global length (length of full unpartitioned dimension)
   I4 GlobalLength;

   /// Local length (size of local partition for this dim)
   I4 LocalLength;

   /// Offset for parallel IO - this is the memory offset in the global
   /// index space for each local address along this dimension
   /// This is typically the 0-based global ID (eg CellID - 1).
   HostArray1DI4 Offset;

 public:
   //---------------------------------------------------------------------------
   /// Creates a distributed dimension given a name, global (unpartitioned)
   /// length, size of local partition (including ghost cells and padding),
   /// and GlobalID/index for each local entry. Entries that should not be
   /// used in IO (eg ghost cells) should be given a GlobalID of -1.
   /// The created dimension will be added to the list of defined dimensions.
   /// If a dimension has already been defined for a dimension of the same
   /// name and global/local length, the existing dimension will be returned.
   static std::shared_ptr<Dimension>
   create(const std::string &Name, ///< [in] name of dimension
          const I4 GlobalLength,   ///< [in] global (unpartitioned) size of dim
          const I4 LocalLength,    ///< [in] size of dim in local partition
          HostArray1DI4 Offset     ///< [in] glob indx offset for each local pt
   );

   //---------------------------------------------------------------------------
   /// Creates a non-distributed dimension given a name and length
   static std::shared_ptr<Dimension>
   create(const std::string &Name, ///< [in] name of dimension
          const I4 GlobalLength   ///< [in] length of dimension
   );

   //---------------------------------------------------------------------------
   // Checks to see if a dim of this name exists
   static bool
   exists(const std::string &Name ///< [in] name of dimension to check
   );

   //---------------------------------------------------------------------------
   // Destroys a dimension
   static void
   destroy(const std::string &Name ///< [in] name of dimension to destroy
   );

   //---------------------------------------------------------------------------
   // Destroys all defined dimensions
   static void clear();

   //---------------------------------------------------------------------------
   // Retrieves full dimension instance by name
   static std::shared_ptr<Dimension>
   get(const std::string &Name // [in] Name of dimension
   );

   //---------------------------------------------------------------------------
   /// Get name of dimension from instance (iterator)
   std::string getName() const;

   //---------------------------------------------------------------------------
   /// Check distributed property of a dimension from instance
   bool isDistributed() const;

   //---------------------------------------------------------------------------
   /// Check distributed property of dimension given its name
   static bool
   isDistributedDim(const std::string &Name ///< [in] name of dimension
   );

   //---------------------------------------------------------------------------
   /// Get global dimension from instance
   I4 getLengthGlobal() const;

   //---------------------------------------------------------------------------
   /// Get global dimension length by name
   static I4
   getDimLengthGlobal(const std::string &Name ///< [in] name of dimension
   );

   //---------------------------------------------------------------------------
   /// Get length of dimension in local partition from instance
   I4 getLengthLocal() const;

   //---------------------------------------------------------------------------
   /// Get length of dimension in local partition by name
   static I4
   getDimLengthLocal(const std::string &Name ///< [in] name of dimension
   );

   //---------------------------------------------------------------------------
   /// Get global offset for each local address from dim instance
   HostArray1DI4 getOffset() const;

   //---------------------------------------------------------------------------
   /// Get global offset for each local address given a dim name
   static HostArray1DI4
   getDimOffset(const std::string &Name ///< [in] name of dimension
   );

   //----------------------------------------------------------------------------//
   // Retrieves the total number of currently defined dimensions
   static int getNumDefinedDims();

   //---------------------------------------------------------------------------
   /// An iterator can be used to loop through all defined dimensions
   using Iter =
       std::map<std::string, std::shared_ptr<Dimension>>::iterator;

   /// Returns an iterator to the first dimension stored in AllDims
   static Iter begin();

   /// Returns an iterator to the last dimension stored in AllDims
   static Iter end();

}; // end class Dimension

} // namespace OMEGA

#endif // OMEGA_DIMENSION_H
