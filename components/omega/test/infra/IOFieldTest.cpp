//===-- Test driver for OMEGA IO Fields --------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA IO Fields
///
/// This driver tests the capabilities for OMEGA to manage the registration
/// and use of fields for IO. This test driver will not actually include IO
/// but only tests the registration and retrieval of fields.
//
//===-----------------------------------------------------------------------===/

#include "IOField.h"
#include "DataTypes.h"
#include "Logging.h"
#include "MetaData.h"
#include "OmegaKokkos.h"
#include "mpi.h"
#include <cstdlib>
#include <vector>

//------------------------------------------------------------------------------
// Initialization routine to create reference IO Fields
int initIOFieldTest() {

   int Err = 0;

   // Define dimensions
   OMEGA::I4 NCellsSize  = 100;
   OMEGA::I4 NVertLevels = 64;
   auto CellDim          = OMEGA::MetaDim::create("NCells", NCellsSize);
   auto VertDim          = OMEGA::MetaDim::create("NVertLevels", NVertLevels);
   std::vector<std::shared_ptr<OMEGA::MetaDim>> Dimensions{CellDim, VertDim};

   // Define field metadata for all fields
   auto FieldMetaI4H = OMEGA::ArrayMetaData::create(
       "FieldI4H",
       "Test IO Field for I4 Host", /// long Name or description
       "unitless",                  /// units
       "I4StdName",                 /// CF standard Name
       0,                           /// min valid value
       100000,                      /// max valid value
       -999,                        /// scalar used for undefined entries
       2,                           /// number of dimensions
       Dimensions                   /// dim pointers
   );

   auto FieldMetaI4D = OMEGA::ArrayMetaData::create(
       "FieldI4D",
       "Test IO Field for I4 Device", /// long Name or description
       "unitless",                    /// units
       "I4StdName",                   /// CF standard Name
       0,                             /// min valid value
       100000,                        /// max valid value
       -999,                          /// scalar used for undefined entries
       2,                             /// number of dimensions
       Dimensions                     /// dim pointers
   );

   auto FieldMetaR8H = OMEGA::ArrayMetaData::create(
       "FieldR8H",
       "Test IO Field for R8 Host", /// long Name or description
       "m",                         /// units
       "R8StdName",                 /// CF standard Name
       -12345.0,                    /// min valid value
       12345.0,                     /// max valid value
       -9.99E+30,                   /// scalar used for undefined entries
       2,                           /// number of dimensions
       Dimensions                   /// dim pointers
   );

   auto FieldMetaR8D = OMEGA::ArrayMetaData::create(
       "FieldR8D",
       "Test IO Field for R8 Device", /// long Name or description
       "m",                           /// units
       "R8StdName",                   /// CF standard Name
       -12345.0,                      /// min valid value
       12345.0,                       /// max valid value
       -9.99E+30,                     /// scalar used for undefined entries
       2,                             /// number of dimensions
       Dimensions                     /// dim pointers
   );

   // Create host data arrays
   OMEGA::HostArray2DI4 DataI4H("FieldI4H", NCellsSize, NVertLevels);
   OMEGA::HostArray2DR8 DataR8H("FieldR8H", NCellsSize, NVertLevels);
   for (int Cell = 0; Cell < NCellsSize; ++Cell) {
      for (int k = 0; k < NVertLevels; ++k) {
         DataI4H(Cell, k) = Cell + k;
         DataR8H(Cell, k) = Cell + k + 1.2345678;
      }
   }

   // Create device data arrays
   OMEGA::Array2DI4 DataI4D("FieldI4D", NCellsSize, NVertLevels);
   OMEGA::Array2DR8 DataR8D("FieldR8D", NCellsSize, NVertLevels);
   OMEGA::parallelFor(
       {NCellsSize, NVertLevels}, KOKKOS_LAMBDA(int Cell, int k) {
          DataI4D(Cell, k) = Cell + k + 1;
          DataR8D(Cell, k) = Cell + k + 2.2345678;
       });

   // Define IO Fields
   int Err1 = OMEGA::IOField::define("FieldI4H");
   int Err2 = OMEGA::IOField::define("FieldI4D");
   if (Err1 == 0 && Err2 == 0) {
      LOG_INFO("IOField: initializing I4 field: PASS");
   } else {
      LOG_ERROR("IOField: initializing I4 field: FAIL");
      Err += std::abs(Err1) + std::abs(Err2);
   }

   Err1 = OMEGA::IOField::attachData<OMEGA::HostArray2DI4>("FieldI4H", DataI4H);
   Err2 = OMEGA::IOField::attachData<OMEGA::Array2DI4>("FieldI4D", DataI4D);
   if (Err1 == 0 && Err2 == 0) {
      LOG_INFO("IOField: attaching I4 data: PASS");
   } else {
      LOG_ERROR("IOField: attaching I4 data: FAIL");
      Err += std::abs(Err1) + std::abs(Err2);
   }

   Err1 = OMEGA::IOField::define("FieldR8H");
   Err2 = OMEGA::IOField::define("FieldR8D");
   if (Err1 == 0 && Err2 == 0) {
      LOG_INFO("IOField: initializing R8 field: PASS");
   } else {
      LOG_ERROR("IOField: initializing R8 field: FAIL");
      Err += std::abs(Err1) + std::abs(Err2);
   }

   Err1 = OMEGA::IOField::attachData<OMEGA::HostArray2DR8>("FieldR8H", DataR8H);
   Err2 = OMEGA::IOField::attachData<OMEGA::Array2DR8>("FieldR8D", DataR8D);
   if (Err1 == 0 && Err2 == 0) {
      LOG_INFO("IOField: attaching R8 data: PASS");
   } else {
      LOG_ERROR("IOField: attaching R8 data: FAIL");
      Err += std::abs(Err1) + std::abs(Err2);
   }

   return Err;

} // End initialization of IO Fields

//------------------------------------------------------------------------------
// We will test the IO Field interfaces by defining an IO Field during init
// and then retrieving the field and comparing results. Because IO Field makes
// use of the Metadata class and mostly points to other arrays/classes, we
// do not test all possible combinations - just two data types and host/device
// arrays to make sure generic pointers work fine.

int main(int argc, char **argv) {

   int Err  = 0;
   int Err1 = 0;
   int Err2 = 0;
   int Err3 = 0;
   int Err4 = 0;

   // Initialize the global MPI environment
   // We do not actually use message passing but need to test the
   // array types and behavior within the distributed environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {
      // Call initialization to create reference IO field
      Err = initIOFieldTest();
      if (Err != 0)
         LOG_ERROR("IOFieldTest: Error in initialization routine");

      // Set reference data - must match the values in the init routine
      std::string RefIUnits = "unitless";
      std::string RefRUnits = "m";
      int NCellsSize        = 100;
      int NVertLevels       = 64;
      OMEGA::HostArray2DI4 RefI4H("RefI4H", NCellsSize, NVertLevels);
      OMEGA::HostArray2DR8 RefR8H("RefR8H", NCellsSize, NVertLevels);
      for (int Cell = 0; Cell < NCellsSize; ++Cell) {
         for (int k = 0; k < NVertLevels; ++k) {
            RefI4H(Cell, k) = Cell + k;
            RefR8H(Cell, k) = Cell + k + 1.2345678;
         }
      }
      OMEGA::Array2DI4 RefI4D("RefI4D", NCellsSize, NVertLevels);
      OMEGA::Array2DR8 RefR8D("RefR8D", NCellsSize, NVertLevels);
      OMEGA::parallelFor(
          {NCellsSize, NVertLevels}, KOKKOS_LAMBDA(int Cell, int k) {
             RefI4D(Cell, k) = Cell + k + 1;
             RefR8D(Cell, k) = Cell + k + 2.2345678;
          });

      // Check existence of fields
      bool FieldExistsI4H = OMEGA::IOField::isDefined("FieldI4H");
      bool FieldExistsI4D = OMEGA::IOField::isDefined("FieldI4D");
      bool FieldExistsR8H = OMEGA::IOField::isDefined("FieldR8H");
      bool FieldExistsR8D = OMEGA::IOField::isDefined("FieldR8D");
      if (FieldExistsI4H && FieldExistsI4D && FieldExistsR8H &&
          FieldExistsR8D) {
         LOG_INFO("IOFieldTest: existence test PASS");
      } else {
         LOG_ERROR("IOFieldTest: existence test FAIL");
      }

      bool FieldExistsJunk = OMEGA::IOField::isDefined("FieldJunk");
      if (!FieldExistsJunk) {
         LOG_INFO("IOFieldTest: non-existence test PASS");
      } else {
         LOG_ERROR("IOFieldTest: non-existence test FAIL");
      }

      // Test retrieval of data and metadata
      // Retrieve metadata first
      std::shared_ptr<OMEGA::MetaData> MetaI4D =
          OMEGA::IOField::getMetaData("FieldI4D");
      std::shared_ptr<OMEGA::MetaData> MetaR8D =
          OMEGA::IOField::getMetaData("FieldR8D");
      std::shared_ptr<OMEGA::MetaData> MetaI4H =
          OMEGA::IOField::getMetaData("FieldI4H");
      std::shared_ptr<OMEGA::MetaData> MetaR8H =
          OMEGA::IOField::getMetaData("FieldR8H");

      std::string NewIUnits;
      std::string NewRUnits;

      Err1 = MetaI4H->getEntry("Units", NewIUnits);
      if (Err1 == 0 && NewIUnits == RefIUnits) {
         LOG_INFO("IOField: Retrieve I4H metadata by name: PASS");
      } else {
         LOG_ERROR("IOField: Retrieve I4H metadata by name: FAIL");
      }

      Err2 = MetaR8H->getEntry("Units", NewRUnits);
      if (Err2 == 0 && NewRUnits == RefRUnits) {
         LOG_INFO("IOField: Retrieve R8H metadata by name: PASS");
      } else {
         LOG_ERROR("IOField: Retrieve R8H metadata by name: FAIL");
      }

      Err3 = MetaI4D->getEntry("Units", NewIUnits);
      if (Err3 == 0 && NewIUnits == RefIUnits) {
         LOG_INFO("IOField: Retrieve I4D metadata by name: PASS");
      } else {
         LOG_ERROR("IOField: Retrieve I4D metadata by name: FAIL");
      }

      Err4 = MetaR8D->getEntry("Units", NewRUnits);
      if (Err4 == 0 && NewRUnits == RefRUnits) {
         LOG_INFO("IOField: Retrieve R8D metadata by name: PASS");
      } else {
         LOG_ERROR("IOField: Retrieve R8D metadata by name: FAIL");
      }
      Err += std::abs(Err1) + std::abs(Err2) + std::abs(Err3) + std::abs(Err4);

      // Now retrieve full data
      OMEGA::HostArray2DI4 NewI4H =
          OMEGA::IOField::getData<OMEGA::HostArray2DI4>("FieldI4H");
      OMEGA::HostArray2DR8 NewR8H =
          OMEGA::IOField::getData<OMEGA::HostArray2DR8>("FieldR8H");
      OMEGA::Array2DI4 NewI4D =
          OMEGA::IOField::getData<OMEGA::Array2DI4>("FieldI4D");
      OMEGA::Array2DR8 NewR8D =
          OMEGA::IOField::getData<OMEGA::Array2DR8>("FieldR8D");

      Err1 = 0;
      Err2 = 0;

      for (int Cell = 0; Cell < NCellsSize; ++Cell) {
         for (int k = 0; k < NVertLevels; ++k) {
            if (NewI4H(Cell, k) != RefI4H(Cell, k))
               ++Err1;
            if (NewR8H(Cell, k) != RefR8H(Cell, k))
               ++Err2;
         }
      }
      if (Err1 == 0) {
         LOG_INFO("IOField: Retrieve I4 host data by name: PASS");
      } else {
         LOG_ERROR("IOField: Retrieve I4 host data by name: FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOField: Retrieve R8 host data by name: PASS");
      } else {
         LOG_ERROR("IOField: Retrieve R8 host data by name: FAIL");
      }

      OMEGA::Array2DI4 ErrArray("ErrorArray", NCellsSize, NVertLevels);
      OMEGA::parallelFor(
          {NCellsSize, NVertLevels}, KOKKOS_LAMBDA(int Cell, int k) {
             if ((NewI4D(Cell, k) != RefI4D(Cell, k)) or
                 (NewR8D(Cell, k) != RefR8D(Cell, k)))
                ErrArray(Cell, k) += 1;
          });

      // Create a reducer
      Kokkos::Sum<OMEGA::I4> reducer(Err3);

      // Perform the reduction
      OMEGA::parallelReduce(
          "SumReduce", {NCellsSize, NVertLevels},
          KOKKOS_LAMBDA(int Cell, int k, OMEGA::I4 &update) {
             update += ErrArray(Cell, k);
          },
          reducer);

      if (Err3 == 0) {
         LOG_INFO("IOField: Retrieve device data by name: PASS");
      } else {
         LOG_ERROR("IOField: Retrieve device data by name: FAIL");
      }

      Err += std::abs(Err1) + std::abs(Err2) + std::abs(Err3);

      // Erase a field and check for non-existence
      OMEGA::IOField::erase("FieldI4D");
      FieldExistsI4D = OMEGA::IOField::isDefined("FieldI4D");
      if (!FieldExistsI4D) {
         LOG_INFO("IOFieldTest: erase field PASS");
      } else {
         LOG_ERROR("IOFieldTest: erase field FAIL");
         ++Err;
      }

      // Clear all fields
      OMEGA::IOField::clear();
      FieldExistsI4H = OMEGA::IOField::isDefined("FieldI4H");
      FieldExistsI4D = OMEGA::IOField::isDefined("FieldI4D");
      FieldExistsR8H = OMEGA::IOField::isDefined("FieldR8H");
      FieldExistsR8D = OMEGA::IOField::isDefined("FieldR8D");
      if (FieldExistsI4H or FieldExistsI4D or FieldExistsR8H or
          FieldExistsR8D) {
         LOG_ERROR("IOFieldTest: clear all data FAIL");
         ++Err;
      } else {
         LOG_INFO("IOFieldTest: clear all data PASS");
      }
   }
   Kokkos::finalize();
   MPI_Finalize();

   // End of testing
   return Err;
}
//===--- End test driver for IO Field --------------------------------------===/
