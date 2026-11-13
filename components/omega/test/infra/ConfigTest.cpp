//===-- Test driver for OMEGA Config -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Config
///
/// This driver tests the OMEGA Configuration module that reads and distributes
/// the model configuration from a YAML input file. It first creates a YAML
/// node, writes to a sample input file, the reads that file to set sample
/// configuration variables.
///
//
//===-----------------------------------------------------------------------===/

#include "Config.h"
#include "Broadcast.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "mpi.h"

#include <iostream>
#include <string>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The test driver for Config. This creates a Configuration in YAML format,
// writes that configuration to a file and then reads it back in, verifying
// that values are the same.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);

   // Initialize the Machine Environment and retrieve the default environment
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv    = OMEGA::MachEnv::getDefault();
   MPI_Comm OmegaComm = DefEnv->getComm();

   // Initialize the Logging and timing systems
   initLogging(DefEnv);
   Pacer::initialize(OmegaComm);
   Pacer::setPrefix("Omega:");
   LOG_INFO("------ Config unit tests ------");

   // Define some variables for a reference configuration
   // These are meant to create a hierachy that tests at least two levels
   // and includes a variable of each supported type.
   // Omega:
   //   Hmix:
   //     HmixOn: true
   //     HmixI4: 3
   //     HmixI8: 123456789
   //     HmixR4: 5.0
   //     HmixR8: 1.234567890
   //     HmixStr: "HmixString"
   //     HmixDel2:
   //       HmixDel2On: true
   //       HmixDel2I4: 5
   //       HmixDel2I8: 999999999
   //       HmixDel2R4: 7.0
   //       HmixDel2R8: 9.999999999
   //       HmixDel2Str: "HmixDel2String"
   //   Vmix:
   //     VmixOn: false
   //     VmixI4: 4
   //     VmixI8: 987654321
   //     VmixR4: 6.0
   //     VmixR8: 9.9876543210;
   //     VmixStr: "VmixString"
   //   VectorI4: [1, 2, 3, 4, 5]
   //   VectorI8: [123456789, 234567890, 345678901, 456789012, 567890123]
   //   VectorR4: [1.2345, 2.3456, 3.4567, 4.5678, 5.6789]
   //   VectorR8: [1.23456789, 2.34567890, 3.45678901, 4.56789012, 5.67890123]
   //   VectorLog: [true, false, true, false, true]
   //   StrList:
   //     - first
   //     - second
   //     - third
   //     - fourth
   //     - fifth
   //

   bool RefHmixOn             = true;
   I4 RefHmixI4               = 3;
   I8 RefHmixI8               = 123456789;
   R4 RefHmixR4               = 5.0;
   R8 RefHmixR8               = 1.234567890;
   std::string RefHmixStr     = "HmixString";
   bool RefHmixDel2On         = true;
   I4 RefHmixDel2I4           = 5;
   I8 RefHmixDel2I8           = 999999999;
   R4 RefHmixDel2R4           = 7.0;
   R8 RefHmixDel2R8           = 9.999999999;
   std::string RefHmixDel2Str = "HmixDel2String";
   bool RefVmixOn             = false;
   I4 RefVmixI4               = 4;
   I8 RefVmixI8               = 987654321;
   R4 RefVmixR4               = 6.0;
   R8 RefVmixR8               = 9.9876543210;
   std::string RefVmixStr     = "VmixString";
   int VecSize                = 5;
   std::vector<I4> RefVecI4{1, 2, 3, 4, 5};
   std::vector<I8> RefVecI8{123456789, 234567890, 345678901, 456789012,
                            567890123};
   std::vector<R4> RefVecR4{1.2345, 2.3456, 3.4567, 4.5678, 5.6789};
   std::vector<R8> RefVecR8{1.23456789, 2.34567890, 3.45678901, 4.56789012,
                            5.67890123};
   std::vector<bool> RefVecLog{true, false, true, false, true};
   std::vector<std::string> RefList{"first", "second", "third", "fourth",
                                    "fifth"};

   // Build up a reference configuration
   Config ConfigOmegaAll("Omegaroot");
   Config ConfigOmegaRef("Omega");
   Config ConfigHmixRef("Hmix");
   Config ConfigHmixDel2Ref("HmixDel2");
   Config ConfigVmixRef("Vmix");

   // Test adding various types to config. These all abort on error so
   // making it through them all is a PASS.
   ConfigHmixDel2Ref.add("HmixDel2On", RefHmixDel2On);
   ConfigHmixDel2Ref.add("HmixDel2I4", RefHmixDel2I4);
   ConfigHmixDel2Ref.add("HmixDel2I8", RefHmixDel2I8);
   ConfigHmixDel2Ref.add("HmixDel2R4", RefHmixDel2R4);
   ConfigHmixDel2Ref.add("HmixDel2R8", RefHmixDel2R8);
   ConfigHmixDel2Ref.add("HmixDel2Str", RefHmixDel2Str);

   ConfigHmixRef.add("HmixOn", RefHmixOn);
   ConfigHmixRef.add("HmixI4", RefHmixI4);
   ConfigHmixRef.add("HmixI8", RefHmixI8);
   ConfigHmixRef.add("HmixR4", RefHmixR4);
   ConfigHmixRef.add("HmixR8", RefHmixR8);
   ConfigHmixRef.add("HmixStr", RefHmixStr);
   ConfigHmixRef.add(ConfigHmixDel2Ref);

   ConfigVmixRef.add("VmixOn", RefVmixOn);
   ConfigVmixRef.add("VmixI4", RefVmixI4);
   ConfigVmixRef.add("VmixI8", RefVmixI8);
   ConfigVmixRef.add("VmixR4", RefVmixR4);
   ConfigVmixRef.add("VmixR8", RefVmixR8);
   ConfigVmixRef.add("VmixStr", RefVmixStr);

   ConfigOmegaRef.add(ConfigHmixRef);
   ConfigOmegaRef.add(ConfigVmixRef);

   ConfigOmegaRef.add("VectorI4", RefVecI4);
   ConfigOmegaRef.add("VectorI8", RefVecI8);
   ConfigOmegaRef.add("VectorR4", RefVecR4);
   ConfigOmegaRef.add("VectorR8", RefVecR8);
   ConfigOmegaRef.add("VectorLog", RefVecLog);
   ConfigOmegaRef.add("StrList", RefList);
   // create the full root node
   ConfigOmegaAll.add(ConfigOmegaRef);
   LOG_INFO("Config: Created test configuration using add functions");

   // Test retrievals by getting from the reference config and
   // checking against reference values.

   // Test logical retrievals for each config level
   // First extract sub-configurations from the top level config
   Config ConfigHmixNew("Hmix");
   Config ConfigVmixNew("Vmix");
   Config ConfigHmixDel2New("HmixDel2");

   Error Err; // default to success
   Err += ConfigOmegaRef.get(ConfigHmixNew);
   Err += ConfigOmegaRef.get(ConfigVmixNew);
   Err += ConfigHmixNew.get(ConfigHmixDel2New);
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of subconfigs - FAIL");

   // Test logical retrievals
   bool HmixOn;
   bool VmixOn;
   bool HmixDel2On;
   Err += ConfigHmixNew.get("HmixOn", HmixOn);
   Err += ConfigVmixNew.get("VmixOn", VmixOn);
   Err += ConfigHmixDel2New.get("HmixDel2On", HmixDel2On);
   if ((HmixOn != RefHmixOn) or (VmixOn != RefVmixOn) or
       (HmixDel2On != RefHmixDel2On))
      Err += Error(ErrorCode::Fail, "Retrieved bool values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of boolean vars - FAIL");

   // Test I4 retrievals
   I4 HmixI4     = 0;
   I4 VmixI4     = 0;
   I4 HmixDel2I4 = 0;
   Err += ConfigHmixNew.get("HmixI4", HmixI4);
   Err += ConfigVmixNew.get("VmixI4", VmixI4);
   Err += ConfigHmixDel2New.get("HmixDel2I4", HmixDel2I4);
   if ((HmixI4 != RefHmixI4) or (VmixI4 != RefVmixI4) or
       (HmixDel2I4 != RefHmixDel2I4))
      Err += Error(ErrorCode::Fail, "Retrieved I4 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of I4 vars - FAIL");

   // Test I8 retrievals
   I8 HmixI8     = 0;
   I8 VmixI8     = 0;
   I8 HmixDel2I8 = 0;
   Err += ConfigHmixNew.get("HmixI8", HmixI8);
   Err += ConfigVmixNew.get("VmixI8", VmixI8);
   Err += ConfigHmixDel2New.get("HmixDel2I8", HmixDel2I8);
   if ((HmixI8 != RefHmixI8) or (VmixI8 != RefVmixI8) or
       (HmixDel2I8 != RefHmixDel2I8))
      Err += Error(ErrorCode::Fail, "Retrieved I8 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of I8 vars - FAIL");

   // Test R4 retrievals
   R4 HmixR4     = 0;
   R4 VmixR4     = 0;
   R4 HmixDel2R4 = 0;
   Err += ConfigHmixNew.get("HmixR4", HmixR4);
   Err += ConfigVmixNew.get("VmixR4", VmixR4);
   Err += ConfigHmixDel2New.get("HmixDel2R4", HmixDel2R4);
   if ((HmixR4 != RefHmixR4) or (VmixR4 != RefVmixR4) or
       (HmixDel2R4 != RefHmixDel2R4))
      Err += Error(ErrorCode::Fail, "Retrieved R4 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of R4 vars - FAIL");

   // Test R8 retrievals
   R8 HmixR8     = 0;
   R8 VmixR8     = 0;
   R8 HmixDel2R8 = 0;
   Err += ConfigHmixNew.get("HmixR8", HmixR8);
   Err += ConfigVmixNew.get("VmixR8", VmixR8);
   Err += ConfigHmixDel2New.get("HmixDel2R8", HmixDel2R8);
   if ((HmixR8 != RefHmixR8) or (VmixR8 != RefVmixR8) or
       (HmixDel2R8 != RefHmixDel2R8))
      Err += Error(ErrorCode::Fail, "Retrieved R8 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of R8 vars - FAIL");

   // Test string retrievals
   std::string HmixStr;
   std::string VmixStr;
   std::string HmixDel2Str;
   Err += ConfigHmixNew.get("HmixStr", HmixStr);
   Err += ConfigVmixNew.get("VmixStr", VmixStr);
   Err += ConfigHmixDel2New.get("HmixDel2Str", HmixDel2Str);
   if ((HmixStr != RefHmixStr) or (VmixStr != RefVmixStr) or
       (HmixDel2Str != RefHmixDel2Str))
      Err += Error(ErrorCode::Fail, "Retrieved string values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of string vars - FAIL");

   // Test vector retrievals
   std::vector<I4> VecI4;
   Err += ConfigOmegaRef.get("VectorI4", VecI4);
   for (int i = 0; i < VecSize; ++i) {
      if (VecI4[i] != RefVecI4[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved I4 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of I4 vector - FAIL");

   std::vector<I8> VecI8;
   Err += ConfigOmegaRef.get("VectorI8", VecI8);
   for (int i = 0; i < VecSize; ++i) {
      if (VecI8[i] != RefVecI8[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved I8 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of I8 vector - FAIL");

   std::vector<R4> VecR4;
   Err += ConfigOmegaRef.get("VectorR4", VecR4);
   for (int i = 0; i < VecSize; ++i) {
      if (VecR4[i] != RefVecR4[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved R4 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of R4 vector - FAIL");

   std::vector<R8> VecR8;
   Err += ConfigOmegaRef.get("VectorR8", VecR8);
   for (int i = 0; i < VecSize; ++i) {
      if (VecR8[i] != RefVecR8[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved R8 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of R8 vector - FAIL");

   std::vector<bool> VecLog;
   Err += ConfigOmegaRef.get("VectorLog", VecLog);
   for (int i = 0; i < VecSize; ++i) {
      if (VecLog[i] != RefVecLog[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved boolean vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of boolean vector - FAIL");

   std::vector<std::string> NewList;
   Err += ConfigOmegaRef.get("StrList", NewList);
   for (int i = 0; i < VecSize; ++i) {
      if (NewList[i] != RefList[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved string list is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Add/Get of string list - FAIL");

   // Test changing values using the set functions
   bool NewHmixOn             = false;
   I4 NewHmixI4               = 4;
   I8 NewHmixI8               = 223456789;
   R4 NewHmixR4               = 6.0;
   R8 NewHmixR8               = 2.234567890;
   std::string NewHmixStr     = "HmixStringNew";
   bool NewHmixDel2On         = false;
   I4 NewHmixDel2I4           = 6;
   I8 NewHmixDel2I8           = 899999999;
   R4 NewHmixDel2R4           = 8.0;
   R8 NewHmixDel2R8           = 8.999999999;
   std::string NewHmixDel2Str = "HmixDel2StringNew";
   bool NewVmixOn             = true;
   I4 NewVmixI4               = 5;
   I8 NewVmixI8               = 887654321;
   R4 NewVmixR4               = 7.0;
   R8 NewVmixR8               = 8.9876543210;
   std::string NewVmixStr     = "VmixStringNew";

   ConfigHmixDel2New.set("HmixDel2On", NewHmixDel2On);
   ConfigHmixDel2New.set("HmixDel2I4", NewHmixDel2I4);
   ConfigHmixDel2New.set("HmixDel2I8", NewHmixDel2I8);
   ConfigHmixDel2New.set("HmixDel2R4", NewHmixDel2R4);
   ConfigHmixDel2New.set("HmixDel2R8", NewHmixDel2R8);
   ConfigHmixDel2New.set("HmixDel2Str", NewHmixDel2Str);

   ConfigHmixNew.set("HmixOn", NewHmixOn);
   ConfigHmixNew.set("HmixI4", NewHmixI4);
   ConfigHmixNew.set("HmixI8", NewHmixI8);
   ConfigHmixNew.set("HmixR4", NewHmixR4);
   ConfigHmixNew.set("HmixR8", NewHmixR8);
   ConfigHmixNew.set("HmixStr", NewHmixStr);

   ConfigVmixNew.set("VmixOn", NewVmixOn);
   ConfigVmixNew.set("VmixI4", NewVmixI4);
   ConfigVmixNew.set("VmixI8", NewVmixI8);
   ConfigVmixNew.set("VmixR4", NewVmixR4);
   ConfigVmixNew.set("VmixR8", NewVmixR8);
   ConfigVmixNew.set("VmixStr", NewVmixStr);
   LOG_INFO("Config: modified all config using set routines");

   // Test logical retrievals
   Err += ConfigHmixNew.get("HmixOn", HmixOn);
   Err += ConfigVmixNew.get("VmixOn", VmixOn);
   Err += ConfigHmixDel2New.get("HmixDel2On", HmixDel2On);
   if ((HmixOn != NewHmixOn) or (VmixOn != NewVmixOn) or
       (HmixDel2On != NewHmixDel2On))
      Err += Error(ErrorCode::Fail, "Retrieved logical values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of boolean vars - FAIL");

   // Test I4 retrievals
   Err += ConfigHmixNew.get("HmixI4", HmixI4);
   Err += ConfigVmixNew.get("VmixI4", VmixI4);
   Err += ConfigHmixDel2New.get("HmixDel2I4", HmixDel2I4);
   if ((HmixI4 != NewHmixI4) or (VmixI4 != NewVmixI4) or
       (HmixDel2I4 != NewHmixDel2I4))
      Err += Error(ErrorCode::Fail, "Retrieved I4 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of I4 vars - FAIL");

   // Test I8 retrievals
   Err += ConfigHmixNew.get("HmixI8", HmixI8);
   Err += ConfigVmixNew.get("VmixI8", VmixI8);
   Err += ConfigHmixDel2New.get("HmixDel2I8", HmixDel2I8);
   if ((HmixI8 != NewHmixI8) or (VmixI8 != NewVmixI8) or
       (HmixDel2I8 != NewHmixDel2I8))
      Err += Error(ErrorCode::Fail, "Retrieved I8 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of I8 vars - FAIL");

   // Test R4 retrievals
   Err += ConfigHmixNew.get("HmixR4", HmixR4);
   Err += ConfigVmixNew.get("VmixR4", VmixR4);
   Err += ConfigHmixDel2New.get("HmixDel2R4", HmixDel2R4);
   if ((HmixR4 != NewHmixR4) or (VmixR4 != NewVmixR4) or
       (HmixDel2R4 != NewHmixDel2R4))
      Err += Error(ErrorCode::Fail, "Retrieved R4 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of R4 vars - FAIL");

   // Test R8 retrievals
   Err += ConfigHmixNew.get("HmixR8", HmixR8);
   Err += ConfigVmixNew.get("VmixR8", VmixR8);
   Err += ConfigHmixDel2New.get("HmixDel2R8", HmixDel2R8);
   if ((HmixR8 != NewHmixR8) or (VmixR8 != NewVmixR8) or
       (HmixDel2R8 != NewHmixDel2R8))
      Err += Error(ErrorCode::Fail, "Retrieved R8 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of R8 vars - FAIL");

   // Test string retrievals
   Err += ConfigHmixNew.get("HmixStr", HmixStr);
   Err += ConfigVmixNew.get("VmixStr", VmixStr);
   Err += ConfigHmixDel2New.get("HmixDel2Str", HmixDel2Str);
   if ((HmixStr != NewHmixStr) or (VmixStr != NewVmixStr) or
       (HmixDel2Str != NewHmixDel2Str))
      Err += Error(ErrorCode::Fail, "Retrieved string values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of string vars - FAIL");

   // Test setting of vectors by modifying one entry, setting and retrieving
   RefVecI4[2]  = -1;
   RefVecI8[2]  = -123456789;
   RefVecR4[2]  = -1.2345;
   RefVecR8[2]  = -1.23456789;
   RefVecLog[2] = false;
   RefList[2]   = "junk";

   ConfigOmegaRef.set("VectorI4", RefVecI4);
   ConfigOmegaRef.set("VectorI8", RefVecI8);
   ConfigOmegaRef.set("VectorR4", RefVecR4);
   ConfigOmegaRef.set("VectorR8", RefVecR8);
   ConfigOmegaRef.set("VectorLog", RefVecLog);
   ConfigOmegaRef.set("StrList", RefList);
   LOG_INFO("Config: modified config vectors using set");

   // Test vector retrievals after set
   Err += ConfigOmegaRef.get("VectorI4", VecI4);
   for (int i = 0; i < VecSize; ++i) {
      if (VecI4[i] != RefVecI4[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved I4 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of I4 vector - FAIL");

   Err += ConfigOmegaRef.get("VectorI8", VecI8);
   for (int i = 0; i < VecSize; ++i) {
      if (VecI8[i] != RefVecI8[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved I8 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of I8 vector - FAIL");

   Err += ConfigOmegaRef.get("VectorR4", VecR4);
   for (int i = 0; i < VecSize; ++i) {
      if (VecR4[i] != RefVecR4[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved R4 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of R4 vector - FAIL");

   Err += ConfigOmegaRef.get("VectorR8", VecR8);
   for (int i = 0; i < VecSize; ++i) {
      if (VecR8[i] != RefVecR8[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved R8 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of R8 vector - FAIL");

   Err += ConfigOmegaRef.get("VectorLog", VecLog);
   for (int i = 0; i < VecSize; ++i) {
      if (VecLog[i] != RefVecLog[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved bool vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of bool vector - FAIL");

   Err += ConfigOmegaRef.get("StrList", NewList);
   for (int i = 0; i < VecSize; ++i) {
      if (NewList[i] != RefList[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved string list is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Set/Get of string list - FAIL");

   // Reset vectors back to original reference values
   RefVecI4[2]  = 3;
   RefVecI8[2]  = 345678901;
   RefVecR4[2]  = 3.4567;
   RefVecR8[2]  = 3.45678901;
   RefVecLog[2] = true;
   RefList[2]   = "third";
   ConfigOmegaRef.set("VectorI4", RefVecI4);
   ConfigOmegaRef.set("VectorI8", RefVecI8);
   ConfigOmegaRef.set("VectorR4", RefVecR4);
   ConfigOmegaRef.set("VectorR8", RefVecR8);
   ConfigOmegaRef.set("VectorLog", RefVecLog);
   ConfigOmegaRef.set("StrList", RefList);

   // Test error modes by trying to retrieve junk values
   Config ConfigJunk("junk");
   bool JunkOn;
   I4 JunkI4;
   I4 JunkI8;
   I4 JunkR4;
   I4 JunkR8;
   std::string JunkStr;

   Err += ConfigOmegaRef.get(ConfigJunk);
   CHECK_ERROR(Err, "Config: Expected error get non-existent config - PASS");
   Err += ConfigHmixRef.get("junk", JunkOn);
   CHECK_ERROR(Err, "Config: Expected error get non-existent bool - PASS");
   Err += ConfigHmixRef.get("junk", JunkI4);
   CHECK_ERROR(Err, "Config: Expected error get non-existent I4 - PASS");
   Err += ConfigHmixRef.get("junk", JunkI8);
   CHECK_ERROR(Err, "Config: Expected error get non-existent I8 - PASS");
   Err += ConfigHmixRef.get("junk", JunkR4);
   CHECK_ERROR(Err, "Config: Expected error get non-existent R4 - PASS");
   Err += ConfigHmixRef.get("junk", JunkR8);
   CHECK_ERROR(Err, "Config: Expected error get non-existent R8 - PASS");
   Err += ConfigHmixRef.get("junk", JunkStr);
   CHECK_ERROR(Err, "Config: Expected error get non-existent string - PASS");

   // Write the reference config to a file
   ConfigOmegaAll.write("omegaConfigTst.yml");

   // Read the environment back in
   Config::readAll("omegaConfigTst.yml");

   // Get the full configuration
   Config *ConfigOmega = Config::getOmegaConfig();

   // Retrieve all values and compare against reference values
   Config ConfigHmixOmega("Hmix");
   Config ConfigVmixOmega("Vmix");
   Config ConfigHmixDel2Omega("HmixDel2");

   Err += ConfigOmega->get(ConfigHmixOmega);
   Err += ConfigOmega->get(ConfigVmixOmega);
   Err += ConfigHmixOmega.get(ConfigHmixDel2Omega);
   CHECK_ERROR_ABORT(Err, "Config: retrieve subgroups from full config FAIL");

   Err += ConfigHmixOmega.get("HmixOn", HmixOn);
   Err += ConfigVmixOmega.get("VmixOn", VmixOn);
   Err += ConfigHmixDel2Omega.get("HmixDel2On", HmixDel2On);
   if ((HmixOn != NewHmixOn) or (VmixOn != NewVmixOn) or
       (HmixDel2On != NewHmixDel2On))
      Err += Error(ErrorCode::Fail, "Retrieved bool values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: retrieve bool from full config FAIL");

   Err += ConfigHmixOmega.get("HmixI4", HmixI4);
   Err += ConfigVmixOmega.get("VmixI4", VmixI4);
   Err += ConfigHmixDel2Omega.get("HmixDel2I4", HmixDel2I4);
   if ((HmixI4 != NewHmixI4) or (VmixI4 != NewVmixI4) or
       (HmixDel2I4 != NewHmixDel2I4))
      Err += Error(ErrorCode::Fail, "Retrieved I4 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: retrieve I4 from full config FAIL");

   Err += ConfigHmixOmega.get("HmixI8", HmixI8);
   Err += ConfigVmixOmega.get("VmixI8", VmixI8);
   Err += ConfigHmixDel2Omega.get("HmixDel2I8", HmixDel2I8);
   if ((HmixI8 != NewHmixI8) or (VmixI8 != NewVmixI8) or
       (HmixDel2I8 != NewHmixDel2I8))
      Err += Error(ErrorCode::Fail, "Retrieved I8 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: retrieve I8 from full config FAIL");

   Err += ConfigHmixOmega.get("HmixR4", HmixR4);
   Err += ConfigVmixOmega.get("VmixR4", VmixR4);
   Err += ConfigHmixDel2Omega.get("HmixDel2R4", HmixDel2R4);
   if ((HmixR4 != NewHmixR4) or (VmixR4 != NewVmixR4) or
       (HmixDel2R4 != NewHmixDel2R4))
      Err += Error(ErrorCode::Fail, "Retrieved R4 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: retrieve R4 from full config FAIL");

   Err += ConfigHmixOmega.get("HmixR8", HmixR8);
   Err += ConfigVmixOmega.get("VmixR8", VmixR8);
   Err += ConfigHmixDel2Omega.get("HmixDel2R8", HmixDel2R8);
   if ((HmixR8 != NewHmixR8) or (VmixR8 != NewVmixR8) or
       (HmixDel2R8 != NewHmixDel2R8))
      Err += Error(ErrorCode::Fail, "Retrieved R8 values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: retrieve R8 from full config FAIL");

   Err += ConfigHmixOmega.get("HmixStr", HmixStr);
   Err += ConfigVmixOmega.get("VmixStr", VmixStr);
   Err += ConfigHmixDel2Omega.get("HmixDel2Str", HmixDel2Str);
   if ((HmixStr != NewHmixStr) or (VmixStr != NewVmixStr) or
       (HmixDel2Str != NewHmixDel2Str))
      Err += Error(ErrorCode::Fail, "Retrieved string values are incorrect");
   CHECK_ERROR_ABORT(Err, "Config: retrieve string from full config FAIL");

   // Vector retrievals
   Err += ConfigOmega->get("VectorI4", VecI4);
   for (int i = 0; i < VecSize; ++i) {
      if (VecI4[i] != RefVecI4[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved I4 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Get I4 vector from full config - FAIL");

   Err += ConfigOmega->get("VectorI8", VecI8);
   for (int i = 0; i < VecSize; ++i) {
      if (VecI8[i] != RefVecI8[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved I8 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Get I8 vector from full config - FAIL");

   Err += ConfigOmega->get("VectorR4", VecR4);
   for (int i = 0; i < VecSize; ++i) {
      if (VecR4[i] != RefVecR4[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved R4 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Get R4 vector from full config - FAIL");

   Err += ConfigOmega->get("VectorR8", VecR8);
   for (int i = 0; i < VecSize; ++i) {
      if (VecR8[i] != RefVecR8[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved R8 vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Get R8 vector from full config - FAIL");

   Err += ConfigOmega->get("VectorLog", VecLog);
   for (int i = 0; i < VecSize; ++i) {
      if (VecLog[i] != RefVecLog[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved bool vector is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Get bool vector from full config - FAIL");

   Err += ConfigOmega->get("StrList", NewList);
   for (int i = 0; i < VecSize; ++i) {
      if (NewList[i] != RefList[i]) {
         Err += Error(ErrorCode::Fail, "Retrieved string list is incorrect");
         break;
      }
   }
   CHECK_ERROR_ABORT(Err, "Config: Get string list from full config - FAIL");

   // Test retrieval of values using an iterator
   int IList = 0;
   for (auto It = ConfigOmega->begin(); It != ConfigOmega->end(); ++It) {
      std::string NodeName;
      Err += Config::getName(It, NodeName);
      // Note that iteration order not guaranteed but so far has
      // been consistent with location in the file.
      if (IList == 0 and NodeName != "Hmix")
         Err += Error(ErrorCode::Fail, "Error retrieving Hmix with iterator");
      if (IList == 1 and NodeName != "Vmix")
         Err += Error(ErrorCode::Fail, "Error retrieving Vmix with iterator");
      if (IList == 2 and NodeName != "VectorI4")
         Err += Error(ErrorCode::Fail, "Error retrieving I4 vec with iterator");
      if (IList == 3 and NodeName != "VectorI8")
         Err += Error(ErrorCode::Fail, "Error retrieving I8 vec with iterator");
      if (IList == 4 and NodeName != "VectorR4")
         Err += Error(ErrorCode::Fail, "Error retrieving R4 vec with iterator");
      if (IList == 5 and NodeName != "VectorR8")
         Err += Error(ErrorCode::Fail, "Error retrieving R8 vec with iterator");
      if (IList == 6 and NodeName != "VectorLog")
         Err += Error(ErrorCode::Fail, "Error retrieving bool vec w/ iterator");
      if (IList == 7 and NodeName != "StrList")
         Err += Error(ErrorCode::Fail, "Error retrieving StrLIst w/ iterator");
      if (IList > 7)
         Err += Error(ErrorCode::Fail, "Too many items with iterator");
      ++IList;
   }
   CHECK_ERROR_ABORT(Err, "Config: Get config list using iter - FAIL");

   // Test removals by removing both variables and sub-configs and
   // checking the further retrievals fail

   ConfigVmixOmega.remove("VmixOn");
   ConfigVmixOmega.remove("VmixI4");
   ConfigVmixOmega.remove("VmixI8");
   ConfigVmixOmega.remove("VmixR4");
   ConfigVmixOmega.remove("VmixR8");
   ConfigVmixOmega.remove("VmixStr");

   // Retrieve a fresh copy to make sure the removal happens to the
   // entire Config, not just local copy. All variable retrievals should fail.

   Config ConfigVmixCheck("Vmix");
   Err += ConfigOmega->get(ConfigVmixCheck);
   CHECK_ERROR_ABORT(Err, "Config: error retrieving revised config - FAIL");
   Err += ConfigVmixCheck.get("VmixOn", NewVmixOn);
   CHECK_ERROR(Err, "Config: expected error retrieving removed bool - PASS");
   Err += ConfigVmixCheck.get("VmixI4", NewVmixI4);
   CHECK_ERROR(Err, "Config: expected error retrieving removed I4 - PASS");
   Err += ConfigVmixCheck.get("VmixI8", NewVmixI8);
   CHECK_ERROR(Err, "Config: expected error retrieving removed I8 - PASS");
   Err += ConfigVmixCheck.get("VmixR4", NewVmixR4);
   CHECK_ERROR(Err, "Config: expected error retrieving removed R4 - PASS");
   Err += ConfigVmixCheck.get("VmixR8", NewVmixR8);
   CHECK_ERROR(Err, "Config: expected error retrieving removed R8 - PASS");
   Err += ConfigVmixCheck.get("VmixStr", NewVmixStr);
   CHECK_ERROR(Err, "Config: expected error retrieving removed string - PASS");

   // Now try removing an entire sub-config
   ConfigHmixOmega.remove("HmixDel2"); // remove entire subconfig
   Config ConfigHmixCheck("Hmix");
   Config ConfigHmixDel2Check("HmixDel2");
   Err += ConfigOmega->get(ConfigHmixCheck);
   CHECK_ERROR(Err, "Config: expected error retrieving removed config - PASS");
   Err += ConfigHmixCheck.get(ConfigHmixDel2Check);
   CHECK_ERROR(Err, "Config: expected error retrieving removed config - PASS");

   // Finalize environments
   Pacer::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0; // successful completion

} // end of main
//===-----------------------------------------------------------------------===/
