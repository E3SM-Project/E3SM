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
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"

#include <iostream>
#include <string>

//------------------------------------------------------------------------------
// The test driver for Config. This creates a Configuration in YAML format,
// writes that configuration to a file and then reads it back in, verifying
// that values are the same.
//
int main(int argc, char *argv[]) {

   int Err = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);

   // Initialize the Machine Environment and retrieve the default environment
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   OMEGA::I4 MyTask       = DefEnv->getMyTask();
   bool IsMaster          = DefEnv->isMasterTask();

   // Define some variables for a reference configuration
   // These are meant to create a hierachy that tests at least two levels
   // and includes a variable of each supported type.
   // omega:
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
   OMEGA::I4 RefHmixI4        = 3;
   OMEGA::I8 RefHmixI8        = 123456789;
   OMEGA::R4 RefHmixR4        = 5.0;
   OMEGA::R8 RefHmixR8        = 1.234567890;
   std::string RefHmixStr     = "HmixString";
   bool RefHmixDel2On         = true;
   OMEGA::I4 RefHmixDel2I4    = 5;
   OMEGA::I8 RefHmixDel2I8    = 999999999;
   OMEGA::R4 RefHmixDel2R4    = 7.0;
   OMEGA::R8 RefHmixDel2R8    = 9.999999999;
   std::string RefHmixDel2Str = "HmixDel2String";
   bool RefVmixOn             = false;
   OMEGA::I4 RefVmixI4        = 4;
   OMEGA::I8 RefVmixI8        = 987654321;
   OMEGA::R4 RefVmixR4        = 6.0;
   OMEGA::R8 RefVmixR8        = 9.9876543210;
   std::string RefVmixStr     = "VmixString";
   int VecSize                = 5;
   std::vector<OMEGA::I4> RefVecI4{1, 2, 3, 4, 5};
   std::vector<OMEGA::I8> RefVecI8{123456789, 234567890, 345678901, 456789012,
                                   567890123};
   std::vector<OMEGA::R4> RefVecR4{1.2345, 2.3456, 3.4567, 4.5678, 5.6789};
   std::vector<OMEGA::R8> RefVecR8{1.23456789, 2.34567890, 3.45678901,
                                   4.56789012, 5.67890123};
   std::vector<bool> RefVecLog{true, false, true, false, true};
   std::vector<std::string> RefList{"first", "second", "third", "fourth",
                                    "fifth"};

   // Build up a reference configuration
   OMEGA::Config ConfigOmegaAll("omegaroot");
   OMEGA::Config ConfigOmegaRef("omega");
   OMEGA::Config ConfigHmixRef("Hmix");
   OMEGA::Config ConfigHmixDel2Ref("HmixDel2");
   OMEGA::Config ConfigVmixRef("Vmix");

   Err = ConfigHmixDel2Ref.add("HmixDel2On", RefHmixDel2On);
   if (Err != 0)
      LOG_INFO("Config test {}: adding bool FAIL", MyTask);
   Err = ConfigHmixDel2Ref.add("HmixDel2I4", RefHmixDel2I4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I4 FAIL", MyTask);
   Err = ConfigHmixDel2Ref.add("HmixDel2I8", RefHmixDel2I8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I8 FAIL", MyTask);
   Err = ConfigHmixDel2Ref.add("HmixDel2R4", RefHmixDel2R4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R4 FAIL", MyTask);
   Err = ConfigHmixDel2Ref.add("HmixDel2R8", RefHmixDel2R8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R8 FAIL", MyTask);
   Err = ConfigHmixDel2Ref.add("HmixDel2Str", RefHmixDel2Str);
   if (Err != 0)
      LOG_INFO("Config test {}: adding string FAIL", MyTask);

   Err = ConfigHmixRef.add("HmixOn", RefHmixOn);
   if (Err != 0)
      LOG_INFO("Config test {}: adding bool FAIL", MyTask);
   Err = ConfigHmixRef.add("HmixI4", RefHmixI4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I4 FAIL", MyTask);
   Err = ConfigHmixRef.add("HmixI8", RefHmixI8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I8 FAIL", MyTask);
   Err = ConfigHmixRef.add("HmixR4", RefHmixR4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R4 FAIL", MyTask);
   Err = ConfigHmixRef.add("HmixR8", RefHmixR8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R8 FAIL", MyTask);
   Err = ConfigHmixRef.add("HmixStr", RefHmixStr);
   if (Err != 0)
      LOG_INFO("Config test {}: adding string FAIL", MyTask);
   Err = ConfigHmixRef.add(ConfigHmixDel2Ref);
   if (Err != 0)
      LOG_INFO("Config test {}: adding sub-config FAIL", MyTask);

   Err = ConfigVmixRef.add("VmixOn", RefVmixOn);
   if (Err != 0)
      LOG_INFO("Config test {}: adding bool FAIL", MyTask);
   Err = ConfigVmixRef.add("VmixI4", RefVmixI4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I4 FAIL", MyTask);
   Err = ConfigVmixRef.add("VmixI8", RefVmixI8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I8 FAIL", MyTask);
   Err = ConfigVmixRef.add("VmixR4", RefVmixR4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R4 FAIL", MyTask);
   Err = ConfigVmixRef.add("VmixR8", RefVmixR8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R8 FAIL", MyTask);
   Err = ConfigVmixRef.add("VmixStr", RefVmixStr);
   if (Err != 0)
      LOG_INFO("Config test {}: adding string FAIL", MyTask);

   Err = ConfigOmegaRef.add(ConfigHmixRef);
   if (Err != 0)
      LOG_INFO("Config test {}: adding sub-config FAIL", MyTask);
   Err = ConfigOmegaRef.add(ConfigVmixRef);
   if (Err != 0)
      LOG_INFO("Config test {}: adding sub-config FAIL", MyTask);

   Err = ConfigOmegaRef.add("VectorI4", RefVecI4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I4 vector FAIL", MyTask);

   Err = ConfigOmegaRef.add("VectorI8", RefVecI8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding I8 vector FAIL", MyTask);

   Err = ConfigOmegaRef.add("VectorR4", RefVecR4);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R4 vector FAIL", MyTask);

   Err = ConfigOmegaRef.add("VectorR8", RefVecR8);
   if (Err != 0)
      LOG_INFO("Config test {}: adding R8 vector FAIL", MyTask);

   Err = ConfigOmegaRef.add("VectorLog", RefVecLog);
   if (Err != 0)
      LOG_INFO("Config test {}: adding bool vector FAIL", MyTask);

   Err = ConfigOmegaRef.add("StrList", RefList);
   if (Err != 0)
      LOG_INFO("Config test {}: adding string list FAIL", MyTask);

   // create the full root node
   Err = ConfigOmegaAll.add(ConfigOmegaRef);
   if (Err != 0)
      LOG_INFO("Config test {}: adding sub-config FAIL", MyTask);

   // Test retrievals by getting from the reference config and
   // checking against reference values.

   // Test logical retrievals for each config level
   // First extract sub-configurations from the top level config
   OMEGA::Config ConfigHmixNew("Hmix");
   OMEGA::Config ConfigVmixNew("Vmix");
   OMEGA::Config ConfigHmixDel2New("HmixDel2");
   int Err1   = ConfigOmegaRef.get(ConfigHmixNew);
   int Err2   = ConfigOmegaRef.get(ConfigVmixNew);
   int Err3   = ConfigHmixNew.get(ConfigHmixDel2New);
   int ErrAll = Err1 + Err2 + Err3;
   if (ErrAll == 0) {
      LOG_INFO("ConfigTest {}: Add/Get of subconfigs - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of subconfigs - FAIL", MyTask);
   }

   // Test logical retrievals
   bool HmixOn;
   bool VmixOn;
   bool HmixDel2On;
   Err1         = ConfigHmixNew.get("HmixOn", HmixOn);
   Err2         = ConfigVmixNew.get("VmixOn", VmixOn);
   Err3         = ConfigHmixDel2New.get("HmixDel2On", HmixDel2On);
   ErrAll       = Err1 + Err2 + Err3;
   bool RefTest = (HmixOn == RefHmixOn) && (VmixOn == RefVmixOn) &&
                  (HmixDel2On == RefHmixDel2On);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of boolean vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of boolean vars - FAIL", MyTask);
   }

   // Test I4 retrievals
   OMEGA::I4 HmixI4     = 0;
   OMEGA::I4 VmixI4     = 0;
   OMEGA::I4 HmixDel2I4 = 0;
   Err1                 = ConfigHmixNew.get("HmixI4", HmixI4);
   Err2                 = ConfigVmixNew.get("VmixI4", VmixI4);
   Err3                 = ConfigHmixDel2New.get("HmixDel2I4", HmixDel2I4);
   ErrAll               = Err1 + Err2 + Err3;
   RefTest              = (HmixI4 == RefHmixI4) && (VmixI4 == RefVmixI4) &&
             (HmixDel2I4 == RefHmixDel2I4);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of I4 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of I4 vars - FAIL", MyTask);
   }

   // Test I8 retrievals
   OMEGA::I8 HmixI8     = 0;
   OMEGA::I8 VmixI8     = 0;
   OMEGA::I8 HmixDel2I8 = 0;
   Err1                 = ConfigHmixNew.get("HmixI8", HmixI8);
   Err2                 = ConfigVmixNew.get("VmixI8", VmixI8);
   Err3                 = ConfigHmixDel2New.get("HmixDel2I8", HmixDel2I8);
   ErrAll               = Err1 + Err2 + Err3;
   RefTest              = (HmixI8 == RefHmixI8) && (VmixI8 == RefVmixI8) &&
             (HmixDel2I8 == RefHmixDel2I8);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of I8 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of I8 vars - FAIL", MyTask);
   }

   // Test R4 retrievals
   OMEGA::R4 HmixR4     = 0;
   OMEGA::R4 VmixR4     = 0;
   OMEGA::R4 HmixDel2R4 = 0;
   Err1                 = ConfigHmixNew.get("HmixR4", HmixR4);
   Err2                 = ConfigVmixNew.get("VmixR4", VmixR4);
   Err3                 = ConfigHmixDel2New.get("HmixDel2R4", HmixDel2R4);
   ErrAll               = Err1 + Err2 + Err3;
   RefTest              = (HmixR4 == RefHmixR4) && (VmixR4 == RefVmixR4) &&
             (HmixDel2R4 == RefHmixDel2R4);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of R4 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of R4 vars - FAIL", MyTask);
   }

   // Test R8 retrievals
   OMEGA::R8 HmixR8     = 0;
   OMEGA::R8 VmixR8     = 0;
   OMEGA::R8 HmixDel2R8 = 0;
   Err1                 = ConfigHmixNew.get("HmixR8", HmixR8);
   Err2                 = ConfigVmixNew.get("VmixR8", VmixR8);
   Err3                 = ConfigHmixDel2New.get("HmixDel2R8", HmixDel2R8);
   ErrAll               = Err1 + Err2 + Err3;
   RefTest              = (HmixR8 == RefHmixR8) && (VmixR8 == RefVmixR8) &&
             (HmixDel2R8 == RefHmixDel2R8);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of R8 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of R8 vars - FAIL", MyTask);
   }

   // Test string retrievals
   std::string HmixStr;
   std::string VmixStr;
   std::string HmixDel2Str;
   Err1    = ConfigHmixNew.get("HmixStr", HmixStr);
   Err2    = ConfigVmixNew.get("VmixStr", VmixStr);
   Err3    = ConfigHmixDel2New.get("HmixDel2Str", HmixDel2Str);
   ErrAll  = Err1 + Err2 + Err3;
   RefTest = ((HmixStr == RefHmixStr) && (VmixStr == RefVmixStr) &&
              (HmixDel2Str == RefHmixDel2Str));
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of string vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of string vars - FAIL", MyTask);
   }

   // Test vector retrievals
   std::vector<OMEGA::I4> VecI4;
   Err1    = ConfigOmegaRef.get("VectorI4", VecI4);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecI4[i] == RefVecI4[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of I4 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of I4 vector - FAIL", MyTask);
   }

   std::vector<OMEGA::I8> VecI8;
   Err1    = ConfigOmegaRef.get("VectorI8", VecI8);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecI8[i] == RefVecI8[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of I8 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of I8 vector - FAIL", MyTask);
   }

   std::vector<OMEGA::R4> VecR4;
   Err1    = ConfigOmegaRef.get("VectorR4", VecR4);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecR4[i] == RefVecR4[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of R4 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of R4 vector - FAIL", MyTask);
   }

   std::vector<OMEGA::R8> VecR8;
   Err1    = ConfigOmegaRef.get("VectorR8", VecR8);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecR8[i] == RefVecR8[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of R8 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of R8 vector - FAIL", MyTask);
   }

   std::vector<bool> VecLog;
   Err1    = ConfigOmegaRef.get("VectorLog", VecLog);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecLog[i] == RefVecLog[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of bool vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of bool vector - FAIL", MyTask);
   }

   std::vector<std::string> NewList;
   Err1    = ConfigOmegaRef.get("StrList", NewList);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (NewList[i] == RefList[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Add/Get of string list - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Add/Get of string list - FAIL", MyTask);
   }

   // Test changing values using the set functions
   bool NewHmixOn             = false;
   OMEGA::I4 NewHmixI4        = 4;
   OMEGA::I8 NewHmixI8        = 223456789;
   OMEGA::R4 NewHmixR4        = 6.0;
   OMEGA::R8 NewHmixR8        = 2.234567890;
   std::string NewHmixStr     = "HmixStringNew";
   bool NewHmixDel2On         = false;
   OMEGA::I4 NewHmixDel2I4    = 6;
   OMEGA::I8 NewHmixDel2I8    = 899999999;
   OMEGA::R4 NewHmixDel2R4    = 8.0;
   OMEGA::R8 NewHmixDel2R8    = 8.999999999;
   std::string NewHmixDel2Str = "HmixDel2StringNew";
   bool NewVmixOn             = true;
   OMEGA::I4 NewVmixI4        = 5;
   OMEGA::I8 NewVmixI8        = 887654321;
   OMEGA::R4 NewVmixR4        = 7.0;
   OMEGA::R8 NewVmixR8        = 8.9876543210;
   std::string NewVmixStr     = "VmixStringNew";

   Err = ConfigHmixDel2New.set("HmixDel2On", NewHmixDel2On);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set bool FAIL", MyTask);
   Err = ConfigHmixDel2New.set("HmixDel2I4", NewHmixDel2I4);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set I4 FAIL", MyTask);
   Err = ConfigHmixDel2New.set("HmixDel2I8", NewHmixDel2I8);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set I8 FAIL", MyTask);
   Err = ConfigHmixDel2New.set("HmixDel2R4", NewHmixDel2R4);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set R4 FAIL", MyTask);
   Err = ConfigHmixDel2New.set("HmixDel2R8", NewHmixDel2R8);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set R8 FAIL", MyTask);
   Err = ConfigHmixDel2New.set("HmixDel2Str", NewHmixDel2Str);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set string FAIL", MyTask);

   Err = ConfigHmixNew.set("HmixOn", NewHmixOn);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set bool FAIL", MyTask);
   Err = ConfigHmixNew.set("HmixI4", NewHmixI4);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set I4 FAIL", MyTask);
   Err = ConfigHmixNew.set("HmixI8", NewHmixI8);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set I8 FAIL", MyTask);
   Err = ConfigHmixNew.set("HmixR4", NewHmixR4);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set R4 FAIL", MyTask);
   Err = ConfigHmixNew.set("HmixR8", NewHmixR8);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set R8 FAIL", MyTask);
   Err = ConfigHmixNew.set("HmixStr", NewHmixStr);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set string FAIL", MyTask);

   Err = ConfigVmixNew.set("VmixOn", NewVmixOn);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set bool FAIL", MyTask);
   Err = ConfigVmixNew.set("VmixI4", NewVmixI4);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set I4 FAIL", MyTask);
   Err = ConfigVmixNew.set("VmixI8", NewVmixI8);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set I8 FAIL", MyTask);
   Err = ConfigVmixNew.set("VmixR4", NewVmixR4);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set R4 FAIL", MyTask);
   Err = ConfigVmixNew.set("VmixR8", NewVmixR8);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set R8 FAIL", MyTask);
   Err = ConfigVmixNew.set("VmixStr", NewVmixStr);
   if (Err != 0)
      LOG_INFO("ConfigTest {}: set string FAIL", MyTask);

   // Test logical retrievals
   Err1    = ConfigHmixNew.get("HmixOn", HmixOn);
   Err2    = ConfigVmixNew.get("VmixOn", VmixOn);
   Err3    = ConfigHmixDel2New.get("HmixDel2On", HmixDel2On);
   ErrAll  = Err1 + Err2 + Err3;
   RefTest = (HmixOn == NewHmixOn) && (VmixOn == NewVmixOn) &&
             (HmixDel2On == NewHmixDel2On);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of boolean vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of boolean vars - FAIL", MyTask);
   }

   // Test I4 retrievals
   Err1    = ConfigHmixNew.get("HmixI4", HmixI4);
   Err2    = ConfigVmixNew.get("VmixI4", VmixI4);
   Err3    = ConfigHmixDel2New.get("HmixDel2I4", HmixDel2I4);
   ErrAll  = Err1 + Err2 + Err3;
   RefTest = (HmixI4 == NewHmixI4) && (VmixI4 == NewVmixI4) &&
             (HmixDel2I4 == NewHmixDel2I4);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of I4 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of I4 vars - FAIL", MyTask);
   }

   // Test I8 retrievals
   Err1    = ConfigHmixNew.get("HmixI8", HmixI8);
   Err2    = ConfigVmixNew.get("VmixI8", VmixI8);
   Err3    = ConfigHmixDel2New.get("HmixDel2I8", HmixDel2I8);
   ErrAll  = Err1 + Err2 + Err3;
   RefTest = (HmixI8 == NewHmixI8) && (VmixI8 == NewVmixI8) &&
             (HmixDel2I8 == NewHmixDel2I8);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of I8 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of I8 vars - FAIL", MyTask);
   }

   // Test R4 retrievals
   Err1    = ConfigHmixNew.get("HmixR4", HmixR4);
   Err2    = ConfigVmixNew.get("VmixR4", VmixR4);
   Err3    = ConfigHmixDel2New.get("HmixDel2R4", HmixDel2R4);
   ErrAll  = Err1 + Err2 + Err3;
   RefTest = (HmixR4 == NewHmixR4) && (VmixR4 == NewVmixR4) &&
             (HmixDel2R4 == NewHmixDel2R4);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of R4 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of R4 vars - FAIL", MyTask);
   }

   // Test R8 retrievals
   Err1    = ConfigHmixNew.get("HmixR8", HmixR8);
   Err2    = ConfigVmixNew.get("VmixR8", VmixR8);
   Err3    = ConfigHmixDel2New.get("HmixDel2R8", HmixDel2R8);
   ErrAll  = Err1 + Err2 + Err3;
   RefTest = (HmixR8 == NewHmixR8) && (VmixR8 == NewVmixR8) &&
             (HmixDel2R8 == NewHmixDel2R8);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of R8 vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of R8 vars - FAIL", MyTask);
   }

   // Test string retrievals
   Err1    = ConfigHmixNew.get("HmixStr", HmixStr);
   Err2    = ConfigVmixNew.get("VmixStr", VmixStr);
   Err3    = ConfigHmixDel2New.get("HmixDel2Str", HmixDel2Str);
   ErrAll  = Err1 + Err2 + Err3;
   RefTest = (HmixStr == NewHmixStr) && (VmixStr == NewVmixStr) &&
             (HmixDel2Str == NewHmixDel2Str);
   if (ErrAll == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of string vars - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of string vars - FAIL", MyTask);
   }

   // Test setting of vectors by modifying one entry, setting and retrieving
   RefVecI4[2]  = -1;
   RefVecI8[2]  = -123456789;
   RefVecR4[2]  = -1.2345;
   RefVecR8[2]  = -1.23456789;
   RefVecLog[2] = false;
   RefList[2]   = "junk";

   Err = ConfigOmegaRef.set("VectorI4", RefVecI4);
   if (Err != 0)
      LOG_INFO("Config test {}: setting I4 vector FAIL", MyTask);

   Err = ConfigOmegaRef.set("VectorI8", RefVecI8);
   if (Err != 0)
      LOG_INFO("Config test {}: setting I8 vector FAIL", MyTask);

   Err = ConfigOmegaRef.set("VectorR4", RefVecR4);
   if (Err != 0)
      LOG_INFO("Config test {}: setting R4 vector FAIL", MyTask);

   Err = ConfigOmegaRef.set("VectorR8", RefVecR8);
   if (Err != 0)
      LOG_INFO("Config test {}: setting R8 vector FAIL", MyTask);

   Err = ConfigOmegaRef.set("VectorLog", RefVecLog);
   if (Err != 0)
      LOG_INFO("Config test {}: setting bool vector FAIL", MyTask);

   Err = ConfigOmegaRef.set("StrList", RefList);
   if (Err != 0)
      LOG_INFO("Config test {}: setting string list FAIL", MyTask);

   // Test vector retrievals after set
   Err1    = ConfigOmegaRef.get("VectorI4", VecI4);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecI4[i] == RefVecI4[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of I4 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of I4 vector - FAIL", MyTask);
   }

   Err1    = ConfigOmegaRef.get("VectorI8", VecI8);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecI8[i] == RefVecI8[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of I8 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of I8 vector - FAIL", MyTask);
   }

   Err1    = ConfigOmegaRef.get("VectorR4", VecR4);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecR4[i] == RefVecR4[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of R4 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of R4 vector - FAIL", MyTask);
   }

   Err1    = ConfigOmegaRef.get("VectorR8", VecR8);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecR8[i] == RefVecR8[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of R8 vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of R8 vector - FAIL", MyTask);
   }

   Err1    = ConfigOmegaRef.get("VectorLog", VecLog);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecLog[i] == RefVecLog[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of bool vector - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of bool vector - FAIL", MyTask);
   }

   Err1    = ConfigOmegaRef.get("StrList", NewList);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (NewList[i] == RefList[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Set/Get of string list - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Set/Get of string list - FAIL", MyTask);
   }

   // Reset vectors back to original reference values
   RefVecI4[2]  = 3;
   RefVecI8[2]  = 345678901;
   RefVecR4[2]  = 3.4567;
   RefVecR8[2]  = 3.45678901;
   RefVecLog[2] = true;
   RefList[2]   = "third";
   Err          = ConfigOmegaRef.set("VectorI4", RefVecI4);
   Err          = ConfigOmegaRef.set("VectorI8", RefVecI8);
   Err          = ConfigOmegaRef.set("VectorR4", RefVecR4);
   Err          = ConfigOmegaRef.set("VectorR8", RefVecR8);
   Err          = ConfigOmegaRef.set("VectorLog", RefVecLog);
   Err          = ConfigOmegaRef.set("StrList", RefList);

   // Test add,set,get error modes
   OMEGA::Config ConfigJunk("junk");
   bool JunkOn;
   OMEGA::I4 JunkI4;
   OMEGA::I4 JunkI8;
   OMEGA::I4 JunkR4;
   OMEGA::I4 JunkR8;
   std::string JunkStr;

   // Try to add a variable or config that already exists
   // Try to set/get a variable or config that does not exist
   Err1    = ConfigOmegaRef.add(ConfigHmixNew);
   Err2    = ConfigOmegaRef.get(ConfigJunk);
   RefTest = (Err1 != 0 && Err2 != 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: SubConfig add/get error modes PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: SubConfig add/get error modes FAIL", MyTask);
   }

   Err1    = ConfigHmixNew.add("HmixOn", NewHmixOn);
   Err2    = ConfigHmixNew.set("junk", JunkOn);
   Err3    = ConfigHmixNew.get("junk", JunkOn);
   RefTest = (Err1 != 0 && Err2 != 0 && Err3 != 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: add/set/get bool error modes PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: add/set/get bool error modes FAIL", MyTask);
   }

   Err1    = ConfigHmixNew.add("HmixI4", NewHmixI4);
   Err2    = ConfigHmixNew.set("junk", JunkI4);
   Err3    = ConfigHmixNew.get("junk", JunkI4);
   RefTest = (Err1 != 0 && Err2 != 0 && Err3 != 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: add/set/get I4 error modes PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: add/set/get I4 error modes FAIL", MyTask);
   }

   Err1    = ConfigHmixNew.add("HmixI8", NewHmixI8);
   Err2    = ConfigHmixNew.set("junk", JunkI8);
   Err3    = ConfigHmixNew.get("junk", JunkI8);
   RefTest = (Err1 != 0 && Err2 != 0 && Err3 != 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: add/set/get I8 error modes PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: add/set/get I8 error modes FAIL", MyTask);
   }

   Err1    = ConfigHmixNew.add("HmixR4", NewHmixR4);
   Err2    = ConfigHmixNew.set("junk", JunkR4);
   Err3    = ConfigHmixNew.get("junk", JunkR4);
   RefTest = (Err1 != 0 && Err2 != 0 && Err3 != 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: add/set/get R4 error modes PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: add/set/get R4 error modes FAIL", MyTask);
   }

   Err1    = ConfigHmixNew.add("HmixR8", NewHmixR8);
   Err2    = ConfigHmixNew.set("junk", JunkR8);
   Err3    = ConfigHmixNew.get("junk", JunkR8);
   RefTest = (Err1 != 0 && Err2 != 0 && Err3 != 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: add/set/get R8 error modes PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: add/set/get R8 error modes FAIL", MyTask);
   }

   Err1    = ConfigHmixNew.add("HmixStr", NewHmixStr);
   Err2    = ConfigHmixNew.set("junk", JunkStr);
   Err3    = ConfigHmixNew.get("junk", JunkStr);
   RefTest = (Err1 != 0 && Err2 != 0 && Err3 != 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: add/set/get string error modes PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: add/set/get string error modes FAIL", MyTask);
   }

   // Write the reference config to a file
   Err = ConfigOmegaAll.write("omegaConfigTst.yml");
   if (Err == 0) {
      LOG_INFO("ConfigTest {}: write successful PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: write successful FAIL", MyTask);
   }

   // Read the environment back in
   Err = OMEGA::Config::readAll("omegaConfigTst.yml");
   if (Err == 0) {
      LOG_INFO("ConfigTest {}: read successful PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: read successful FAIL", MyTask);
   }

   // Get the full configuration
   OMEGA::Config *ConfigOmega = OMEGA::Config::getOmegaConfig();

   // Retrieve all values and compare against reference values
   OMEGA::Config ConfigHmixOmega("Hmix");
   OMEGA::Config ConfigVmixOmega("Vmix");
   OMEGA::Config ConfigHmixDel2Omega("HmixDel2");

   Err1 = ConfigOmega->get(ConfigHmixOmega);
   Err2 = ConfigOmega->get(ConfigVmixOmega);
   Err3 = ConfigHmixOmega.get(ConfigHmixDel2Omega);

   RefTest = (Err1 == 0 && Err2 == 0 && Err3 == 0);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: retrieve subgroups from full config PASS",
               MyTask);
   } else {
      LOG_INFO("ConfigTest {}: retrieve subgroups from full config FAIL",
               MyTask);
   }

   Err1    = ConfigHmixOmega.get("HmixOn", HmixOn);
   Err2    = ConfigVmixOmega.get("VmixOn", VmixOn);
   Err3    = ConfigHmixDel2Omega.get("HmixDel2On", HmixDel2On);
   RefTest = (Err1 == 0 && Err2 == 0 && Err3 == 0 && HmixOn == NewHmixOn &&
              VmixOn == NewVmixOn && HmixDel2On == NewHmixDel2On);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: retrieve bool from full config PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: retrieve bool from full config FAIL", MyTask);
   }

   Err1    = ConfigHmixOmega.get("HmixI4", HmixI4);
   Err2    = ConfigVmixOmega.get("VmixI4", VmixI4);
   Err3    = ConfigHmixDel2Omega.get("HmixDel2I4", HmixDel2I4);
   RefTest = (Err1 == 0 && Err2 == 0 && Err3 == 0 && HmixI4 == NewHmixI4 &&
              VmixI4 == NewVmixI4 && HmixDel2I4 == NewHmixDel2I4);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: retrieve I4 from full config PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: retrieve I4 from full config FAIL", MyTask);
   }

   Err1    = ConfigHmixOmega.get("HmixI8", HmixI8);
   Err2    = ConfigVmixOmega.get("VmixI8", VmixI8);
   Err3    = ConfigHmixDel2Omega.get("HmixDel2I8", HmixDel2I8);
   RefTest = (Err1 == 0 && Err2 == 0 && Err3 == 0 && HmixI8 == NewHmixI8 &&
              VmixI8 == NewVmixI8 && HmixDel2I8 == NewHmixDel2I8);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: retrieve I8 from full config PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: retrieve I8 from full config FAIL", MyTask);
   }

   Err1    = ConfigHmixOmega.get("HmixR4", HmixR4);
   Err2    = ConfigVmixOmega.get("VmixR4", VmixR4);
   Err3    = ConfigHmixDel2Omega.get("HmixDel2R4", HmixDel2R4);
   RefTest = (Err1 == 0 && Err2 == 0 && Err3 == 0 && HmixR4 == NewHmixR4 &&
              VmixR4 == NewVmixR4 && HmixDel2R4 == NewHmixDel2R4);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: retrieve R4 from full config PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: retrieve R4 from full config FAIL", MyTask);
   }

   Err1    = ConfigHmixOmega.get("HmixR8", HmixR8);
   Err2    = ConfigVmixOmega.get("VmixR8", VmixR8);
   Err3    = ConfigHmixDel2Omega.get("HmixDel2R8", HmixDel2R8);
   RefTest = (Err1 == 0 && Err2 == 0 && Err3 == 0 && HmixR8 == NewHmixR8 &&
              VmixR8 == NewVmixR8 && HmixDel2R8 == NewHmixDel2R8);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: retrieve R8 from full config PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: retrieve R8 from full config FAIL", MyTask);
   }

   Err1    = ConfigHmixOmega.get("HmixStr", HmixStr);
   Err2    = ConfigVmixOmega.get("VmixStr", VmixStr);
   Err3    = ConfigHmixDel2Omega.get("HmixDel2Str", HmixDel2Str);
   RefTest = (Err1 == 0 && Err2 == 0 && Err3 == 0 && HmixStr == NewHmixStr &&
              VmixStr == NewVmixStr && HmixDel2Str == NewHmixDel2Str);
   if (RefTest) {
      LOG_INFO("ConfigTest {}: retrieve string from full config PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: retrieve string from full config FAIL", MyTask);
   }

   // Vector retrievals
   Err1    = ConfigOmega->get("VectorI4", VecI4);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecI4[i] == RefVecI4[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Get I4 vector from full config - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Get I4 vector from full config - FAIL", MyTask);
   }

   Err1    = ConfigOmega->get("VectorI8", VecI8);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecI8[i] == RefVecI8[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Get I8 vector from full config - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Get I8 vector from full config - FAIL", MyTask);
   }

   Err1    = ConfigOmega->get("VectorR4", VecR4);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecR4[i] == RefVecR4[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Get R4 vector from full config - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Get R4 vector from full config - FAIL", MyTask);
   }

   Err1    = ConfigOmega->get("VectorR8", VecR8);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecR8[i] == RefVecR8[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Get R8 vector from full config - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Get R8 vector from full config - FAIL", MyTask);
   }

   Err1    = ConfigOmega->get("VectorLog", VecLog);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (VecLog[i] == RefVecLog[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Get bool vec from full config - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Get bool vec from full config - FAIL", MyTask);
   }

   Err1    = ConfigOmega->get("StrList", NewList);
   RefTest = true;
   for (int i = 0; i < VecSize; ++i) {
      RefTest = RefTest && (NewList[i] == RefList[i]);
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Get string list full config - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Get string list full config - FAIL", MyTask);
   }

   // Test retrieval of values using an iterator
   Err1      = 0;
   RefTest   = true;
   int IList = 0;
   for (auto It = ConfigOmega->begin(); It != ConfigOmega->end(); ++It) {
      std::string NodeName;
      Err1 = OMEGA::Config::getName(It, NodeName);
      // Note that iteration order not guaranteed but so far has
      // been consistent with location in the file.
      if (IList == 0 and NodeName != "Hmix")
         RefTest = false;
      if (IList == 1 and NodeName != "Vmix")
         RefTest = false;
      if (IList == 2 and NodeName != "VectorI4")
         RefTest = false;
      if (IList == 3 and NodeName != "VectorI8")
         RefTest = false;
      if (IList == 4 and NodeName != "VectorR4")
         RefTest = false;
      if (IList == 5 and NodeName != "VectorR8")
         RefTest = false;
      if (IList == 6 and NodeName != "VectorLog")
         RefTest = false;
      if (IList == 7 and NodeName != "StrList")
         RefTest = false;
      if (IList > 7)
         RefTest = false;
      ++IList;
   }
   if (Err1 == 0 && RefTest) {
      LOG_INFO("ConfigTest {}: Get config list using iter - PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: Get config list using iter - FAIL", MyTask);
   }

   // Test removals by removing both variables and sub-configs and
   // checking the further retrievals fail

   int Err4;
   int Err5;
   int Err6;
   Err1 = ConfigVmixOmega.remove("VmixOn");
   Err2 = ConfigVmixOmega.remove("VmixI4");
   Err3 = ConfigVmixOmega.remove("VmixI8");
   Err4 = ConfigVmixOmega.remove("VmixR4");
   Err5 = ConfigVmixOmega.remove("VmixR8");
   Err6 = ConfigVmixOmega.remove("VmixStr");

   if (Err1 == 0 && Err2 == 0 && Err3 == 0 && Err4 == 0 && Err5 == 0 &&
       Err6 == 0) {
      LOG_INFO("ConfigTest {}: removal call PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: removal call FAIL", MyTask);
   }

   // Retrieve a fresh copy to make sure the removal happens to the
   // entire Config, not just local copy. All variable retrievals should fail.

   OMEGA::Config ConfigVmixCheck("Vmix");
   Err  = ConfigOmega->get(ConfigVmixCheck);
   Err1 = ConfigVmixCheck.get("VmixOn", NewVmixOn);
   Err2 = ConfigVmixCheck.get("VmixI4", NewVmixI4);
   Err3 = ConfigVmixCheck.get("VmixI8", NewVmixI8);
   Err4 = ConfigVmixCheck.get("VmixR4", NewVmixR4);
   Err5 = ConfigVmixCheck.get("VmixR8", NewVmixR8);
   Err6 = ConfigVmixCheck.get("VmixStr", NewVmixStr);

   if (Err1 != 0 && Err2 != 0 && Err3 != 0 && Err4 != 0 && Err5 != 0 &&
       Err6 != 0 && Err == 0) {
      LOG_INFO("ConfigTest {}: variable removal test PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: variable removal test FAIL", MyTask);
   }

   // Now try removing an entire sub-config
   Err1 = ConfigHmixOmega.remove("HmixDel2"); // remove entire subconfig
   OMEGA::Config ConfigHmixCheck("Hmix");
   OMEGA::Config ConfigHmixDel2Check("HmixDel2");
   Err2 = ConfigOmega->get(ConfigHmixCheck);
   Err3 = ConfigHmixCheck.get(ConfigHmixDel2Check);
   if (Err1 == 0 && Err2 == 0 && Err3 != 0) {
      LOG_INFO("ConfigTest {}: subconfig removal test PASS", MyTask);
   } else {
      LOG_INFO("ConfigTest {}: subconfig removal test FAIL", MyTask);
   }

   // Finalize environments
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/
