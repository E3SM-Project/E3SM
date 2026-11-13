//===-- Test driver for Omega CF metadata ---------------------*- C++ -*-===//
//
/// \file
/// \brief Test driver for CF compliance of Omega field and file metadata
///
/// Initializes every Omega module so that all fields are defined, along with
/// the GlobalStats and MonthlyAverages analysis groups so that the fields
/// their operators derive are defined too, then checks the metadata that ends
/// up as netCDF attributes against the rules the CF conventions impose and
/// that a CF checker enforces:
///   1. Every units string is written in the plain udunits form used by the CF
///      standard name table (m s-1, m2, N m-2, degree_C, 1) and none of the
///      invalid spellings (carets, slashes, braces, "dimensionless", "none")
///   2. No field carries an empty units or standard_name entry and every
///      standard_name has the form of a CF standard name
///   3. valid_min does not exceed valid_max and is not a tiny positive number
///      (the smallest positive Real, which excludes every negative value)
///   4. No field carries the retired non-CF duplicates of the standard
///      attributes (Name, name, Description, Units, StdName, ValidMin,
///      ValidMax)
///   5. The global code metadata declares the Conventions attribute
///   6. A file written through IOStream carries the Conventions attribute, no
///      empty or malformed units or standard_name attributes and none of the
///      retired duplicates
///
/// The standard names themselves can only be validated against the CF
/// standard name table, which is left to an external CF checker.
//
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "Forcing.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "KPPMix.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "Pacer.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "VertMix.h"
#include "mpi.h"
#include "yaml-cpp/yaml.h"

#include <algorithm>
#include <fstream>
#include <netcdf.h>
#include <regex>
#include <string>
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Name and file of the stream written for the file-level checks
static const std::string CFStreamName = "CFMetadataCheck";
static const std::string CFFileName   = "ocn.cfmeta.nc";

//------------------------------------------------------------------------------
// Returns true if a units string is written in the CF plain form: either the
// dimensionless "1" or factors separated by single spaces, each a unit name
// (letters and underscores, eg m, kg, degree_C, radians) optionally followed
// by a signed integer exponent with no caret (m2, s-1), optionally preceded by
// the positive scale factor of a scaled unit (1e+06 m3 s-1). Time coordinate
// units of the form "<unit> since <reference time>" are also accepted.
bool isCFPlainFormUnits(const std::string &Units) {

   static const std::regex PlainForm(
       "1|([0-9]*\\.?[0-9]+(e[+-]?[0-9]+)? )?"
       "[A-Za-z_]+(-?[1-9][0-9]*)?( [A-Za-z_]+(-?[1-9][0-9]*)?)*");
   static const std::regex TimeForm("[A-Za-z_]+ since .+");

   // udunits accepts these words but CF spells a dimensionless quantity "1"
   // and "none" is not a unit at all
   if (Units == "dimensionless" or Units == "none")
      return false;

   return std::regex_match(Units, PlainForm) or
          std::regex_match(Units, TimeForm);
}

//------------------------------------------------------------------------------
// Returns true if a standard name has the form used by every entry of the CF
// standard name table (lower-case letters, digits and underscores)
bool isCFStandardNameForm(const std::string &StdName) {
   static const std::regex NameForm("[a-z][a-z0-9_]*");
   return std::regex_match(StdName, NameForm);
}

//------------------------------------------------------------------------------
// Returns true if an attribute name is one of the non-CF duplicates of the
// standard attributes that Omega used to write alongside the CF ones
bool isRetiredAttribute(const std::string &AttName) {
   static const std::vector<std::string> Retired = {
       "Name",    "name",     "Description", "Units",
       "StdName", "ValidMin", "ValidMax"};
   return std::find(Retired.begin(), Retired.end(), AttName) != Retired.end();
}

//------------------------------------------------------------------------------
// Extracts a metadata value of any supported numeric type as a double.
// Returns false if the value is not numeric.
bool metadataAsDouble(const std::any &Value, double &Result) {
   if (Value.type() == typeid(R8)) {
      Result = std::any_cast<R8>(Value);
   } else if (Value.type() == typeid(R4)) {
      Result = std::any_cast<R4>(Value);
   } else if (Value.type() == typeid(I4)) {
      Result = std::any_cast<I4>(Value);
   } else if (Value.type() == typeid(I8)) {
      Result = std::any_cast<I8>(Value);
   } else {
      return false;
   }
   return true;
}

//------------------------------------------------------------------------------
// Extracts a metadata string value, whether stored as std::string or as a
// string literal (const char *). Returns false if the value is not a string.
bool metadataAsString(const std::any &Value, std::string &Result) {
   if (Value.type() == typeid(std::string)) {
      Result = std::any_cast<std::string>(Value);
   } else if (Value.type() == typeid(const char *)) {
      Result = std::any_cast<const char *>(Value);
   } else {
      return false;
   }
   return true;
}

//------------------------------------------------------------------------------
// Writes a copy of the input config that adds an OnDemand write stream with a
// fixed file name and a broad set of field groups and enables the analysis
// groups that do not need region masks, then returns the new file name. Only
// task 0 writes the file; all tasks then read the same copy.
std::string addCFStreamConfig(const std::string &InFile, MachEnv *Env) {

   std::string OutFile = "omega.cfmeta.yml";

   if (Env->getMyTask() == 0) {
      YAML::Node Root       = YAML::LoadFile(InFile);
      YAML::Node Streams    = Root["Omega"]["IOStreams"];
      YAML::Node CFStream   = YAML::Clone(Streams["History"]);
      CFStream["Filename"]  = CFFileName;
      CFStream["FreqUnits"] = "OnDemand";
      YAML::Node Contents(YAML::NodeType::Sequence);
      for (const std::string Group :
           {"Tracers", "State", "SshCell", "AuxiliaryState", "VertAdv",
            "VertMix", "KPPMix", "Eos", "Tendencies", "Forcing"})
         Contents.push_back(Group);
      CFStream["Contents"]  = Contents;
      Streams[CFStreamName] = CFStream;

      YAML::Node Analysis                   = Root["Omega"]["Analysis"];
      Analysis["GlobalStats"]["Enable"]     = true;
      Analysis["MonthlyAverages"]["Enable"] = true;

      std::ofstream OutStream(OutFile);
      OutStream << Root;
   }
   MPI_Barrier(Env->getComm());

   return OutFile;

} // End addCFStreamConfig

//------------------------------------------------------------------------------
// Full initialization of all Omega modules so that every field is defined,
// following the order used by the ocean driver (initOmegaModules)
Clock *initCFMetadataTest() {

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ CF Metadata Tests ------");

   Config("Omega");
   std::string ConfigFile = addCFStreamConfig("omega.yml", DefEnv);
   Config::readAll(ConfigFile);

   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   IOStream::init(ModelClock);
   Field::init(ModelClock);
   Decomp::init();

   int Err = Halo::init();
   if (Err != 0)
      ABORT_ERROR("CFMetadataTest: error initializing default halo");

   HorzMesh::init(ModelClock);
   VertCoord::init();
   Tracers::init();
   VertAdv::init();
   Forcing::init();
   AuxiliaryState::init();
   Eos::init();
   PressureGrad::init();
   VertMix::init();
   KPPMix::init();
   Tendencies::init();
   TimeStepper::init2();

   Err = OceanState::init();
   if (Err != 0)
      ABORT_ERROR("CFMetadataTest: error initializing OceanState");

   Analysis::init();

   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid)
      ABORT_ERROR("CFMetadataTest: error validating IO streams");

   return ModelClock;
}

//------------------------------------------------------------------------------
// Checks the in-memory metadata of every defined field. Returns the number
// of fields with problems; each problem is logged.
int checkFieldMetadata() {

   int NBad = 0;

   for (const std::string &FieldName : Field::getAllFieldNames()) {

      std::shared_ptr<Metadata> AllMeta = Field::getFieldMetadata(FieldName);
      bool Bad                          = false;

      // Units: must be absent or in the CF plain form
      auto UnitsIt = AllMeta->find("units");
      if (UnitsIt != AllMeta->end()) {
         std::string Units;
         if (!metadataAsString(UnitsIt->second, Units) or Units.empty() or
             !isCFPlainFormUnits(Units)) {
            LOG_ERROR("CFMetadataTest: field {} has units \"{}\" that are not "
                      "in the CF plain form (eg m s-1, m2, N m-2, 1)",
                      FieldName, Units);
            Bad = true;
         }
      }

      // Standard name: must be absent or shaped like a CF standard name
      auto StdNameIt = AllMeta->find("standard_name");
      if (StdNameIt != AllMeta->end()) {
         std::string StdName;
         if (!metadataAsString(StdNameIt->second, StdName) or
             !isCFStandardNameForm(StdName)) {
            LOG_ERROR("CFMetadataTest: field {} has standard_name \"{}\" "
                      "that is empty or not of CF form",
                      FieldName, StdName);
            Bad = true;
         }
      }

      // Retired duplicates of the standard attributes must not be present
      for (const auto &MetaPair : *AllMeta) {
         if (isRetiredAttribute(MetaPair.first)) {
            LOG_ERROR("CFMetadataTest: field {} has retired attribute {}",
                      FieldName, MetaPair.first);
            Bad = true;
         }
      }

      // Valid range: min must not exceed max, and a tiny positive min is the
      // signature of numeric_limits::min() used where lowest() was meant
      auto MinIt = AllMeta->find("valid_min");
      auto MaxIt = AllMeta->find("valid_max");
      if (MinIt != AllMeta->end() and MaxIt != AllMeta->end()) {
         double ValidMin, ValidMax;
         if (metadataAsDouble(MinIt->second, ValidMin) and
             metadataAsDouble(MaxIt->second, ValidMax)) {
            if (ValidMin > ValidMax) {
               LOG_ERROR("CFMetadataTest: field {} has valid_min {} greater "
                         "than valid_max {}",
                         FieldName, ValidMin, ValidMax);
               Bad = true;
            }
            if (ValidMin > 0.0 and ValidMin < 1.0e-30) {
               LOG_ERROR("CFMetadataTest: field {} has valid_min {} - use "
                         "numeric_limits::lowest() rather than min() for a "
                         "signed field",
                         FieldName, ValidMin);
               Bad = true;
            }
         }
      }

      if (Bad)
         ++NBad;
   }

   return NBad;
}

//------------------------------------------------------------------------------
// Reads a string attribute of a netCDF variable (or NC_GLOBAL). Returns false
// if the attribute does not exist.
bool readStringAttr(int NcID, int VarID, const char *AttName,
                    std::string &Value) {
   size_t Len = 0;
   if (nc_inq_attlen(NcID, VarID, AttName, &Len) != NC_NOERR)
      return false;
   std::vector<char> Buf(Len + 1, '\0');
   if (nc_get_att_text(NcID, VarID, AttName, Buf.data()) != NC_NOERR)
      return false;
   Value = std::string(Buf.data(), Len);
   return true;
}

//------------------------------------------------------------------------------
// Checks the attributes of the file written by the CF stream, on task 0 only.
// Returns the number of problems found; each problem is logged.
int checkFileMetadata(const std::string &FileName) {

   int NBad  = 0;
   int NcID  = 0;
   int NcErr = nc_open(FileName.c_str(), NC_NOWRITE, &NcID);
   if (NcErr != NC_NOERR) {
      LOG_ERROR("CFMetadataTest: cannot open {}: {}", FileName,
                nc_strerror(NcErr));
      return 1;
   }

   // Global Conventions attribute
   std::string Conventions;
   if (!readStringAttr(NcID, NC_GLOBAL, "Conventions", Conventions)) {
      LOG_ERROR("CFMetadataTest: {} has no Conventions attribute", FileName);
      ++NBad;
   } else if (Conventions != CFConventions) {
      LOG_ERROR("CFMetadataTest: {} has Conventions \"{}\", expected \"{}\"",
                FileName, Conventions, CFConventions);
      ++NBad;
   }

   // Per-variable units and standard_name attributes
   int NVars = 0;
   nc_inq_nvars(NcID, &NVars);
   for (int VarID = 0; VarID < NVars; ++VarID) {
      char VarNameBuf[NC_MAX_NAME + 1];
      nc_inq_varname(NcID, VarID, VarNameBuf);
      std::string VarName(VarNameBuf);

      std::string Units;
      if (readStringAttr(NcID, VarID, "units", Units) and
          (Units.empty() or !isCFPlainFormUnits(Units))) {
         LOG_ERROR("CFMetadataTest: {} variable {} has units \"{}\" that are "
                   "empty or not in the CF plain form",
                   FileName, VarName, Units);
         ++NBad;
      }

      std::string StdName;
      if (readStringAttr(NcID, VarID, "standard_name", StdName) and
          !isCFStandardNameForm(StdName)) {
         LOG_ERROR("CFMetadataTest: {} variable {} has standard_name \"{}\" "
                   "that is empty or not of CF form",
                   FileName, VarName, StdName);
         ++NBad;
      }

      int NAtts = 0;
      nc_inq_varnatts(NcID, VarID, &NAtts);
      for (int AttID = 0; AttID < NAtts; ++AttID) {
         char AttNameBuf[NC_MAX_NAME + 1];
         nc_inq_attname(NcID, VarID, AttID, AttNameBuf);
         if (isRetiredAttribute(AttNameBuf)) {
            LOG_ERROR("CFMetadataTest: {} variable {} has retired attribute {}",
                      FileName, VarName, AttNameBuf);
            ++NBad;
         }
      }
   }

   nc_close(NcID);
   return NBad;
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[]) {

   Error ErrAll;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      Clock *ModelClock = initCFMetadataTest();
      MachEnv *DefEnv   = MachEnv::getDefault();

      // ------------------------------------------------------------------
      // Test 1: the units-form checker itself accepts the CF plain form and
      //         rejects the spellings the CF checker flags
      // ------------------------------------------------------------------
      {
         int Err = 0;
         for (const char *Good :
              {"m", "m s-1", "m2 s-2", "m-2 s-1", "N m-2", "kg m-2 s-1",
               "m3 kg-1 K-1", "g kg-1", "degree_C", "radians s-1", "Pa", "1",
               "1e+06 m3 s-1", "0.001 kg", "seconds since 0001-01-01 00:00:00"})
            if (!isCFPlainFormUnits(Good))
               ++Err;
         for (const char *Bad :
              {"", "m/s", "m s^-1", "N m^{-2}", "m^2 s^-2", "dimensionless",
               "none", "kg/m^3/s", "m  s-1", " m", "m ", "1e+06", "-2 m"})
            if (isCFPlainFormUnits(Bad))
               ++Err;

         if (Err == 0) {
            LOG_INFO("CFMetadataTest: units form checker PASS");
         } else {
            ErrAll += Error(ErrorCode::Fail,
                            "CFMetadataTest: units form checker FAIL");
         }
      }

      // ------------------------------------------------------------------
      // Test 2: metadata of every defined field
      // ------------------------------------------------------------------
      {
         int NBad = checkFieldMetadata();
         if (NBad == 0) {
            LOG_INFO("CFMetadataTest: field metadata PASS");
         } else {
            ErrAll +=
                Error(ErrorCode::Fail,
                      "CFMetadataTest: field metadata FAIL ({} fields)", NBad);
         }
      }

      // ------------------------------------------------------------------
      // Test 3: global code metadata declares the conventions
      // ------------------------------------------------------------------
      {
         std::shared_ptr<Field> CodeField = Field::get(CodeMeta);
         std::string Conventions;
         Error Err = CodeField->getMetadata("Conventions", Conventions);
         if (Err.isSuccess() and Conventions == CFConventions) {
            LOG_INFO("CFMetadataTest: Conventions metadata PASS");
         } else {
            ErrAll += Error(ErrorCode::Fail,
                            "CFMetadataTest: Conventions metadata FAIL");
         }
      }

      // ------------------------------------------------------------------
      // Test 4: attributes of a file written through IOStream
      // ------------------------------------------------------------------
      {
         IOStream::write(CFStreamName, ModelClock, true);

         int NBad = 0;
         if (DefEnv->getMyTask() == 0)
            NBad = checkFileMetadata(CFFileName);
         MPI_Bcast(&NBad, 1, MPI_INT, 0, DefEnv->getComm());

         if (NBad == 0) {
            LOG_INFO("CFMetadataTest: written file metadata PASS");
         } else {
            ErrAll += Error(ErrorCode::Fail,
                            "CFMetadataTest: written file metadata FAIL "
                            "({} problems)",
                            NBad);
         }
      }

      // ------------------------------------------------------------------
      // Finalize
      // ------------------------------------------------------------------
      Analysis::finalize();
      OceanState::clear();
      Tendencies::clear();
      Eos::destroyInstance();
      KPPMix::destroyInstance();
      VertMix::destroyInstance();
      PressureGrad::clear();
      AuxiliaryState::clear();
      Forcing::clear();
      Tracers::clear();
      IOStream::finalize();
      TimeStepper::clear();
      VertAdv::clear();
      HorzMesh::clear();
      VertCoord::clear();
      Field::clear();
      Dimension::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
   }
   Pacer::finalize();
   Kokkos::finalize();
   int RetVal = ErrAll.isFail() ? 1 : 0;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return RetVal;
}
//===----------------------------------------------------------------------===//
