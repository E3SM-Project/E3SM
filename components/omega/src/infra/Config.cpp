//===- infra/Config.cpp - omega configuration implementation ---*- C++
////-*-===//
//
// The Config capability manages all the input parameters needed to configure
// and run a simulation using OMEGA. It reads an input configuration file
// in YAML format and stores it in an efficient way for retrieval by
// each module in OMEGA. It is based on the yaml-cpp library and relies on
// that library for all functionality.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "Broadcast.h"
#include "DataTypes.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "mpi.h"
#include "yaml-cpp/yaml.h"

#include <fstream>
#include <iostream>
#include <string>

namespace OMEGA {

// Declare some of the Config static variables
const int Config::ReadGroupSize = 20;
bool Config::NotInitialized     = true;
int Config::ReadGroupID         = -1;
int Config::NumReadGroups       = 0;
MPI_Comm Config::ConfigComm;
Config Config::ConfigAll;

void Config::Initialize() {
   if (NotInitialized) {
      // Determine MPI variables
      MachEnv *DefEnv = MachEnv::getDefault();
      I4 NumTasks     = DefEnv->getNumTasks();
      I4 MyTask       = DefEnv->getMyTask();
      ConfigComm      = DefEnv->getComm();

      // Set number of groups to read and assign tasks to each
      // group in a round-robin way
      NumReadGroups = (NumTasks - 1) / ReadGroupSize + 1;
      ReadGroupID   = MyTask % NumReadGroups;

      NotInitialized = false; // now initialized for future calls
   }
}
//------------------------------------------------------------------------------
// Constructor that creates configuration with a given name and an emtpy
// YAML node that will be filled later with either a readAll or get call.
Config::Config(const std::string &InName // [in] name of config, node
) {

   // If this is the first Config created, set variables for reading the
   // input configuration file. In particular, to avoid too many MPI tasks
   // reading the same input stream, we divide the tasks into groups and
   // only one group at a time loads the full configuration from a stream
   Initialize();

   // Set the name - this is also used as the name of the root YAML node
   // (a map node) in this configuration.
   Name = InName;

   // Create an empty node
   Node[Name];

} // end Config constructor

//------------------------------------------------------------------------------
// Constructor that creates a completely empty configuration. This should
// not be used except to initialize the static variable ConfigAll that
// is created before it can be properly filled.
Config::Config() {} // end Config constructor

//------------------------------------------------------------------------------
// Reads the full configuration for omega and stores in it a static
// YAML node for later use.  The file must be in YAML format and must be in
// the same directory as the executable, though Unix soft links can be used
// to point to a file in an alternate location.  The full configuration is
// only read by the master task.  The accessor functions (get/set) will
// broadcast or communicate the data as needed to/from other MPI tasks.

void Config::readAll(const std::string &ConfigFile // [in] input YAML file
) {

   // Start a timer for config file reads
   Pacer::start("ConfigReadAll", 0);

   // Now give the full config the omega name and extract the
   // top-level omega node from the Root.
   ConfigAll.Name = "Omega";

   for (int ReadGroup = 0; ReadGroup < Config::NumReadGroups; ++ReadGroup) {

      // If it is this tasks turn, read the configuration file
      if (Config::ReadGroupID == ReadGroup) {
         // Read temporary root node
         YAML::Node RootNode = YAML::LoadFile(ConfigFile);
         // Extract Omega node
         ConfigAll.Node = RootNode["Omega"];
      }
      MPI_Barrier(ConfigComm);
   }

   Pacer::stop("ConfigReadAll", 0);
   return;

} // end Config::readAll

// Iterator for Config
//------------------------------------------------------------------------------
// Returns an iterator to the first item in the contained YAML::Node
Config::Iter Config::begin() const { return Node.begin(); }

// Returns an iterator to the last item in the contained YAML::Node
Config::Iter Config::end() const { return Node.end(); }

//------------------------------------------------------------------------------
// Retrieval (get) functions
//------------------------------------------------------------------------------
// Retrieves the top-level OMEGA config
Config *Config::getOmegaConfig() { return &ConfigAll; }

//------------------------------------------------------------------------------
// Retrieves a sub-configuration from a parent Config. An empty SubConfig must
// have already been created with a name.
// Returns a fail error code if the group does not exist
Error Config::get(Config &SubConfig // [inout] sub-configuration to retrieve
) {

   Error Err; // default success error code

   std::string GroupName = SubConfig.Name;
   if (Node[GroupName]) { // the group exists
      SubConfig.Node = Node[GroupName];
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get group: could not find group {}", GroupName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves the name of a configuration or a key in a key:value pair. It is
// used when iterating through a configuration so takes the iterator as input.
// This assumes the Config is a name-value pair (map).

Error Config::getName(Config::Iter ConfigIter, // [in] input iterator
                      std::string &ConfigName  // [out] name of a configuration
) {

   // Initialize return values
   Error Err;
   ConfigName = "unknown";

   ConfigName = ConfigIter->first.as<std::string>();

   if (ConfigName == "unknown")
      RETURN_ERROR(Err, ErrorCode::Fail, "Could not retrieve Config name");

   return Err;

} // end getName

//------------------------------------------------------------------------------
// Retrieves a 4-byte integer value from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  I4 &Value                   // [out] value of the variable
) {
   Error Err; // success error code

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<I4>();
   } else {
      // Do not modify value if not found to preserve any default value set
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get I4: could not find variable {}", VarName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves an 8-byte integer value from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  I8 &Value                   // [out] value of the variable
) {
   Error Err; // success error code

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<I8>();
   } else {
      // Do not modify value if not found to preserve any default value set
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get I8: could not find variable {}", VarName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a 4-byte real value from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  R4 &Value                   // [out] value of the variable
) {
   Error Err; // success error code

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<R4>();
   } else {
      // Do not modify value if not found to preserve any default value set
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get R4: could not find variable {}", VarName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves an 8-byte real value from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  R8 &Value                   // [out] value of the variable
) {
   Error Err; // success error code

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<R8>();
   } else {
      // Do not modify value if not found to preserve any default value set
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get R8: could not find variable {}", VarName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a logical/boolean value from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  bool &Value                 // [out] value of the variable
) {
   Error Err; // success error code

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<bool>();
   } else {
      // Do not modify value if not found to preserve any default value set
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get bool: could not find variable {}", VarName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a string value from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  std::string &Value          // [out] value of the variable
) {
   Error Err; // success error code

   // Extract variable from config on master task
   std::string TmpVal;
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<std::string>();
   } else {
      // Do not modify value if not found to preserve any default value set
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get string: could not find variable {}", VarName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a 4-byte integer vector from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  std::vector<I4> &Vector     // [out] vector to retrieve
) {
   Error Err; // success error code

   // Extract variable from config
   // First check if it exists and verify that it is a sequence node
   if (Node[VarName]) {                    // the variable exists
      YAML::Node Sequence = Node[VarName]; // extract as a node

      // if it is a sequence node, copy the sequence into a vector
      if (Sequence.IsSequence()) {

         // Determine size and resize vector to fit
         int VecSize = Sequence.size();
         Vector.resize(VecSize);

         // Now copy the sequence into the vector
         for (int i = 0; i < VecSize; ++i) {
            Vector[i] = Sequence[i].as<I4>();
         }

      } else { // not a sequence (vector) node so log an error
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "Config get I4 vector: entry not a sequence {}", VarName);
      }
   } else { // node with that name does not exist
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get I4 vector: could not find variable {}", VarName);
   }

   return Err;

} // End get I4 vector

//------------------------------------------------------------------------------
// Retrieves an 8-byte integer vector from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  std::vector<I8> &Vector     // [out] vector to retrieve
) {
   Error Err; // success error code

   // Extract variable from config
   // First check if it exists and verify that it is a sequence node
   if (Node[VarName]) {                    // the variable exists
      YAML::Node Sequence = Node[VarName]; // extract as a node

      // if it is a sequence node, copy the sequence into a vector
      if (Sequence.IsSequence()) {

         // Determine size and resize vector to fit
         int VecSize = Sequence.size();
         Vector.resize(VecSize);

         // Now copy the sequence into the vector
         for (int i = 0; i < VecSize; ++i) {
            Vector[i] = Sequence[i].as<I8>();
         }

      } else { // not a sequence (vector) so log an error
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "Config get I8 vector: entry not a sequence {}", VarName);
      }
   } else { // Node with that name does not exist
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get I8 vector: could not find variable {}", VarName);
   }

   return Err;

} // End get I8 vector

//------------------------------------------------------------------------------
// Retrieves a 4-byte real vector from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  std::vector<R4> &Vector     // [out] vector to retrieve
) {
   Error Err; // success error code

   // Extract variable from config
   // First check if it exists and verify that it is a sequence node
   if (Node[VarName]) {                    // the variable exists
      YAML::Node Sequence = Node[VarName]; // extract as a node

      // if it is a sequence node, copy the sequence into a vector
      if (Sequence.IsSequence()) {

         // Determine size and resize vector to fit
         int VecSize = Sequence.size();
         Vector.resize(VecSize);

         // Now copy the sequence into the vector
         for (int i = 0; i < VecSize; ++i) {
            Vector[i] = Sequence[i].as<R4>();
         }

      } else { // not a sequence (vector) so log an error
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "Config get R4 vector: entry not a sequence {}", VarName);
      }
   } else { // Node with that name does not exist
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get R4 vector: could not find variable {}", VarName);
   }

   return Err;

} // End get R4 vector

//------------------------------------------------------------------------------
// Retrieves an 8-byte real vector from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  std::vector<R8> &Vector     // [out] vector to retrieve
) {
   Error Err; // success error code

   // Extract variable from config
   // First check if it exists and verify that it is a sequence node
   if (Node[VarName]) {                    // the variable exists
      YAML::Node Sequence = Node[VarName]; // extract as a node

      // if it is a sequence node, copy the sequence into a vector
      if (Sequence.IsSequence()) {

         // Determine size and resize vector to fit
         int VecSize = Sequence.size();
         Vector.resize(VecSize);

         // Now copy the sequence into the vector
         for (int i = 0; i < VecSize; ++i) {
            Vector[i] = Sequence[i].as<R8>();
         }

      } else { // not a sequence (vector) so log an error
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "Config get R8 vector: entry not a sequence {}", VarName);
      }
   } else { // Node with that name does not exist
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get R8 vector: could not find variable {}", VarName);
   }

   return Err;

} // End get R8 vector

//------------------------------------------------------------------------------
// Retrieves a boolean vector from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName, // [in] name of variable to get
                  std::vector<bool> &Vector   // [out] vector to retrieve
) {
   Error Err; // success error code

   // Extract variable from config
   // First check if it exists and verify that it is a sequence node
   if (Node[VarName]) {                    // the variable exists
      YAML::Node Sequence = Node[VarName]; // extract as a node

      // if it is a sequence node, copy the sequence into a vector
      if (Sequence.IsSequence()) {

         // Determine size and resize vector to fit
         int VecSize = Sequence.size();
         Vector.resize(VecSize);

         // Now copy the sequence into the vector
         for (int i = 0; i < VecSize; ++i) {
            Vector[i] = Sequence[i].as<bool>();
         }

      } else { // not a sequence (vector) so log an error
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "Config get bool vector: entry not a sequence {}",
                      VarName);
      }
   } else { // Node with that name does not exist
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get bool vector: could not find variable {}",
                   VarName);
   }

   return Err;

} // End get bool vector

//------------------------------------------------------------------------------
// Retrieves a list of strings from the Config based on name
// Returns a fail error code if the variable does not exist
Error Config::get(const std::string &VarName,    // [in] name of variable to get
                  std::vector<std::string> &List // [out] string list retrieved
) {
   Error Err; // success error code

   // Extract variable from config
   // First check if it exists and verify that it is a sequence node
   if (Node[VarName]) {                    // the variable exists
      YAML::Node Sequence = Node[VarName]; // extract as a node

      // if it is a sequence node, copy the sequence into a vector
      if (Sequence.IsSequence()) {

         // Determine size and resize vector to fit
         int VecSize = Sequence.size();
         List.resize(VecSize);

         // Now copy the sequence into the vector
         for (int i = 0; i < VecSize; ++i) {
            List[i] = Sequence[i].as<std::string>();
         }

      } else { // not a sequence (list) so log an error
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "Config get string list: entry not a sequence {}",
                      VarName);
      }
   } else { // Node with that name does not exist
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Config get string list: could not find variable {}",
                   VarName);
   }

   return Err;

} // End get string list

//------------------------------------------------------------------------------
// Set functions
//------------------------------------------------------------------------------
// Resets the value of a variable in the config.
// Aborts if the variable is not found.
void Config::set(const std::string &VarName, // [in] name of variable to set
                 I4 Value                    // [in] new value of the variable
) {
   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      ABORT_ERROR("Config set I4: could not find variable {}", VarName);
   }
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName, // [in] name of variable to set
                 I8 Value                    // [in] new value of the variable
) {
   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      ABORT_ERROR("Config set I8: could not find variable {}", VarName);
   }
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName, // [in] name of variable to set
                 R4 Value                    // [in] new value of the variable
) {
   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      ABORT_ERROR("Config set R4: could not find variable {}", VarName);
   }
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName, // [in] name of variable to set
                 R8 Value                    // [in] new value of the variable
) {
   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      ABORT_ERROR("Config set R8: could not find variable {}", VarName);
   }
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName, // [in] name of variable to set
                 bool Value                  // [in] new value of the variable
) {
   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      ABORT_ERROR("Config set bool: could not find variable {}", VarName);
   }
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName, // [in] name of variable to set
                 const std::string &Value    // [in] new value of the variable
) {
   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      ABORT_ERROR("Config set string: could not find variable {}", VarName);
   }
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName,   // [in] name of variable to set
                 const std::vector<I4> &Vector // [in] new value of the vector
) {
   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   this->remove(VarName);
   this->add(VarName, Vector);
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName,   // [in] name of variable to set
                 const std::vector<I8> &Vector // [in] new value of the vector
) {
   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   this->remove(VarName);
   this->add(VarName, Vector);
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName,   // [in] name of variable to set
                 const std::vector<R4> &Vector // [in] new value of the vector
) {
   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   this->remove(VarName);
   this->add(VarName, Vector);
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName,   // [in] name of variable to set
                 const std::vector<R8> &Vector // [in] new value of the vector
) {
   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   this->remove(VarName);
   this->add(VarName, Vector);
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName,     // [in] name of variable to set
                 const std::vector<bool> &Vector // [in] new vector values
) {
   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   this->remove(VarName);
   this->add(VarName, Vector);
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Aborts if the variable is not found.
void Config::set(const std::string &VarName,            // [in] name of var
                 const std::vector<std::string> &Vector // [in] new values
) {
   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   this->remove(VarName);
   this->add(VarName, Vector);
}

//------------------------------------------------------------------------------
// Add functions
//------------------------------------------------------------------------------
// Adds a complete subconfiguration to an existing configuration
// Aborts if the group already exists
void Config::add(const Config SubConfig // [in] Config to add
) {
   std::string LocName = SubConfig.Name;
   if (Node[LocName]) { // the variable exists
      ABORT_ERROR("Config add group: cannot add, group {} already exists",
                  LocName);
   } else {
      Node[LocName] = SubConfig.Node;
   }
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName, // [in] name of variable to add
                 I4 Value                    // [in] value of the variable
) {
   if (Node[VarName]) { // the variable exists
      ABORT_ERROR("Config add I4: variable {} already exists use set instead",
                  VarName);
   } else {
      Node[VarName] = Value;
   }
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName, // [in] name of variable to add
                 I8 Value                    // [in] value of the variable
) {
   if (Node[VarName]) { // the variable exists
      ABORT_ERROR("Config add I8: variable {} already exists use set instead",
                  VarName);
   } else {
      Node[VarName] = Value;
   }
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName, // [in] name of variable to add
                 R4 Value                    // [in] value of the variable
) {
   if (Node[VarName]) { // the variable exists
      ABORT_ERROR("Config add R4: variable {} already exists use set instead",
                  VarName);
   } else {
      Node[VarName] = Value;
   }
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName, // [in] name of variable to add
                 R8 Value                    // [in] value of the variable
) {
   if (Node[VarName]) { // the variable exists
      ABORT_ERROR("Config add R8: variable {} already exists use set instead",
                  VarName);
   } else {
      Node[VarName] = Value;
   }
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName, // [in] name of variable to add
                 bool Value                  // [in] value of the variable
) {
   if (Node[VarName]) { // the variable exists
      ABORT_ERROR("Config add bool: variable {} already exists use set instead",
                  VarName);
   } else {
      Node[VarName] = Value;
   }
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName, // [in] name of variable to add
                 const std::string &Value    // [in] value of the variable
) {
   if (Node[VarName]) { // the variable exists
      ABORT_ERROR("Config add string: variable {} already exists - use set",
                  VarName);
   } else {
      Node[VarName] = Value;
   }
}

//------------------------------------------------------------------------------
// Adds a new I4 vector to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName,   // [in] name of vector to add
                 const std::vector<I4> &Vector // [in] vector to add
) {
   if (Node[VarName]) { // the variable already exists
      ABORT_ERROR(
          "Config add I4 vector: variable {} already exists use set instead",
          VarName);

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }
}

//------------------------------------------------------------------------------
// Adds a new I8 vector to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName,   // [in] name of vector to add
                 const std::vector<I8> &Vector // [in] vector to add
) {
   if (Node[VarName]) { // the variable already exists
      ABORT_ERROR(
          "Config add I8 vector: variable {} already exists use set instead",
          VarName);

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }
}

//------------------------------------------------------------------------------
// Adds a new R4 vector to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName,   // [in] name of vector to add
                 const std::vector<R4> &Vector // [in] vector to add
) {
   if (Node[VarName]) { // the variable already exists
      ABORT_ERROR(
          "Config add R4 vector: variable {} already exists use set instead",
          VarName);

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }
}

//------------------------------------------------------------------------------
// Adds a new R8 vector to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName,   // [in] name of vector to add
                 const std::vector<R8> &Vector // [in] vector to add
) {
   if (Node[VarName]) { // the variable already exists
      ABORT_ERROR(
          "Config add R8 vector: variable {} already exists use set instead",
          VarName);

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }
}

//------------------------------------------------------------------------------
// Adds a new bool vector to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName,     // [in] name of vector to add
                 const std::vector<bool> &Vector // [in] vector to add
) {
   if (Node[VarName]) { // the variable already exists
      ABORT_ERROR(
          "Config add bool vector: variable {} already exists - use set",
          VarName);

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         bool Tmp = Vector[i];
         Node[VarName].push_back(Tmp);
      }
   }
}

//------------------------------------------------------------------------------
// Adds a new string list to a configuration
// Aborts if the variable already exists
void Config::add(const std::string &VarName,            // [in] name of vector
                 const std::vector<std::string> &Vector // [in] vector to add
) {
   if (Node[VarName]) { // the variable already exists
      ABORT_ERROR(
          "Config add string list: variable {} already exists - use set",
          VarName);

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }
}

//------------------------------------------------------------------------------
// Remove functions
//------------------------------------------------------------------------------
// Removes a variable from a configuration
void Config::remove(const std::string &VarName // [in] name of var to remove
) {
   if (Node[VarName]) { // the variable exists
      Node.remove(VarName);
   }
}

//------------------------------------------------------------------------------
// Inquiry (existence) function
//------------------------------------------------------------------------------

// Checks to see if a configuration group exists in the configuration
bool Config::existsGroup(const std::string &GroupName // [in] group to find
) {
   bool result = false;
   if (Node[GroupName])
      result = true;
   return result;
}

// Checks to see if a variable exists in the configuration
bool Config::existsVar(const std::string &VarName // [in] name of var to find
) {
   bool result = false;
   if (Node[VarName])
      result = true;
   return result;
}

//------------------------------------------------------------------------------
// Write function
//------------------------------------------------------------------------------
// Writes a configuration to a file
// Returns a fail error code if not successful
void Config::write(const std::string &FileName /// name of file for config
) {
   int Err = 0; // temp integer flag for error broadcast

   MachEnv *DefEnv = MachEnv::getDefault();
   if (DefEnv->isMasterTask()) {
      std::ofstream Outfile(FileName);
      if (Outfile.good()) {
         Outfile << Node;
         Outfile.close();
      } else {
         Err = -1;
      }
   }

   Broadcast(Err);

   if (Err != 0)
      ABORT_ERROR("Error writing ConfigFile {}", FileName);

   return;
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
