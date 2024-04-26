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
#include "Logging.h"
#include "MachEnv.h"
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

//------------------------------------------------------------------------------
// Constructor that creates configuration with a given name and an emtpy
// YAML node that will be filled later with either a readAll or get call.
Config::Config(const std::string &InName // [in] name of config, node
) {

   // If this is the first Config created, set variables for reading the
   // input configuration file. In particular, to avoid too many MPI tasks
   // reading the same input stream, we divide the tasks into groups and
   // only one group at a time loads the full configuration from a stream

   if (NotInitialized) {
      // Determine MPI variables
      MachEnv *DefEnv = MachEnv::getDefaultEnv();
      I4 NumTasks     = DefEnv->getNumTasks();
      I4 MyTask       = DefEnv->getMyTask();
      ConfigComm      = DefEnv->getComm();

      // Set number of groups to read and assign tasks to each
      // group in a round-robin way
      NumReadGroups = (NumTasks - 1) / ReadGroupSize + 1;
      ReadGroupID   = MyTask % NumReadGroups;

      NotInitialized = false; // now initialized for future calls
   }

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

int Config::readAll(const std::string ConfigFile // [in] input YAML config file
) {
   int Err = 0;

   // Now give the full config the omega name and extract the
   // top-level omega node from the Root.
   ConfigAll.Name = "omega";

   for (int ReadGroup = 0; ReadGroup < Config::NumReadGroups; ++ReadGroup) {

      // If it is this tasks turn, read the configuration file
      if (Config::ReadGroupID == ReadGroup) {
         // Read temporary root node
         YAML::Node RootNode = YAML::LoadFile(ConfigFile);
         // Extract Omega node
         ConfigAll.Node = RootNode["omega"];
      }
      MPI_Barrier(ConfigComm);
   }

   return Err;

} // end Config::readAll

//------------------------------------------------------------------------------
// Retrieval (get) functions
//------------------------------------------------------------------------------
// Retrieves the top-level OMEGA config
Config *Config::getOmegaConfig() { return &ConfigAll; }

//------------------------------------------------------------------------------
// Retrieves a sub-configuration from a parent Config. An empty SubConfig must
// have already been created with a name.
// Returns an error code that is non-zero if the group does not exist
int Config::get(Config &SubConfig // [inout] sub-configuration to retrieve
) {

   int Err = 0;

   std::string GroupName = SubConfig.Name;
   if (Node[GroupName]) { // the group exists
      SubConfig.Node = Node[GroupName];
   } else {
      LOG_ERROR("Config get group: could not find group {}", GroupName);
      Err = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a 4-byte integer value from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                I4 &Value                  // [out] value of the variable
) {
   int Err = 0;

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<OMEGA::I4>();
   } else {
      LOG_ERROR("Config get I4: could not find variable {}", VarName);
      Value = -999;
      Err   = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves an 8-byte integer value from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                I8 &Value                  // [out] value of the variable
) {
   int Err = 0;

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<OMEGA::I8>();
   } else {
      LOG_ERROR("Config get I8: could not find variable {}", VarName);
      Value = -999;
      Err   = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a 4-byte real value from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                R4 &Value                  // [out] value of the variable
) {
   int Err = 0;

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<OMEGA::R4>();
   } else {
      LOG_ERROR("Config get R4: could not find variable {}", VarName);
      Value = -999.999;
      Err   = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves an 8-byte real value from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                R8 &Value                  // [out] value of the variable
) {
   int Err = 0;

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<OMEGA::R8>();
   } else {
      LOG_ERROR("Config get R8: could not find variable {}", VarName);
      Value = -99999.999;
      Err   = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a logical/boolean value from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                bool &Value                // [out] value of the variable
) {
   int Err = 0;

   // Extract variable from config
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<bool>();
   } else {
      LOG_ERROR("Config get bool: could not find variable {}", VarName);
      Value = false;
      Err   = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a string value from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                std::string &Value         // [out] value of the variable
) {
   int Err = 0;

   // Extract variable from config on master task
   std::string TmpVal;
   if (Node[VarName]) { // the variable exists
      Value = Node[VarName].as<std::string>();
   } else {
      LOG_ERROR("Config get string: could not find variable {}", VarName);
      Value = "ConfigError";
      Err   = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Retrieves a 4-byte integer vector from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                std::vector<I4> &Vector    // [out] vector of values retrieved
) {
   int Err = 0;

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
         LOG_ERROR("Config get I4 vector: entry not a sequence {}", VarName);
         Err = -2;
      }
   } else { // node with that name does not exist
      LOG_ERROR("Config get I4 vector: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;

} // End get I4 vector

//------------------------------------------------------------------------------
// Retrieves an 8-byte integer vector from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                std::vector<I8> &Vector    // [out] vector of values retrieved
) {
   int Err = 0;

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
         LOG_ERROR("Config get I8 vector: entry not a sequence {}", VarName);
         Err = -2;
      }
   } else { // Node with that name does not exist
      LOG_ERROR("Config get I8 vector: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;

} // End get I8 vector

//------------------------------------------------------------------------------
// Retrieves a 4-byte real vector from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                std::vector<R4> &Vector    // [out] vector of values retrieved
) {
   int Err = 0;

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
         LOG_ERROR("Config get R4 vector: entry not a sequence {}", VarName);
         Err = -2;
      }
   } else { // Node with that name does not exist
      LOG_ERROR("Config get R4 vector: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;

} // End get R4 vector

//------------------------------------------------------------------------------
// Retrieves an 8-byte real vector from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                std::vector<R8> &Vector    // [out] vector of values retrieved
) {
   int Err = 0;

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
         LOG_ERROR("Config get R8 vector: entry not a sequence {}", VarName);
         Err = -2;
      }
   } else { // Node with that name does not exist
      LOG_ERROR("Config get R8 vector: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;

} // End get R8 vector

//------------------------------------------------------------------------------
// Retrieves a boolean vector from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName, // [in] name of variable to get
                std::vector<bool> &Vector  // [out] vector of values retrieved
) {
   int Err = 0;

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
         LOG_ERROR("Config get bool vector: entry not a sequence {}", VarName);
         Err = -2;
      }
   } else { // Node with that name does not exist
      LOG_ERROR("Config get bool vector: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;

} // End get bool vector

//------------------------------------------------------------------------------
// Retrieves a list of strings from the Config based on name
// Returns a non-zero error code if the variable does not exist
int Config::get(const std::string VarName,     // [in] name of variable to get
                std::vector<std::string> &List // [out] string list retrieved
) {
   int Err     = 0;
   int VecSize = 0;

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
         LOG_ERROR("Config get string list: entry not a sequence {}", VarName);
         Err = -2;
      }
   } else { // Node with that name does not exist
      LOG_ERROR("Config get string list: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;

} // End get string list

//------------------------------------------------------------------------------
// Set functions
//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                I4 Value                   // [in] new value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      LOG_ERROR("Config set I4: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                I8 Value                   // [in] new value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      LOG_ERROR("Config set I8: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                R4 Value                   // [in] new value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      LOG_ERROR("Config set R4: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                R8 Value                   // [in] new value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      LOG_ERROR("Config set R8: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                bool Value                 // [in] new value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      LOG_ERROR("Config set bool: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a variable in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                std::string Value          // [in] new value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      Node[VarName] = Value;
   } else {
      LOG_ERROR("Config set string: could not find variable {}", VarName);
      Err = -1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                std::vector<I4> Vector     // [in] new value of the vector
) {
   int Err = 0;

   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   Err = this->remove(VarName);
   if (Err != 0)
      LOG_ERROR("Config set I4 vector: could not remove old variable {}",
                VarName);
   Err = this->add(VarName, Vector);
   if (Err != 0)
      LOG_ERROR("Config set I4 vector: could not re-add variable {}", VarName);

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                std::vector<I8> Vector     // [in] new value of the vector
) {
   int Err = 0;

   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   Err = this->remove(VarName);
   if (Err != 0)
      LOG_ERROR("Config set I8 vector: could not remove old variable {}",
                VarName);

   Err = this->add(VarName, Vector);
   if (Err != 0)
      LOG_ERROR("Config set I8 vector: could not re-add variable {}", VarName);

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                std::vector<R4> Vector     // [in] new value of the vector
) {
   int Err = 0;

   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   Err = this->remove(VarName);
   if (Err != 0)
      LOG_ERROR("Config set R4 vector: could not remove old variable {}",
                VarName);

   Err = this->add(VarName, Vector);
   if (Err != 0)
      LOG_ERROR("Config set R4 vector: could not re-add variable {}", VarName);

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                std::vector<R8> Vector     // [in] new value of the vector
) {
   int Err = 0;

   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   Err = this->remove(VarName);
   if (Err != 0)
      LOG_ERROR("Config set R8 vector: could not remove old variable {}",
                VarName);

   Err = this->add(VarName, Vector);
   if (Err != 0)
      LOG_ERROR("Config set R8 vector: could not re-add variable {}", VarName);

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName, // [in] name of variable to set
                std::vector<bool> Vector   // [in] new value of the vector
) {
   int Err = 0;

   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   Err = this->remove(VarName);
   if (Err != 0)
      LOG_ERROR("Config set bool vector: could not remove old variable {}",
                VarName);

   Err = this->add(VarName, Vector);
   if (Err != 0)
      LOG_ERROR("Config set bool vector: could not re-add variable {}",
                VarName);

   return Err;
}

//------------------------------------------------------------------------------
// Resets the value of a vector/list in the config
// Returns a non-zero error code if variable is not found.
int Config::set(const std::string VarName,      // [in] name of variable to set
                std::vector<std::string> Vector // [in] new value of the vector
) {
   int Err = 0;

   // Because the vector length may also change, it is best to just
   // remove and replace the entire vector
   Err = this->remove(VarName);
   if (Err != 0)
      LOG_ERROR("Config set string list: could not remove old variable {}",
                VarName);

   Err = this->add(VarName, Vector);
   if (Err != 0)
      LOG_ERROR("Config set string list: could not re-add variable {}",
                VarName);

   return Err;
}

//------------------------------------------------------------------------------
// Add functions
//------------------------------------------------------------------------------
// Adds a complete subconfiguration to an existing configuration
// Returns a non-zero error code if group already exists
int Config::add(const Config SubConfig // [in] Config to add
) {
   int Err = 0;

   std::string LocName = SubConfig.Name;
   if (Node[LocName]) { // the variable exists
      LOG_ERROR("Config add group: cannot add, group {} already exists",
                LocName);
      Err = -1;
   } else {
      Node[LocName] = SubConfig.Node;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of variable to add
                I4 Value                   // [in] value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      LOG_ERROR("Config add I4: variable {} already exists use set instead",
                VarName);
      Err = -1;
   } else {
      Node[VarName] = Value;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of variable to add
                I8 Value                   // [in] value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      LOG_ERROR("Config add I8: variable {} already exists use set instead",
                VarName);
      Err = -1;
   } else {
      Node[VarName] = Value;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of variable to add
                R4 Value                   // [in] value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      LOG_ERROR("Config add R4: variable {} already exists use set instead",
                VarName);
      Err = -1;
   } else {
      Node[VarName] = Value;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of variable to add
                R8 Value                   // [in] value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      LOG_ERROR("Config add R8: variable {} already exists use set instead",
                VarName);
      Err = -1;
   } else {
      Node[VarName] = Value;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of variable to add
                bool Value                 // [in] value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      LOG_ERROR("Config add bool: variable {} already exists use set instead",
                VarName);
      Err = -1;
   } else {
      Node[VarName] = Value;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new variable to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of variable to add
                std::string Value          // [in] value of the variable
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      LOG_ERROR("Config add string: variable {} already exists - use set",
                VarName);
      Err = -1;
   } else {
      Node[VarName] = Value;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new I4 vector to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of vector to add
                std::vector<I4> Vector     // [in] vector to add
) {
   int Err = 0;

   if (Node[VarName]) { // the variable already exists
      LOG_ERROR(
          "Config add I4 vector: variable {} already exists use set instead",
          VarName);
      Err = -1;

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new I8 vector to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of vector to add
                std::vector<I8> Vector     // [in] vector to add
) {
   int Err = 0;

   if (Node[VarName]) { // the variable already exists
      LOG_ERROR(
          "Config add I8 vector: variable {} already exists use set instead",
          VarName);
      Err = -1;

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new R4 vector to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of vector to add
                std::vector<R4> Vector     // [in] vector to add
) {
   int Err = 0;

   if (Node[VarName]) { // the variable already exists
      LOG_ERROR(
          "Config add R4 vector: variable {} already exists use set instead",
          VarName);
      Err = -1;

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new R8 vector to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of vector to add
                std::vector<R8> Vector     // [in] vector to add
) {
   int Err = 0;

   if (Node[VarName]) { // the variable already exists
      LOG_ERROR(
          "Config add R8 vector: variable {} already exists use set instead",
          VarName);
      Err = -1;

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new bool vector to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName, // [in] name of vector to add
                std::vector<bool> Vector   // [in] vector to add
) {
   int Err = 0;

   if (Node[VarName]) { // the variable already exists
      LOG_ERROR("Config add bool vector: variable {} already exists - use set",
                VarName);
      Err = -1;

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         bool Tmp = Vector[i];
         Node[VarName].push_back(Tmp);
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Adds a new string list to a configuration
// Returns a non-zero error code if variable already exists
int Config::add(const std::string VarName,      // [in] name of vector to add
                std::vector<std::string> Vector // [in] vector to add
) {
   int Err = 0;

   if (Node[VarName]) { // the variable already exists
      LOG_ERROR("Config add string list: variable {} already exists - use set",
                VarName);
      Err = -1;

   } else { // create new sequence node for vector

      int VecSize = Vector.size();
      for (int i = 0; i < VecSize; ++i) {
         Node[VarName].push_back(Vector[i]);
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Remove functions
//------------------------------------------------------------------------------
// Removes a variable from a configuration
int Config::remove(const std::string VarName // [in] name of variable to remove
) {
   int Err = 0;

   if (Node[VarName]) { // the variable exists
      Node.remove(VarName);
   }

   return Err;
}

//------------------------------------------------------------------------------
// Inquiry (existence) function
//------------------------------------------------------------------------------

// Checks to see if a configuration group exists in the configuration
bool Config::existsGroup(std::string GroupName // [in] name of group to find
) {
   bool result = false;
   if (Node[GroupName])
      result = true;
   return result;
}

// Checks to see if a variable exists in the configuration
bool Config::existsVar(std::string VarName // [in] name of variable to find
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
// Returns a non-zero error code if not successful
int Config::write(std::string FileName /// name of file for config
) {
   int Err = 0;

   MachEnv *DefEnv = MachEnv::getDefaultEnv();
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
   return Err;
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
