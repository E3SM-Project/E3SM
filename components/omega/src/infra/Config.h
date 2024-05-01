#ifndef OMEGA_CONFIG_H
#define OMEGA_CONFIG_H
//===-- infra/Config.h - omega configuration definitions --------*- C++ -*-===//
//
/// \file
/// \brief Omega model configuration
///
/// The Config capability manages all the input parameters needed to configure
/// and run a simulation using OMEGA. It reads an input configuration file
/// in YAML format and stores it in an efficient way for retrieval by
/// each module in OMEGA. It is based on the yaml-cpp library and relies on
/// that library for all functionality.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "MachEnv.h"
#include "mpi.h"
#include "yaml-cpp/yaml.h"

#include <string>

namespace OMEGA {

/// The Config class reads and maintains the OMEGA model configuration
/// for a given run. The configuration is a hierarchy of maps (keyword
/// value pairs). In YAML format, these are considered nodes. The
/// class is built around a yaml-cpp YAML Node and uses most of
/// its functionality while adding set/get/add/remove functions
/// for each supported type to provide a cleaner interface.
class Config {

 private:
   /// The full configuration is read once and stored as a static
   /// variable for later use.
   static Config ConfigAll;

   /// The name of the configuration that also serves as the key
   /// for the root map
   std::string Name;

   /// The YAML node containing the configuration.
   YAML::Node Node;

   /// All MPI tasks must read the initial input file but there may
   /// be a limit to the number of tasks who can simultaneously read
   /// so only a subset of tasks reads at the same time
   static const int ReadGroupSize;
   static int ReadGroupID;
   static int NumReadGroups;
   static MPI_Comm ConfigComm;

   /// We do not use an initialization routine, so we include this
   /// initialization flag so that the first Config constructed will
   /// set the master task flag for the config routines.
   static bool NotInitialized;

 public:
   // Methods

   /// Constructor that creates a Configuration with a given name and an
   /// empty YAML::Node. The node will be filled later with either a read
   /// or get operation.
   Config(const std::string &InName ///< [in] name of config (and node)
   );

   /// Constructor that creates a completely empty configuration. This
   /// should not be used generally but is needed for the static ConfigAll that
   /// is created before it can be properly initialized.
   Config();

   /// An iterator is sometimes needed to iterate through the YAML::Node
   using Iter = YAML::const_iterator;

   /// Returns an iterator to the first item in the contained YAML::Node
   Iter begin() const;

   /// Returns an iterator to the last item in the contained YAML::Node
   Iter end() const;

   /// Reads the full configuration for omega and stores in it a static
   /// YAML node for later use.  The file must be in YAML format and must be in
   /// the same directory as the executable, though Unix soft links can be used
   /// to point to a file in an alternate location.  The full configuration is
   /// only read by the master task.  The accessor functions (get/set) will
   /// broadcast or communicate the data as needed to/from other MPI tasks.
   static int readAll(std::string FileName ///< [in] input omega config file
   );

   // Retrieval (get) functions
   // -------------------------
   /// Retrieves a pointer to the full top-level Omega configuration
   static Config *getOmegaConfig();

   /// Retrieves a sub-configuration from a parent Config. We support a
   /// hierarchy of maps so that configuration for each module can
   /// be a sub-configuration under the overall Omega configuration.
   /// An empty SubConfig should have already been created with a
   /// constructor before this call so that we can extract the name and
   /// then fill the YAML node.
   /// Returns an error code that is non-zero if the group does not exist
   int get(Config &SubConfig ///< [inout] sub-configuration to retrieve
   );

   /// Retrieves the name of a configuration (or the key in the key:value
   /// pair). This is used when iterating through configuration entries so
   /// takes the iterator as an argument. This assumes the Config is a
   /// name-value pair (map) and returns an error code if not.
   static int getName(Iter ConfigIter,        ///< [in] iterator for this entry
                      std::string &ConfigName ///< [out] name of a configuration
   );
   
   /// Retrieves a 4-byte integer value from the Config based on name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of variable to get
           I4 &Value                  ///< [out] value of the variable
   );

   /// Retrieves an 8-byte integer value from the Config based on name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of variable to get
           I8 &Value                  ///< [out] value of the variable
   );

   /// Retrieves a 4-byte real value from the Config based on name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of variable to get
           R4 &Value                  ///< [out] value of the variable
   );

   /// Retrieves an 8-byte real value from the Config based on name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of variable to get
           R8 &Value                  ///< [out] value of the variable
   );

   /// Retrieves a logical/boolean value from the Config based on name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of variable to get
           bool &Value                ///< [out] value of the variable
   );

   /// Retrieves a string value from the Config based on name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of variable to get
           std::string &Value         ///< [out] value of the variable
   );

   /// Retrieves a vector of 4-byte integer values from the Config by name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of vector to get
           std::vector<I4> &Values    ///< [out] vector of values to retrieve
   );

   /// Retrieves a vector of 8-byte integer values from the Config by name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of vector to get
           std::vector<I8> &Values    ///< [out] vector of values to retrieve
   );

   /// Retrieves a vector of 4-byte real values from the Config by name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of vector to get
           std::vector<R4> &Values    ///< [out] vector of values to retrieve
   );

   /// Retrieves a vector of 8-byte real values from the Config by name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of vector to get
           std::vector<R8> &Values    ///< [out] vector of values to retrieve
   );

   /// Retrieves a vector of logical/boolean values from the Config by name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName, ///< [in] name of vector to get
           std::vector<bool> &Values  ///< [out] vector of values to retrieve
   );

   /// Retrieves a list of strings from the Config based on name
   /// Returns a non-zero error code if the variable does not exist
   int get(const std::string VarName,       ///< [in] name of list to get
           std::vector<std::string> &Values ///< [out] list of strings
   );

   // Set functions
   // -------------

   /// Resets the value of a variable in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of variable to set
           I4 Value                   ///< [in] value of the variable
   );

   /// Resets the value of a variable in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of variable to set
           I8 Value                   ///< [in] value of the variable
   );

   /// Resets the value of a variable in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of variable to set
           R4 Value                   ///< [in] value of the variable
   );

   /// Resets the value of a variable in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of variable to set
           R8 Value                   ///< [in] value of the variable
   );

   /// Resets the value of a variable in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of variable to set
           bool Value                 ///< [in] value of the variable
   );

   /// Resets the value of a variable in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of variable to set
           const std::string Value    ///< [in] value of the variable
   );

   /// Resets the value of a vector in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of vector to set
           std::vector<I4> Vector     ///< [in] vector to reset
   );

   /// Resets the value of a vector in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of vector to set
           std::vector<I8> Vector     ///< [in] vector to reset
   );

   /// Resets the value of a vector in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of vector to set
           std::vector<R4> Vector     ///< [in] vector to reset
   );

   /// Resets the value of a vector in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of vector to set
           std::vector<R8> Vector     ///< [in] vector to reset
   );

   /// Resets the value of a vector in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName, ///< [in] name of vector to set
           std::vector<bool> Vector   ///< [in] vector to reset
   );

   /// Resets the value of a list in the config
   /// Returns a non-zero error code if variable is not found.
   int set(const std::string VarName,          ///< [in] name of list to set
           const std::vector<std::string> List ///< [in] list to reset
   );

   // Add functions
   // -------------

   /// Adds a complete subconfiguration to an existing configuration
   /// Returns a non-zero error code if group already exists
   int add(const Config SubConfig ///< [in] Config to add
   );

   /// Adds a new variable to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName, ///< [in] name of variable to add
           I4 Value                   ///< [in] value of the variable
   );

   /// Adds a new variable to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName, ///< [in] name of variable to add
           I8 Value                   ///< [in] value of the variable
   );

   /// Adds a new variable to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName, ///< [in] name of variable to add
           R4 Value                   ///< [in] value of the variable
   );

   /// Adds a new variable to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName, ///< [in] name of variable to add
           R8 Value                   ///< [in] value of the variable
   );

   /// Adds a new variable to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName, ///< [in] name of variable to add
           bool Value                 ///< [in] value of the variable
   );

   /// Adds a new variable to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName, ///< [in] name of variable to add
           const std::string Value    ///< [in] value of the variable
   );

   /// Adds a new vector to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName,   ///< [in] name of vector to add
           const std::vector<I4> Values ///< [in] vector to be added
   );

   /// Adds a new vector to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName,   ///< [in] name of vector to add
           const std::vector<I8> Values ///< [in] vector to be added
   );

   /// Adds a new vector to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName,   ///< [in] name of vector to add
           const std::vector<R4> Values ///< [in] vector to be added
   );

   /// Adds a new vector to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName,   ///< [in] name of vector to add
           const std::vector<R8> Values ///< [in] vector to be added
   );

   /// Adds a new vector to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName,     ///< [in] name of vector to add
           const std::vector<bool> Values ///< [in] vector to be added
   );

   /// Adds a new list to a configuration
   /// Returns a non-zero error code if variable already exists
   int add(const std::string VarName,            ///< [in] name of list to add
           const std::vector<std::string> Values ///< [in] list to be added
   );

   // Remove function
   // ---------------

   /// Removes a variable from a configuration
   int remove(const std::string VarName ///< [in] name of variable to remove
   );

   // Inquiry (existence) function
   // ----------------------------

   /// Checks to see if a configuration group exists in the configuration
   bool existsGroup(const std::string GroupName ///< [in] name of group to find
   );

   /// Checks to see if a variable exists in the configuration
   bool existsVar(const std::string VarName ///< [in] name of variable find
   );

   // Write function
   // --------------

   /// Writes a configuration to a file
   /// Returns a non-zero error code if not successful
   int write(const std::string FileName /// name of file for config
   );

}; // end class Config

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_CONFIG_H
