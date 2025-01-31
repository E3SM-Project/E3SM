#ifndef OMEGA_AUXSTATE_H
#define OMEGA_AUXSTATE_H

#include "Config.h"
#include "DataTypes.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "Tracers.h"
#include "auxiliaryVars/KineticAuxVars.h"
#include "auxiliaryVars/LayerThicknessAuxVars.h"
#include "auxiliaryVars/TracerAuxVars.h"
#include "auxiliaryVars/VelocityDel2AuxVars.h"
#include "auxiliaryVars/VorticityAuxVars.h"

#include <memory>
#include <string>

namespace OMEGA {

/// A class for the ocean auxiliary variables.
/// The AuxiliaryState class groups together all of the
/// model auxiliary variables. It handles IO and contains methods that
/// compute the variables over a mesh.

class AuxiliaryState {
 public:
   // Member variables

   // Names of the auxiliary state and its FieldGroup
   std::string Name;
   std::string GroupName;

   // Auxiliary variables
   KineticAuxVars KineticAux;
   LayerThicknessAuxVars LayerThicknessAux;
   TracerAuxVars TracerAux;
   VorticityAuxVars VorticityAux;
   VelocityDel2AuxVars VelocityDel2Aux;

   ~AuxiliaryState();

   // Methods

   // Initialize the default auxiliary state
   static int init();

   // Create a non-default auxiliary state
   static AuxiliaryState *create(const std::string &Name, const HorzMesh *Mesh,
                                 int NVertLevels, int NTracers);

   /// Get the default auxiliary state
   static AuxiliaryState *getDefault();

   /// Get auxiliary state by name
   static AuxiliaryState *get(const std::string &Name);

   /// Remove auxiliary state by name
   static void erase(const std::string &Name);

   /// Remove all auxiliary states
   static void clear();

   /// Read and set config options
   int readConfigOptions(Config *OmegaConfig);

   // Compute all auxiliary variables needed for momentum equation
   void computeMomAux(const OceanState *State, int ThickTimeLevel,
                      int VelTimeLevel) const;

   /// Compute all auxiliary variables based on an ocean state at a given time
   /// level
   void computeAll(const OceanState *State, const Array3DReal &TracerArray,
                   int ThickTimeLevel, int VelTimeLevel) const;
   void computeAll(const OceanState *State, const Array3DReal &TracerArray,
                   int TimeLevel) const;

 private:
   AuxiliaryState(const std::string &Name, const HorzMesh *Mesh,
                  int NVertLevels, int NTracers);

   AuxiliaryState(const AuxiliaryState &) = delete;
   AuxiliaryState(AuxiliaryState &&)      = delete;

   const HorzMesh *Mesh;
   static AuxiliaryState *DefaultAuxState;
   static std::map<std::string, std::unique_ptr<AuxiliaryState>> AllAuxStates;
};

} // namespace OMEGA
#endif
