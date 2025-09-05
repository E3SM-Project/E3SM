//===-- ocn/VertMix.cpp - Vertical Mixing Coefficients -----------*- C++ -*-===//
//
// The VertMix class is responsible for managing the calculation of the
// vertical diffusivity and viscosity needed for the vertical mixing. 
// It currently has background, convective, and shear mixing options, and 
// they can be additively combined, depending on configuration options. It 
// contains arrays that store the vertical top-of-cell diffusivity and 
// viscosity values for each cell and vertical level.
//
//===----------------------------------------------------------------------===//

#include "VertMix.h"
#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {
      
PPShearMix::PPShearMix(const HorzMesh *Mesh) 
    : NVertLevels(Mesh->NVertLevels), NEdgesOnCell(Mesh->NEdgesOnCell), 
      EdgesOnCell(Mesh->EdgesOnCell), DvEdge(Mesh->DvEdge), DcEdge(Mesh->DvEdge), 
      AreaCell(Mesh->AreaCell) {} 

ConvectiveMix::ConvectiveMix() {} 

/// Constructor for VertMix
VertMix::VertMix(const std::string &Name_,  ///< [in] Name for VertMix object
         const HorzMesh *Mesh,     ///< [in] Horizontal mesh
         int NVertLevels           ///< [in] Number of vertical levels
         )
    : ComputeVertMixShear(Mesh), ComputeVertMixConv() {
   VertDiff  = Array2DReal("VertDiff", Mesh->NCellsAll, NVertLevels);
   VertVisc  = Array2DReal("VertVisc", Mesh->NCellsAll, NVertLevels);
   NCellsAll = Mesh->NCellsAll;
   Name      = Name_;

   defineFields();
}

/// Destructor for VertMix 
VertMix::~VertMix() {}

/// Instance management
VertMix *VertMix::Instance = nullptr;

/// Get instance of VertMix
VertMix* VertMix::getInstance() { return Instance; }

/// Destroy instance of VertMix
void VertMix::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initializes the VertMix (Vertical Mixing Coefficients) class and its options.
/// it ASSUMES that HorzMesh was initialized and initializes the VertMix class by
/// using the default mesh, reading the config file, and setting parameters
/// for the background, convective, and/or shear mixing routines.
/// Returns 0 on success, or an error code if any required option is missing.
void VertMix::init() {

   if (!Instance) {
      Instance = new VertMix("Default", HorzMesh::getDefault(),
                         HorzMesh::getDefault()->NVertLevels);
   }

   Error Err; // error code
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   I4 NVertLevels        = DefHorzMesh->NVertLevels;

   /// Retrieve default VertMix
   VertMix *vertmix = VertMix::getInstance();

   /// Get VertMixConfig group from Omega config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config VertMixConfig("VertMix");
   Err += OmegaConfig->get(VertMixConfig);
   CHECK_ERROR_ABORT(Err, "VertMix::init: VertMix group not found in Config");

   /// Get VertMixType from VertMixConfig
   /// and set the VertMixChoice accordingly
   std::string VertMixTypeStr;
   Err += VertMixConfig.get("VertMixType", VertMixTypeStr);
   CHECK_ERROR_ABORT(Err, "VertMix::init: VertMixType subgroup not found in VertMixConfig");

   if (VertMixTypeStr == "PP" or VertMixTypeStr == "pp") {
      vertmix->VertMixChoice = VertMixType::PP;
      LOG_INFO("VertMix::init: Using Pacanowski and Philander (1981) vertical mixing.");
   } else if (VertMixTypeStr == "KPP" or VertMixTypeStr == "kpp") {
      vertmix->VertMixChoice = VertMixType::KPP;
      LOG_INFO("VertMix::init: Using K-Profile Parameterization (Large et al., 1994) vertical mixing.");
   } else {
      LOG_ERROR("VertMix::init: Invalid VertMixType specified in configuration: " + VertMixTypeStr);
   }

   /// Get Background from VertMixConfig
   /// and set associated parameters
   Config BackConfig("Background");
   Err += VertMixConfig.get(BackConfig);
   CHECK_ERROR_ABORT(Err,"VertMix::init: Background subgroup not found in VertMixConfig");

   /// Get diffusivity and viscosity parameters
   Err += BackConfig.get("Viscosity", vertmix->BackVisc);
   CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Background:Viscosity not found in BackConfig");

   Err += BackConfig.get("Diffusivity", vertmix->BackDiff);
   CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Background:Diffusivity not found in BackConfig");

   /// Get Convective from VertMixConfig
   Config ConvConfig("Convective");
   Err += VertMixConfig.get(ConvConfig);
   CHECK_ERROR_ABORT(Err,"VertMix::init: Convective subgroup not found in VertMixConfig");

   /// Get convective diffusivity and viscosity parameters
   Err += ConvConfig.get("Enable", vertmix->ComputeVertMixConv.Enabled);
   CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Convective:Enable not found in ConvConfig");

   if (!vertmix->ComputeVertMixConv.Enabled) {
      LOG_INFO("VertMix::init: Convective mixing is disabled.");
   } else {
      LOG_INFO("VertMix::init: Convective mixing is enabled.");
      Err += ConvConfig.get("Diffusivity", vertmix->ComputeVertMixConv.ConvDiff);
      CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Convective:Diffusivity not found in ConvConfig");

      Err += ConvConfig.get("TriggerBVF", vertmix->ComputeVertMixConv.ConvTriggerBVF);
      CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Convective:TriggerBVF not found in ConvConfig");
   }

   /// Get Shear from VertMixConfig
   Config ShearConfig("Shear");
   Err += VertMixConfig.get(ShearConfig);
   CHECK_ERROR_ABORT(Err,"VertMix::init: Shear subgroup not found in VertMixConfig");

   /// Get shear diffusivity and viscosity parameters
   Err += ShearConfig.get("Enable", vertmix->ComputeVertMixShear.Enabled);
   CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Shear:Enable not found in ShearConfig");

   if (!vertmix->ComputeVertMixShear.Enabled) {
      LOG_INFO("VertMix::init: Shear mixing is disabled.");
   } else {
      LOG_INFO("VertMix::init: Shear mixing is enabled.");
      Err += ShearConfig.get("NuZero", vertmix->ComputeVertMixShear.ShearNuZero);
      CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Shear:NuZero not found in ShearConfig");

      Err += ShearConfig.get("Alpha", vertmix->ComputeVertMixShear.ShearAlpha);
      CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Shear:Alpha not found in ShearConfig");

      Err += ShearConfig.get("Exponent", vertmix->ComputeVertMixShear.ShearExponent);
      CHECK_ERROR_ABORT(Err,"VertMix::init: Parameter Shear:Exponent not found in ShearConfig");
   }
} // end init

/// Compute diffusivity and viscosity for all cells/levels (no displacement)
void VertMix::computeVertMix(const Array2DReal &NormalVelocity,
                             const Array2DReal &TangentialVelocity,
                             const Array2DReal &BruntVaisalaFreq,
                             const Array2DReal &ZMid) {
   OMEGA_SCOPE(LocVertDiff, VertDiff); /// Create a local view for computation
   OMEGA_SCOPE(LocVertVisc, VertVisc); /// Create a local view for computation
   OMEGA_SCOPE(LocComputeVertMixConv,
               ComputeVertMixConv);    /// Local view for PP VertMix computation
   OMEGA_SCOPE(LocComputeVertMixShear,
               ComputeVertMixShear);   /// Local view for PP VertMix computation

   deepCopy(LocVertDiff, BackDiff);
   deepCopy(LocVertVisc, BackVisc);

   /// Dispatch to the correct VertMix calculation
   if (VertMixChoice == VertMixType::PP) {
      if (LocComputeVertMixShear.Enabled && LocComputeVertMixConv.Enabled) {
         parallelFor(
            "VertMix-PP", {NCellsAll, NVertLevels},
            KOKKOS_LAMBDA(I4 ICell, I4 K) {
               //LocComputeVertMixBack(LocVertDiff, LocVertVisc, ICell, K);
               LocComputeVertMixConv(LocVertDiff, LocVertVisc, ICell, K, 
                  BruntVaisalaFreq);
               LocComputeVertMixShear(LocVertDiff, LocVertVisc, ICell, K, 
                  NormalVelocity, TangentialVelocity,
                  BruntVaisalaFreq, ZMid);
         });
      } else if (LocComputeVertMixShear.Enabled) {
         parallelFor(
            "VertMix-PP", {NCellsAll, NVertLevels},
            KOKKOS_LAMBDA(I4 ICell, I4 K) {
               //LocComputeVertMixBack(LocVertDiff, LocVertVisc, ICell, K);
               LocComputeVertMixShear(LocVertDiff, LocVertVisc, ICell, K, 
                  NormalVelocity, TangentialVelocity,
                  BruntVaisalaFreq, ZMid);
         });
      } else if (LocComputeVertMixConv.Enabled) {
         parallelFor(
            "VertMix-PP", {NCellsAll, NVertLevels},
            KOKKOS_LAMBDA(I4 ICell, I4 K) {
               //LocComputeVertMixBack(LocVertDiff, LocVertVisc, ICell, K);
               LocComputeVertMixConv(LocVertDiff, LocVertVisc, ICell, K, 
                  BruntVaisalaFreq);
         });
      }
   } else if (VertMixChoice == VertMixType::KPP) {
      LOG_ERROR("VertMix: VertMixType = KPP is not supported yet.");
      //parallelFor(
      //   "VertMix-KPP", {NCellsAll, NVertLevels},
      //   KOKKOS_LAMBDA(I4 ICell, I4 K) {
      //      LocComputeVertMixCoeffKPP(LocVertDiff, LocVertVisc, ICell, K, 
      //         NormalVelocity, TangentialVelocity, BruntVaisalaFreq);
      //   });
   }
}

/// Define IO fields and metadata for output
void VertMix::defineFields() {

   I4 Err = 0;

   /// Set field names (append Name if not default)
   VertDiffFldName = "VertDiff";
   VertViscFldName = "VertVisc";
   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
   }

   /// Create fields for state variables
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";

   /// Create and register the Diffusivity field
   auto VertDiffField =
       Field::create(VertDiffFldName,                       // Field name
                     "Vertical diffusivity at center and"
                     " top of cell",                        // Long Name
                     "m2 s-1",                              // Units
                     "vertical_diffusivity",                // CF-ish Name
                     0.0,                                   // Min valid value
                     9.99E+30,                              // Max valid value
                     -9.99E+30,                             // Scalar used for undefined entries
                     NDims,                                 // Number of dimensions
                     DimNames                               // Dimension names
       );
   /// Create and register the VertVisc field
   auto VertViscField =
       Field::create(VertViscFldName,                     // Field name
                     "Vertical viscosity at center and"
                     " top of cell",                      // Long Name
                     "m2 s-1",                            // Units
                     "vertical_viscosity",                // CF-ish Name
                     0.0,                                 // Min valid value
                     9.99E+30,                            // Max valid value
                     -9.99E+30,                           // Scalar used for undefined entried
                     NDims,                               // Number of dimensions
                     DimNames                             // Dimension names
       );

   // Create a field group for the vertmix-specific state fields
   VertMixGroupName = "VertMix";
   if (Name != "Default") {
      VertMixGroupName.append(Name);
   }
   auto VertMixGroup = FieldGroup::create(VertMixGroupName);

   // Add fields to the VertMix group
   Err = VertMixGroup->addField(VertDiffFldName);
   if (Err != 0)
      LOG_ERROR("VertMix::defineFields: Error adding {} to field group {}", 
                VertDiffFldName, VertMixGroupName);
   Err = VertMixGroup->addField(VertViscFldName);
   if (Err != 0)
      LOG_ERROR("VertMix::defineFields: Error adding {} to field group {}", 
                VertViscFldName, VertMixGroupName);

   // Attach Kokkos views to the fields
   Err = VertDiffField->attachData<Array2DReal>(VertDiff);
   if (Err != 0)
      LOG_ERROR("VertMix::defineFields: Error attaching data array to field {}",
                VertDiffFldName);
   Err = VertViscField->attachData<Array2DReal>(VertVisc);
   if (Err != 0)
      LOG_ERROR("VertMix::defineFields: Error attaching data array to field {}", 
                VertViscFldName);

} // end defineIOFields

} // namespace OMEGA