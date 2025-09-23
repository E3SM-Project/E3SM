//===-- ocn/Eos.cpp - Equation of State ------------------*- C++ -*-===//
//
// The Eos class is responsible for managing the equation of state. It
// has a linear EOS and TEOS-10 EOS option, which is determined at
// initialization. It contains arrays that store the specific volume and
// displaced specific volume data.
//
//===----------------------------------------------------------------------===//

#include "Eos.h"
#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {

/// Constructor for Teos10Eos
Teos10Eos::Teos10Eos(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {
   SpecVolPCoeffs = Array2DReal("SpecVolPCoeffs", 6, VecLength);
}

/// Constructor for LinearEos
LinearEos::LinearEos(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

/// Constructor for Teos10 Brunt-Vaisala frequency
Teos10BruntVaisalaFreq::Teos10BruntVaisalaFreq(const HorzMesh *Mesh) 
    : NVertLevels(Mesh->NVertLevels), LatCell(Mesh->LatCellH) {}

/// Constructor for Simple Brunt-Vaisala frequency
SimpleBruntVaisalaFreq::SimpleBruntVaisalaFreq(int NVertLevels) : NVertLevels(NVertLevels) {}

/// Constructor for Eos
Eos::Eos(const std::string &Name, ///< [in] Name for eos object
         const HorzMesh *Mesh,    ///< [in] Horizontal mesh
         const VertCoord *VCoord  ///< [in] Vertical coordinate
         )
    : ComputeSpecVolLinear(VCoord), ComputeSpecVolTeos10(VCoord), 
      ComputeBruntVaisalaFreqSimple(VCoord), ComputeBruntVaisalaFreqTeos10(Mesh),
      Name(Name), Mesh(Mesh), VCoord(VCoord) {
   SpecVol = Array2DReal("SpecVol", Mesh->NCellsAll, VCoord->NVertLayers);
   SpecVolDisplaced =
       Array2DReal("SpecVolDisplaced", Mesh->NCellsAll, VCoord->NVertLayers);
   BruntVaisalaFreq =
       Array2DReal("BruntVaisalaFreq", Mesh->NCellsAll, VCoord->NVertLayers);

   defineFields();
}

/// Destructor for Eos
Eos::~Eos() {}

/// Instance management
Eos *Eos::Instance = nullptr;

/// Get instance of Eos
Eos *Eos::getInstance() { return Instance; }

/// Destroy instance of Eos
void Eos::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initializes the Eos (Equation of State) class and its options.
/// it ASSUMES that HorzMesh was initialized and initializes the Eos class by
/// using the default mesh, reading the config file, and setting parameters
/// for either a Linear or TEOS-10 equation.
void Eos::init() {

   if (!Instance) {
      Instance =
          new Eos("Default", HorzMesh::getDefault(), VertCoord::getDefault());
   }

   Error Err; // error code

   /// Retrieve default eos
   Eos *eos = Eos::getInstance();

   /// Get EosConfig group from Omega config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config EosConfig("Eos");
   Err += OmegaConfig->get(EosConfig);
   CHECK_ERROR_ABORT(Err, "Eos::init: Eos group not found in Config");

   /// Get EosType from EosConfig
   /// and set the EosChoice accordingly
   std::string EosTypeStr;
   Err += EosConfig.get("EosType", EosTypeStr);
   CHECK_ERROR_ABORT(Err, "Eos::init: EosType subgroup not found in EosConfig");

   /// Set EosChoice and parameters based on EosTypeStr
   if (EosTypeStr == "Linear" or EosTypeStr == "linear") {
      Config EosLinConfig("Linear");
      Err += EosConfig.get(EosLinConfig);

      eos->EosChoice = EosType::LinearEos;

      CHECK_ERROR_ABORT(Err,
                        "Eos::init: Linear subgroup not found in EosConfig");
      Err += EosLinConfig.get("DRhoDT", eos->ComputeSpecVolLinear.DRhodT);
      CHECK_ERROR_ABORT(
          Err, "Eos::init: Parameter Linear:DRhodT not found in EosLinConfig");

      Err += EosLinConfig.get("DRhoDS", eos->ComputeSpecVolLinear.DRhodS);
      CHECK_ERROR_ABORT(
          Err, "Eos::init: Parameter Linear:DRhodS not found in EosLinConfig");

      Err += EosLinConfig.get("RhoT0S0", eos->ComputeSpecVolLinear.RhoT0S0);
      CHECK_ERROR_ABORT(
          Err, "Eos::init: Parameter Linear:RhoT0S0 not found in EosLinConfig");
   } else if ((EosTypeStr == "teos10") or (EosTypeStr == "teos-10") or
              (EosTypeStr == "TEOS-10")) {
      eos->EosChoice = EosType::Teos10Eos;
   } else {
      LOG_ERROR("Eos::init: Unknown EosType requested");
   }
} // end init

/// Compute specific volume for all cells/layers (no displacement)
void Eos::computeSpecVol(const Array2DReal &ConservTemp,
                         const Array2DReal &AbsSalinity,
                         const Array2DReal &Pressure) {
   OMEGA_SCOPE(LocSpecVol, SpecVol); /// Create a local view for computation
   OMEGA_SCOPE(LocComputeSpecVolLinear,
               ComputeSpecVolLinear); /// Local view for linear EOS computation
   OMEGA_SCOPE(LocComputeSpecVolTeos10,
               ComputeSpecVolTeos10); /// Local view for TEOS-10 computation
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   deepCopy(LocSpecVol, 0); /// Initialize local specific volume to zero

   I4 KDisp = 0; /// No displacement in this case

   /// Dispatch to the correct EOS calculation
   if (EosChoice == EosType::LinearEos) {
      parallelForOuter(
          "eos-linear", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolLinear(LocSpecVol, ICell, KChunk,
                                            ConservTemp, AbsSalinity);
                 });
          });
   } else if (EosChoice == EosType::Teos10Eos) {
      parallelForOuter(
          "eos-teos10", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolTeos10(LocSpecVol, ICell, KChunk,
                                            ConservTemp, AbsSalinity, Pressure,
                                            KDisp);
                 });
          });
   }
}

/// Compute displaced specific volume (for vertical displacement)
void Eos::computeSpecVolDisp(const Array2DReal &ConservTemp,
                             const Array2DReal &AbsSalinity,
                             const Array2DReal &Pressure, I4 KDisp) {
   OMEGA_SCOPE(LocSpecVolDisplaced,
               SpecVolDisplaced); /// Local view for computation
   OMEGA_SCOPE(LocComputeSpecVolLinear,
               ComputeSpecVolLinear); /// Local view for linear EOS computation
   OMEGA_SCOPE(LocComputeSpecVolTeos10,
               ComputeSpecVolTeos10); /// Local view for TEOS-10 computation
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   deepCopy(LocSpecVolDisplaced,
            0); /// Initialize local specific volume to zero

   /// Dispatch to the correct EOS calculation
   /// If EosChoice is Linear, the displaced specific
   /// volume is the same as the specific volume
   if (EosChoice == EosType::LinearEos) {
      LOG_INFO("Eos::computeSpecVolDisp called with Linear EOS. "
               "SpecVol is independent of pressure/depth, so the "
               "displaced value will be the same as SpecVol.");
      parallelForOuter(
          "eos-linear", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolLinear(LocSpecVolDisplaced, ICell, KChunk,
                                            ConservTemp, AbsSalinity);
                 });
          });
   } else if (EosChoice == EosType::Teos10Eos) {
      parallelForOuter(
          "eos-teos10", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolTeos10(LocSpecVolDisplaced, ICell, KChunk,
                                            ConservTemp, AbsSalinity, Pressure,
                                            KDisp);
                 });
          });
   }
}

/// Compute Brunt-Vaisala frequency (NSquared) for all cells/levels
void Eos::computeBruntVaisalaFreq(const Array2DReal &ConservTemp,
                                  const Array2DReal &AbsSalinity,
                                  const Array2DReal &Pressure,
                                  const Array2DReal &SpecVol,
                                  const Array2DReal &ZMid) {
   OMEGA_SCOPE(LocBruntVaisalaFreq,
               BruntVaisalaFreq); /// Local view for computation
   OMEGA_SCOPE(LocComputeBruntVaisalaFreqSimple,
               ComputeBruntVaisalaFreqSimple); /// Local view for simple computation
   OMEGA_SCOPE(LocComputeBruntVaisalaFreqTeos10,
               ComputeBruntVaisalaFreqTeos10); /// Local view for TEOS-10 computation
   OMEGA_SCOPE(LocComputeSpecVolTeos10,
               ComputeSpecVolTeos10); /// Local view for TEOS-10 computation
   deepCopy(LocBruntVaisalaFreq,
            0); /// Initialize local Brunt-Vaisala frequency to zero
   //OMEGA_SCOPE(LocLatCell, Mesh->LatCellH);

   /// Dispatch to the correct Brunt-Vaisala frequency calculation
   if (EosChoice == EosType::LinearEos) {
      /// If Linear EOS, use simple Brunt-Vaisala frequency calculation
      parallelFor(
          "bvf-linear", {NCellsAll},
          KOKKOS_LAMBDA(I4 ICell) {
             LocComputeBruntVaisalaFreqSimple(LocBruntVaisalaFreq, ICell,
                                     SpecVol, ZMid);
          });
   } else if (EosChoice == EosType::Teos10Eos) {
      /// If TEOS-10 EOS, use TEOS-10 Brunt-Vaisala frequency calculation
      parallelFor(
          "bvf-teos10", {NCellsAll},
          KOKKOS_LAMBDA(I4 ICell) {
             /// Compute Brunt-Vaisala frequency
             LocComputeBruntVaisalaFreqTeos10(LocBruntVaisalaFreq, ICell,
                                     ConservTemp, AbsSalinity, Pressure,
                                     SpecVol, ZMid);
          });
   }
}

/// Define IO fields and metadata for output
void Eos::defineFields() {

   /// Set field names (append Name if not default)
   SpecVolFldName          = "SpecVol";
   SpecVolDisplacedFldName = "SpecVolDisplaced";
   BruntVaisalaFreqFldName = "BruntVaisalaFreq";
   if (Name != "Default") {
      SpecVolFldName.append(Name);
      SpecVolDisplacedFldName.append(Name);
      BruntVaisalaFreqFldName.append(Name);
   }

   /// Create fields for state variables
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLayers";

   /// Create and register the specific volume field
   auto SpecVolField =
       Field::create(SpecVolFldName,                   // Field name
                     "Layer-averaged Specific Volume", // Long Name
                     "m3 kg-1",                        // Units
                     "sea_water_specific_volume",      // CF-ish Name
                     0.0,                              // Min valid value
                     9.99E+30,                         // Max valid value
                     -9.99E+30, // Scalar used for undefined entries
                     NDims,     // Number of dimensions
                     DimNames   // Dimension names
       );
   /// Create and register the displaced specific volume field
   auto SpecVolDisplacedField =
       Field::create(SpecVolDisplacedFldName, // Field name
                     "Specific Volume displaced adiabatically "
                     "to specified layer",                  // long Name
                     "m3 kg-1",                             // Units
                     "sea_water_specific_volume_displaced", // CF-ish Name
                     0.0,                                   // Min valid value
                     9.99E+30,                              // Max valid value
                     -9.99E+30, // Scalar used for undefined entried
                     NDims,     // Number of dimensions
                     DimNames   // Dimension names
       );
   /// Create and register the BruntVaisalaFreq field
   auto BruntVaisalaFreqField =
         Field::create(BruntVaisalaFreqFldName, // Field name
                         "Brunt-Vaisala frequency squared", // Long Name
                         "s-2",                             // Units
                         "sea_water_brunt_vaisala_frequency_squared", // CF-ish Name
                         0.0,                                // Min valid value
                         9.99E+30,                           // Max valid value
                         -9.99E+30, // Scalar used for undefined entries
                         NDims,     // Number of dimensions
                         DimNames   // Dimension names
       );

   // Create a field group for the eos-specific state fields
   EosGroupName = "Eos";
   if (Name != "Default") {
      EosGroupName.append(Name);
   }
   auto EosGroup = FieldGroup::create(EosGroupName);

   // Add fields to the EOS group
   EosGroup->addField(SpecVolDisplacedFldName);
   EosGroup->addField(SpecVolFldName);
   EosGroup->addField(BruntVaisalaFreqFldName);
   
   // Attach Kokkos views to the fields
   SpecVolDisplacedField->attachData<Array2DReal>(SpecVolDisplaced);
   SpecVolField->attachData<Array2DReal>(SpecVol);
   BruntVaisalaFreqField->attachData<Array2DReal>(BruntVaisalaFreq);

} // end defineIOFields

} // namespace OMEGA
