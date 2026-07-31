//===-- ocn/VertMix.cpp - Vertical Mixing Coefficients -----------*- C++
//-*-===//
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
#include "Eos.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "TimeStepper.h"
#include "TriDiagSolvers.h"

namespace OMEGA {

ShearMix::ShearMix(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

ConvectiveMix::ConvectiveMix(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

GradRichardsonNum::GradRichardsonNum(const HorzMesh *Mesh,
                                     const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), GeomZMid(VCoord->GeomZMid),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnCell(Mesh->CellsOnCell), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop),
      MaxLayerEdgeBot(VCoord->MaxLayerEdgeBot), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge) {}

OneTwoOneFilter::OneTwoOneFilter(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

VelVertMixSetupOnEdge::VelVertMixSetupOnEdge(const HorzMesh *Mesh,
                                             const VertCoord *VCoord)
    : Enabled(false), ImplicitBottomDragEnabled(false),
      BottomDragCoeff(0.0_Real), LocRhoSw(RhoSw),
      NVertLayers(VCoord->NVertLayers), CellsOnEdge(Mesh->CellsOnEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerVertMixSetupOnCell::TracerVertMixSetupOnCell(const HorzMesh *Mesh,
                                                   const VertCoord *VCoord)
    : Enabled(false), LocRhoSw(RhoSw), NVertLayers(VCoord->NVertLayers),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

/// Constructor for VertMix
VertMix::VertMix(const std::string &Name, ///< [in] Name for VertMix object
                 const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                 const VertCoord *VCoord  ///< [in] Vertical coordinate
                 )
    : ComputeVertMixConv(VCoord), ComputeVertMixShear(VCoord),
      ComputeGradRichardsonNum(Mesh, VCoord), ComputeOneTwoOneFilter(VCoord),
      Name(Name), Mesh(Mesh), VCoord(VCoord), VelVertMixSetup(Mesh, VCoord),
      TracerVertMixSetup(Mesh, VCoord) {
   VertDiff = Array2DReal("VertDiff", Mesh->NCellsSize, VCoord->NVertLayersP1);
   VertVisc = Array2DReal("VertVisc", Mesh->NCellsSize, VCoord->NVertLayersP1);
   GradRichNum =
       Array2DReal("GradRichNum", Mesh->NCellsSize, VCoord->NVertLayersP1);
   GradRichNumSmoothed = Array2DReal("GradRichNumSmoothed", Mesh->NCellsSize,
                                     VCoord->NVertLayersP1);

   // TODO: Temporary handling of TangentialVelocity
   TangentialVelocity =
       Array2DReal("TangentialVelocity", Mesh->NEdgesSize, VCoord->NVertLayers);

   defineFields();
}

/// Destructor for VertMix
VertMix::~VertMix() {}

/// Instance management
VertMix *VertMix::Instance = nullptr;

/// Get instance of VertMix
VertMix *VertMix::getInstance() { return Instance; }

/// Destroy instance of VertMix
void VertMix::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initializes the VertMix (Vertical Mixing Coefficients) class and its
/// options. It ASSUMES that HorzMesh was initialized and initializes the
/// VertMix class by using the default mesh, reading the config file, and
/// setting parameters for the background, convective, and/or shear mixing
/// routines. Returns 0 on success, or an error code if any required option is
/// missing.
void VertMix::init() {

   if (!Instance) {
      Instance = new VertMix("Default", HorzMesh::getDefault(),
                             VertCoord::getDefault());
   }

   Error Err; // error code

   /// Retrieve default VertMix
   VertMix *DefVertMix = VertMix::getInstance();

   /// Get VertMixConfig group from Omega config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config VertMixConfig("VertMix");
   Err += OmegaConfig->get(VertMixConfig);
   CHECK_ERROR_ABORT(Err, "VertMix::init: VertMix group not found in Config");

   /// Get Background from VertMixConfig
   /// and set associated parameters
   Config BackConfig("Background");
   Err += VertMixConfig.get(BackConfig);
   CHECK_ERROR_ABORT(
       Err, "VertMix::init: Background subgroup not found in VertMixConfig");

   /// Get diffusivity and viscosity parameters
   Err += BackConfig.get("Viscosity", DefVertMix->BackVisc);
   CHECK_ERROR_ABORT(
       Err,
       "VertMix::init: Parameter Background:Viscosity not found in BackConfig");

   Err += BackConfig.get("Diffusivity", DefVertMix->BackDiff);
   CHECK_ERROR_ABORT(Err, "VertMix::init: Parameter Background:Diffusivity not "
                          "found in BackConfig");

   /// Get Convective from VertMixConfig
   Config ConvConfig("Convective");
   Err += VertMixConfig.get(ConvConfig);
   CHECK_ERROR_ABORT(
       Err, "VertMix::init: Convective subgroup not found in VertMixConfig");

   /// Get convective diffusivity and viscosity parameters
   Err += ConvConfig.get("Enable", DefVertMix->ComputeVertMixConv.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "VertMix::init: Parameter Convective:Enable not found in ConvConfig");

   if (!DefVertMix->ComputeVertMixConv.Enabled) {
      LOG_INFO("VertMix::init: Convective mixing is disabled.");
   } else {
      LOG_INFO("VertMix::init: Convective mixing is enabled.");
      Err += ConvConfig.get("Diffusivity",
                            DefVertMix->ComputeVertMixConv.ConvDiff);
      CHECK_ERROR_ABORT(Err, "VertMix::init: Parameter Convective:Diffusivity "
                             "not found in ConvConfig");

      Err += ConvConfig.get("TriggerBVF",
                            DefVertMix->ComputeVertMixConv.ConvTriggerBVF);
      CHECK_ERROR_ABORT(Err, "VertMix::init: Parameter Convective:TriggerBVF "
                             "not found in ConvConfig");
   }

   /// Get Shear from VertMixConfig
   Config ShearConfig("Shear");
   Err += VertMixConfig.get(ShearConfig);
   CHECK_ERROR_ABORT(
       Err, "VertMix::init: Shear subgroup not found in VertMixConfig");

   /// Get shear diffusivity and viscosity parameters
   Err += ShearConfig.get("Enable", DefVertMix->ComputeVertMixShear.Enabled);
   CHECK_ERROR_ABORT(
       Err, "VertMix::init: Parameter Shear:Enable not found in ShearConfig");

   if (!DefVertMix->ComputeVertMixShear.Enabled) {
      LOG_INFO("VertMix::init: Shear mixing is disabled.");
   } else {
      LOG_INFO("VertMix::init: Shear mixing is enabled.");
      Err += ShearConfig.get("BaseShearValue",
                             DefVertMix->ComputeVertMixShear.BaseShearValue);
      CHECK_ERROR_ABORT(Err, "VertMix::init: Parameter Shear:BaseShearValue "
                             "not found in ShearConfig");

      Err += ShearConfig.get("RiCrit",
                             DefVertMix->ComputeVertMixShear.ShearRiCrit);
      CHECK_ERROR_ABORT(
          Err,
          "VertMix::init: Parameter Shear:RiCrit not found in ShearConfig");

      Err += ShearConfig.get("Exponent",
                             DefVertMix->ComputeVertMixShear.ShearExponent);
      CHECK_ERROR_ABORT(
          Err,
          "VertMix::init: Parameter Shear:Exponent not found in ShearConfig");

      Err += ShearConfig.get("RiSmoothLoops",
                             DefVertMix->ComputeVertMixShear.RiSmoothLoops);
      CHECK_ERROR_ABORT(Err, "VertMix::init: Parameter Shear:RiSmoothLoops not "
                             "found in ShearConfig");
   }
} // end init

/// Compute diffusivity and viscosity for all cells/layers (no displacement)
void VertMix::computeVertMix(const Array2DReal &NormalVelocity,
                             const Array2DReal &TangentialVelocity,
                             const Array2DReal &BruntVaisalaFreqSq) {
   OMEGA_SCOPE(LocVertDiff, VertDiff); /// Create a local view for computation
   OMEGA_SCOPE(LocVertVisc, VertVisc); /// Create a local view for computation
   OMEGA_SCOPE(LocGradRichNum,
               GradRichNum); /// Local view for computation
   OMEGA_SCOPE(LocGradRichNumSmoothed,
               GradRichNumSmoothed); /// Local view for computation
   OMEGA_SCOPE(
       LocComputeVertMixConv,
       ComputeVertMixConv); /// Local view for Convective VertMix computation
   OMEGA_SCOPE(
       LocComputeVertMixShear,
       ComputeVertMixShear); /// Local view for Shear VertMix computation
   OMEGA_SCOPE(
       LocComputeGradRichardsonNum,
       ComputeGradRichardsonNum); /// Local view for GradRichNum computation
   OMEGA_SCOPE(LocOneTwoOneFilter,
               ComputeOneTwoOneFilter); /// Local view for 1-2-1 filter
   OMEGA_SCOPE(LocBackDiff, BackDiff);
   OMEGA_SCOPE(LocBackVisc, BackVisc);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   /// First, initialize VertDiff and VertVisc to background values
   parallelForOuter(
       "VertMix-BackAndRich", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const int KStart = chunkStart(KChunk, KMin);
                 const int KLen   = chunkLength(KChunk, KStart, KMax);
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const int K           = KStart + KVec;
                    LocVertDiff(ICell, K) = LocBackDiff;
                    LocVertVisc(ICell, K) = LocBackVisc;
                    LocGradRichNum(ICell, K) =
                        LocComputeGradRichardsonNum.RiInitValue;
                    LocGradRichNumSmoothed(ICell, K) =
                        LocComputeGradRichardsonNum.RiInitValue;
                 }
              });
       });
   /// Second, compute shear mixing if enabled
   if (LocComputeVertMixShear.Enabled) {
      /// Compute Richardson number
      parallelForOuter(
          "VertMix-ComputeRi", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell) + 1;
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeGradRichardsonNum(
                        LocGradRichNum, ICell, KChunk, NormalVelocity,
                        TangentialVelocity, BruntVaisalaFreqSq);
                 });

             teamBarrier(Team);

             // Fill Richardson number at vertical boundaries using the
             // closest valid value. This is equivalent to doing one-sided
             // differencing at the boundary.
             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    LocGradRichNum(ICell, MinLayerCell(ICell)) =
                        LocGradRichNum(ICell, KMin);
                    LocGradRichNum(ICell, MaxLayerCell(ICell) + 1) =
                        LocGradRichNum(ICell, KMax);
                 });
          });
      /// Smooth Richardson number with 1-2-1 filter the number of times
      /// specified by RiSmoothLoops
      deepCopy(LocGradRichNumSmoothed, LocGradRichNum);
      for (int SmoothLoop = 0;
           SmoothLoop < LocComputeVertMixShear.RiSmoothLoops; ++SmoothLoop) {
         parallelForOuter(
             "VertMix-RiSmooth", {Mesh->NCellsAll},
             KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
                const int KMin   = MinLayerCell(ICell);
                const int KMax   = MaxLayerCell(ICell);
                const int KRange = vertRangeChunked(KMin, KMax);
                parallelForInner(
                    Team, KRange, INNER_LAMBDA(int KChunk) {
                       LocOneTwoOneFilter(LocGradRichNumSmoothed, ICell, KChunk,
                                          LocGradRichNumSmoothed);
                    });
             });
      }
      /// Compute shear mixing using smoothed Richardson number
      parallelForOuter(
          "VertMix-Shear", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell) + 1;
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeVertMixShear(LocVertDiff, LocVertVisc, ICell,
                                           KChunk, LocGradRichNumSmoothed);
                 });
          });
   }
   /// Third, compute convective mixing if enabled
   if (LocComputeVertMixConv.Enabled) {
      parallelForOuter(
          "VertMix-Conv", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell) + 1;
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeVertMixConv(LocVertDiff, LocVertVisc, ICell,
                                          KChunk, BruntVaisalaFreqSq);
                 });
          });
   }
   /// Finally, zero viscosity/diffusivity at surface and bottom boundaries
   parallelFor(
       "VertMix-Boundaries", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          const int KMin           = MinLayerCell(ICell);
          const int KMax           = MaxLayerCell(ICell) + 1;
          LocVertDiff(ICell, KMin) = 0.0_Real;
          LocVertVisc(ICell, KMin) = 0.0_Real;
          LocVertDiff(ICell, KMax) = 0.0_Real;
          LocVertVisc(ICell, KMax) = 0.0_Real;
       });
}

/// Define IO fields and metadata for output
void VertMix::defineFields() {

   /// Set field names (append Name if not default)
   VertDiffFldName            = "VertDiff";
   VertViscFldName            = "VertVisc";
   GradRichNumFldName         = "GradRichNum";
   GradRichNumSmoothedFldName = "GradRichNumSmoothed";
   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
      GradRichNumFldName.append(Name);
      GradRichNumSmoothedFldName.append(Name);
   }

   /// Create fields for state variables
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLayersP1";

   /// Create and register the Diffusivity field
   auto VertDiffField =
       Field::create(VertDiffFldName, // Field name
                     "Vertical diffusivity at center of"
                     " cell and top of layer",         // Long Name
                     "m2 s-1",                         // Units
                     "vertical_diffusivity",           // CF-ish Name
                     0.0,                              // Min valid value
                     std::numeric_limits<Real>::max(), // Max valid value
                     NDims,                            // Number of dimensions
                     DimNames                          // Dimension names
       );
   /// Create and register the VertVisc field
   auto VertViscField =
       Field::create(VertViscFldName, // Field name
                     "Vertical viscosity at center of"
                     " cell and top of layer",         // Long Name
                     "m2 s-1",                         // Units
                     "vertical_viscosity",             // CF-ish Name
                     0.0,                              // Min valid value
                     std::numeric_limits<Real>::max(), // Max valid value
                     NDims,                            // Number of dimensions
                     DimNames                          // Dimension names
       );
   /// Create and register the GradRichNum field
   auto GradRichNumField =
       Field::create(GradRichNumFldName,                     // Field name
                     "Gradient Richardson number",           // Long Name
                     "dimensionless",                        // Units
                     "sea_water_gradient_richardson_number", // CF-ish Name
                     std::numeric_limits<Real>::min(),       // Min valid value
                     std::numeric_limits<Real>::max(),       // Max valid value
                     NDims,   // Number of dimensions
                     DimNames // Dimension names
       );
   /// Create and register the GradRichNumSmoothed field
   auto GradRichNumSmoothedField = Field::create(
       GradRichNumSmoothedFldName,                      // Field name
       "Smoothed Gradient Richardson number",           // Long Name
       "dimensionless",                                 // Units
       "sea_water_gradient_richardson_number_smoothed", // CF-ish Name
       std::numeric_limits<Real>::min(),                // Min valid value
       std::numeric_limits<Real>::max(),                // Max valid value
       NDims,                                           // Number of dimensions
       DimNames                                         // Dimension names
   );

   // Create a field group for the vertmix-specific state fields
   VertMixGroupName = "VertMix";
   if (Name != "Default") {
      VertMixGroupName.append(Name);
   }
   auto VertMixGroup = FieldGroup::create(VertMixGroupName);

   // Add fields to the VertMix group
   VertMixGroup->addField(VertDiffFldName);
   VertMixGroup->addField(VertViscFldName);
   VertMixGroup->addField(GradRichNumFldName);
   VertMixGroup->addField(GradRichNumSmoothedFldName);

   // Attach Kokkos views to the fields
   VertDiffField->attachData<Array2DReal>(VertDiff);
   VertViscField->attachData<Array2DReal>(VertVisc);
   GradRichNumField->attachData<Array2DReal>(GradRichNum);
   GradRichNumSmoothedField->attachData<Array2DReal>(GradRichNumSmoothed);

} // end defineIOFields

// Apply implicit velocity vertical mixing
void VertMix::applyVelVertMixImplicit(
    OceanState *State,              ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel                ///< [in] Time level
) {

   OMEGA_SCOPE(LocNEdgesAll, Mesh->NEdgesAll);
   OMEGA_SCOPE(LocVelVertMixSetup, VelVertMixSetup);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   const Array2DReal &NormalVelEdge = State->getNormalVelocity(VelTimeLevel);
   const Array2DReal &PseudoThickCell =
       State->getPseudoThickness(ThickTimeLevel);

   // Compute velocity vertical mixing
   if (LocVelVertMixSetup.Enabled) {
      Pacer::start("Tend:velocityVertMix", 1);

      Eos *EosInstance         = Eos::getInstance();
      VertMix *VertMixInstance = VertMix::getInstance();

      // Obtain TimeStep
      const auto *DefTimeStepper  = TimeStepper::getDefault();
      const TimeInterval TimeStep = DefTimeStepper->getTimeStep();
      R8 DT;
      TimeStep.get(DT, TimeUnits::Seconds);

      const auto &SpecVol           = EosInstance->SpecVol;
      const auto &VertVisc          = VertMixInstance->VertVisc;
      const auto &KineticEnergyCell = AuxState->KineticAux.KineticEnergyCell;

      const int NVertLayers  = VCoord->NVertLayers;
      const int LocVecLength = VecLength;
      auto LConfig =
          TriDiagSolver::makeLaunchConfig(Mesh->NEdgesAll, NVertLayers);

      parallelForOuter(
          LConfig, KOKKOS_LAMBDA(int, const TeamMember &Team) {
             const int IStart = Team.league_rank() * LocVecLength;
             const int ILen   = Kokkos::max(
                 0, Kokkos::min(LocVecLength, LocNEdgesAll - IStart));

             TriDiagDiffScratch Scratch(Team, NVertLayers);

             // Construct a tri-diag diffusion matrix and RHS
             parallelForInner(Team, NVertLayers, [=](int K) {
                for (int IVec = 0; IVec < LocVecLength; ++IVec) {
                   const int IEdge = IStart + IVec;

                   if (IEdge >= LocNEdgesAll) {
                      // Fill values
                      Scratch.G(K, IVec) = 0._Real;
                      Scratch.H(K, IVec) = 1._Real;
                      Scratch.X(K, IVec) = 0._Real;
                      continue;
                   }

                   const int KMin = MinLayerEdgeBot(IEdge);
                   const int KMax = MaxLayerEdgeTop(IEdge);

                   if (K < KMin || K > KMax) {
                      // Fill values
                      Scratch.G(K, IVec) = 0._Real;
                      Scratch.H(K, IVec) = 1._Real;
                      Scratch.X(K, IVec) = 0._Real;
                      continue;
                   }

                   Real G, H, X;
                   LocVelVertMixSetup(IEdge, K, KMin, KMax, DT, SpecVol,
                                      KineticEnergyCell, PseudoThickCell,
                                      VertVisc, NormalVelEdge, G, H, X);

                   Scratch.G(K, IVec) = G;
                   Scratch.H(K, IVec) = H;
                   Scratch.X(K, IVec) = X;
                }
             });

             // Solve the tri-diag diffusion system
             Team.team_barrier();
             TriDiagDiffSolver::solve(Team, Scratch);
             Team.team_barrier();

             // Store the solution vector X
             parallelForInner(Team, NVertLayers, [=](int K) {
                for (int IVec = 0; IVec < ILen; ++IVec) {
                   const int IEdge = IStart + IVec;

                   if (K >= MinLayerEdgeBot(IEdge) &&
                       K <= MaxLayerEdgeTop(IEdge)) {
                      NormalVelEdge(IEdge, K) = Scratch.X(K, IVec);
                   }
                }
             });
          });

      Pacer::stop("Tend:velocityVertMix", 1);
   }
} // applyVelVertMixImplicit

// Apply implicit tracer vertical mixing
void VertMix::applyTracerVertMixImplicit(
    OceanState *State,              ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    Array3DReal &TracerArray,       ///< [in] Tracer array
    int NTracers,                   ///< [in] Number of tracers
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel                ///< [in] Time level
) {

   OMEGA_SCOPE(LocNCellsAll, Mesh->NCellsAll);
   OMEGA_SCOPE(LocTracerVertMixSetup, TracerVertMixSetup);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   const Array2DReal &PseudoThickCell =
       State->getPseudoThickness(ThickTimeLevel);

   if (LocTracerVertMixSetup.Enabled) {
      Pacer::start("Tend:tracerVertMix", 1);

      Eos *EosInstance         = Eos::getInstance();
      VertMix *VertMixInstance = VertMix::getInstance();

      // Obtain TimeStep
      const auto *DefTimeStepper  = TimeStepper::getDefault();
      const TimeInterval TimeStep = DefTimeStepper->getTimeStep();
      R8 DT;
      TimeStep.get(DT, TimeUnits::Seconds);

      const auto &SpecVol  = EosInstance->SpecVol;
      const auto &VertDiff = VertMixInstance->VertDiff;

      const int NVertLayers = VCoord->NVertLayers;
      auto LConfig =
          TriDiagSolver::makeLaunchConfig(Mesh->NCellsAll, NVertLayers);
      const int LocVecLength = VecLength;

      for (int L = 0; L < NTracers; ++L) {
         parallelForOuter(
             LConfig, KOKKOS_LAMBDA(int, const TeamMember &Team) {
                const int IStart = Team.league_rank() * LocVecLength;
                const int ILen   = Kokkos::max(
                    0, Kokkos::min(LocVecLength, LocNCellsAll - IStart));

                TriDiagDiffScratch Scratch(Team, NVertLayers);

                // Construct a tri-diag diffusion matrix and RHS
                parallelForInner(Team, NVertLayers, [=](int K) {
                   for (int IVec = 0; IVec < LocVecLength; ++IVec) {
                      const int ICell = IStart + IVec;

                      if (ICell >= LocNCellsAll) {
                         // Fill values
                         Scratch.G(K, IVec) = 0._Real;
                         Scratch.H(K, IVec) = 1._Real;
                         Scratch.X(K, IVec) = 0._Real;
                         continue;
                      }

                      const int KMin = MinLayerCell(ICell);
                      const int KMax = MaxLayerCell(ICell);

                      if (K < KMin || K > KMax) {
                         // Fill values
                         Scratch.G(K, IVec) = 0._Real;
                         Scratch.H(K, IVec) = 1._Real;
                         Scratch.X(K, IVec) = 0._Real;
                         continue;
                      }

                      Real G, H, X;
                      LocTracerVertMixSetup(L, ICell, K, KMin, KMax, DT,
                                            SpecVol, PseudoThickCell, VertDiff,
                                            TracerArray, G, H, X);
                      Scratch.G(K, IVec) = G;
                      Scratch.H(K, IVec) = H;
                      Scratch.X(K, IVec) = X;
                   }
                });

                // Solve the tri-diag diffusion system
                Team.team_barrier();
                TriDiagDiffSolver::solve(Team, Scratch);
                Team.team_barrier();

                // Store the solution vector X
                parallelForInner(Team, NVertLayers, [=](int K) {
                   for (int IVec = 0; IVec < ILen; ++IVec) {
                      const int ICell = IStart + IVec;

                      if (K >= MinLayerCell(ICell) &&
                          K <= MaxLayerCell(ICell)) {
                         TracerArray(L, ICell, K) = Scratch.X(K, IVec);
                      }
                   }
                });
             });

      } // for L

      Pacer::stop("Tend:tracerVertMix", 1);
   }

} // applyTracerVertMixImplicit

/// Apply implicit vertical mixing to velocities and tracers
void VertMix::VertMixImplicit(OceanState *State, AuxiliaryState *AuxState,
                              Array3DReal &TracerArray, int NTracers,
                              int TimeLevel) {

   // get NormalVelocity
   Array2DReal NormalVelEdge = State->getNormalVelocity(TimeLevel);

   // get temperature and salinity
   I4 ConservTempIdx;
   I4 AbsSalinityIdx;
   Tracers::getIndex(ConservTempIdx, "Temperature");
   Tracers::getIndex(AbsSalinityIdx, "Salinity");

   const auto ConservTemp =
       Kokkos::subview(TracerArray, ConservTempIdx, Kokkos::ALL, Kokkos::ALL);
   const auto AbsSalinity =
       Kokkos::subview(TracerArray, AbsSalinityIdx, Kokkos::ALL, Kokkos::ALL);

   // get an instance of equation of state
   Eos *EqState = Eos::getInstance();

   // TODO: Temporary handling of computation of tangential velocity
   // Compute tangential velocity
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(LocTangentialVelocity, TangentialVelocity);

   TangentialReconOnEdge TanReconEdge(Mesh);

   parallelForOuter(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeBot(IEdge);
          const int KMax   = MaxLayerEdgeTop(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 TanReconEdge(LocTangentialVelocity, IEdge, KChunk,
                              NormalVelEdge);
              });
       });

   // Refresh kinetic energy for the implicit bottom drag
   // from the pre-vmix velocity
   if (VelVertMixSetup.ImplicitBottomDragEnabled) {
      OMEGA_SCOPE(LocKineticAux, AuxState->KineticAux);
      OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
      OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

      Pacer::start("VertMix:computeKineticAuxForBottomDrag", 2);
      parallelForOuter(
          "computeKineticAuxForBottomDrag", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocKineticAux.computeVarsOnCell(ICell, KChunk,
                                                    NormalVelEdge);
                 });
          });
      Pacer::stop("VertMix:computeKineticAuxForBottomDrag", 2);
   }

   // Update Pressure, SpecVol
   AuxState->computeMomVertAux(State, TracerArray, TimeLevel, TimeLevel);

   // Compute Brunt-Vaisala frequency squared
   EqState->computeBruntVaisalaFreqSq(
       ConservTemp, AbsSalinity, VCoord->PressureInterface, EqState->SpecVol);

   // Compute vertical mixing coefficients
   computeVertMix(NormalVelEdge, LocTangentialVelocity,
                  EqState->BruntVaisalaFreqSq);

   // Apply implicit mixing to velocities
   applyVelVertMixImplicit(State, AuxState, TimeLevel, TimeLevel);

   // Apply implicit mixing to tracers
   applyTracerVertMixImplicit(State, AuxState, TracerArray, NTracers, TimeLevel,
                              TimeLevel);

} // VertMixImplicit

} // namespace OMEGA
