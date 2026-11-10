#include "SubmesoEddies.h"
#include "Field.h"
#include "GlobalConstants.h"

namespace OMEGA {

static KOKKOS_FUNCTION Real shapeFunction(Real Z, Real H) {
   const Real Tmp = (2 * Z) / H + 1;
   return Kokkos::max(0._Real, (1 - Tmp * Tmp) * (1 + 5 * Tmp * Tmp / 21));
}

/// Instance management
SubmesoEddies *SubmesoEddies::Instance = nullptr;

/// Get instance of SubmesoEddies
SubmesoEddies *SubmesoEddies::getInstance() { return Instance; }

/// Destroy instance of SubmesoEddies
void SubmesoEddies::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initializes the SubmesoEddies class and its options.
/// It assumes that HorzMesh and VertCoord were initialized and
/// initializes the SubmesoEddies class by using the default mesh and vertical
/// coordinate, reading the config file, and setting the parametrization
/// parameters.
void SubmesoEddies::init() {

   if (!Instance) {
      Instance =
          new SubmesoEddies(HorzMesh::getDefault(), VertCoord::getDefault());
   }

   Error Err; // error code

   /// Retrieve default eos
   SubmesoEddies *SubEddies = SubmesoEddies::getInstance();

   /// Get Submeso group from Omega config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config SubmesoConfig("Submeso");
   Err += OmegaConfig->get(SubmesoConfig);
   CHECK_ERROR_ABORT(Err,
                     "SubmesoEddies::init: Submeso group not found in Config");

   Err += SubmesoConfig.get("Enable", SubEddies->Enable);
   CHECK_ERROR_ABORT(Err,
                     "SubmesoEddies::init: Enable not found in SubmesoConfig");

   Err += SubmesoConfig.get("Tau", SubEddies->Tau);
   CHECK_ERROR_ABORT(Err,
                     "SubmesoEddies::init: Tau not found in SubmesoConfig");

   Err += SubmesoConfig.get("Ce", SubEddies->Ce);
   CHECK_ERROR_ABORT(Err, "SubmesoEddies::init: Ce not found in SubmesoConfig");

   Err += SubmesoConfig.get("LfMin", SubEddies->LfMin);
   CHECK_ERROR_ABORT(Err,
                     "SubmesoEddies::init: LfMin not found in SubmesoConfig");

   Err += SubmesoConfig.get("DsMax", SubEddies->DsMax);
   CHECK_ERROR_ABORT(Err,
                     "SubmesoEddies::init: DsMax not found in SubmesoConfig");

   // Precompute time scale
   SubEddies->computeTimeScale();

} // end init

SubmesoEddies::SubmesoEddies(const HorzMesh *Mesh, const VertCoord *VCoord)
    : Mesh(Mesh), VCoord(VCoord), TimeScale("TimeScale", Mesh->NEdgesSize),
      GradBuoyEdgeInterface("GradBuoyEdgeInterface", Mesh->NEdgesSize,
                            VCoord->NVertLayersP1),
      EddyVelocity("EddyVelocity", Mesh->NEdgesSize, VCoord->NVertLayers) {

   // define fields for IO
   defineFields();
}

void SubmesoEddies::defineFields() {

   // Create a group for the submesoscale eddy parametrization fields
   auto SubmesoGroup = FieldGroup::create("Submeso");

   // Create and add buoyancy gradient field
   {
      int NDims = 2;
      std::vector<std::string> DimNames(NDims);
      DimNames[0] = "NEdges";
      DimNames[1] = "NVertLayersP1";

      auto BuoyancyGradientInterfaceField =
          Field::create(GradBuoyEdgeInterface.label(),       // Field name
                        "Buoyancy Gradient",                 // Long Name
                        "1/s^2",                             // Units
                        "",                                  // CF-ish Name
                        std::numeric_limits<Real>::lowest(), // Min valid value
                        std::numeric_limits<Real>::max(),    // Max valid value
                        NDims,   // Number of dimensions
                        DimNames // Dimension names
          );

      BuoyancyGradientInterfaceField->attachData<Array2DReal>(
          GradBuoyEdgeInterface);

      SubmesoGroup->addField(GradBuoyEdgeInterface.label());
   }

   // Create and add eddy velocity field
   {
      int NDims = 2;
      std::vector<std::string> DimNames(NDims);
      DimNames[0] = "NEdges";
      DimNames[1] = "NVertLayers";

      auto EddyVelocityField =
          Field::create(EddyVelocity.label(),                // Field name
                        "Eddy Velocity",                     // Long Name
                        "m/s",                               // Units
                        "",                                  // CF-ish Name
                        std::numeric_limits<Real>::lowest(), // Min valid value
                        std::numeric_limits<Real>::max(),    // Max valid value
                        NDims,   // Number of dimensions
                        DimNames // Dimension names
          );

      EddyVelocityField->attachData<Array2DReal>(EddyVelocity);

      SubmesoGroup->addField(EddyVelocity.label());
   }
}

void SubmesoEddies::computeTimeScale() {
   OMEGA_SCOPE(TimeScale, this->TimeScale);
   OMEGA_SCOPE(Tau, this->Tau);

   const auto &FEdge = Mesh->FEdge;

   parallelFor(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          TimeScale(IEdge) =
              Kokkos::sqrt(FEdge(IEdge) * FEdge(IEdge) + 1._Real / (Tau * Tau));
       });
}

void SubmesoEddies::computeBuoyGrad(const Array2DReal &SpecVol,
                                    const Array2DReal &MeanPseudoThickEdge,
                                    const Array2DReal &GeomZMid,
                                    const Array2DReal &BruntVaisalaFreqSq) {
   OMEGA_SCOPE(GradBuoyEdgeInterface, this->GradBuoyEdgeInterface);

   const auto &DcEdge      = Mesh->DcEdge;
   const auto &CellsOnEdge = Mesh->CellsOnEdge;

   const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   const auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   const auto NVertLayers      = VCoord->NVertLayers;
   const auto NVertLayersP1    = VCoord->NVertLayersP1;

   // We use six scratch arrays in the loop below
   // Three for layer data and three for interface data
   const auto ScratchMem =
       TeamScratch<Real>(3 * NVertLayers + 3 * NVertLayersP1);
   parallelForOuter(
       LaunchConfig({Mesh->NEdgesAll}, ScratchMem),
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          int MinLyrEdgeBot = MinLayerEdgeBot(IEdge);
          int MaxLyrEdgeTop = MaxLayerEdgeTop(IEdge);

          ScratchArray1DReal GradBuoyEdge(teamScratch(Team), NVertLayers);
          ScratchArray1DReal SpecVolEdge(teamScratch(Team), NVertLayers);
          ScratchArray1DReal GradGeomZMidEdge(teamScratch(Team), NVertLayers);

          ScratchArray1DReal BVFSqEdge(teamScratch(Team), NVertLayersP1);
          ScratchArray1DReal GradGeomZMidEdgeInterface(teamScratch(Team),
                                                       NVertLayersP1);
          ScratchArray1DReal SpecVolEdgeInterface(teamScratch(Team),
                                                  NVertLayersP1);

          // Horizontal interpolations and gradients
          parallelForInner(
              Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                 const int JCell0 = CellsOnEdge(IEdge, 0);
                 const int JCell1 = CellsOnEdge(IEdge, 1);

                 const Real InvDcEdge = 1._Real / DcEdge(IEdge);

                 // interpolate Brunt-Vaisala freq to edges
                 BVFSqEdge(K) = 0.5_Real * (BruntVaisalaFreqSq(JCell1, K) +
                                            BruntVaisalaFreqSq(JCell0, K));

                 // interpolate SpecVol to Edges
                 SpecVolEdge(K) =
                     0.5_Real * (SpecVol(JCell1, K) + SpecVol(JCell0, K));

                 // Compute grad of GeomZMid
                 GradGeomZMidEdge(K) =
                     InvDcEdge * (GeomZMid(JCell1, K) - GeomZMid(JCell0, K));

                 // Compute grad of buoyancy
                 GradBuoyEdge(K) = -Gravity * RhoSw * InvDcEdge *
                                   (SpecVol(JCell1, K) - SpecVol(JCell0, K));
              });

          // Brunt-Vaisala freq is located on interfaces
          // Handle one remaining interpolation to edges at the ocean bottom
          Kokkos::single(
              PerTeam(Team), INNER_LAMBDA() {
                 const int K      = MaxLyrEdgeTop + 1;
                 const int JCell0 = CellsOnEdge(IEdge, 0);
                 const int JCell1 = CellsOnEdge(IEdge, 1);
                 BVFSqEdge(K)     = 0.5_Real * (BruntVaisalaFreqSq(JCell1, K) +
                                            BruntVaisalaFreqSq(JCell0, K));
              });

          teamBarrier(Team);

          // Vertical interpolations

          // This interpolation can only be carried out on non-boundary edges
          if (MaxLyrEdgeTop >= MinLyrEdgeBot) {

             parallelForInner(
                 Team, Range{MinLyrEdgeBot + 1, MaxLyrEdgeTop},
                 INNER_LAMBDA(int K) {
                    const Real PseudoThickKm1 =
                        MeanPseudoThickEdge(IEdge, K - 1);
                    const Real PseudoThickK = MeanPseudoThickEdge(IEdge, K);

                    const Real CoeffKm1 =
                        PseudoThickKm1 / (PseudoThickKm1 + PseudoThickK);
                    const Real CoeffK =
                        PseudoThickK / (PseudoThickKm1 + PseudoThickK);

                    // interpolate GeomZMid gradient to interfaces
                    GradGeomZMidEdgeInterface(K) =
                        CoeffKm1 * GradGeomZMidEdge(K - 1) +
                        CoeffK * GradGeomZMidEdge(K);

                    // interpolate SpecVol to interfaces
                    SpecVolEdgeInterface(K) =
                        CoeffKm1 * SpecVolEdge(K - 1) + CoeffK * SpecVolEdge(K);

                    // interpolate GradBuoyEdge to interfaces
                    GradBuoyEdgeInterface(IEdge, K) =
                        CoeffKm1 * GradBuoyEdge(K - 1) +
                        CoeffK * GradBuoyEdge(K);
                 });

             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    SpecVolEdgeInterface(MinLyrEdgeBot) =
                        SpecVolEdge(MinLyrEdgeBot);
                    SpecVolEdgeInterface(MaxLyrEdgeTop + 1) =
                        SpecVolEdge(MaxLyrEdgeTop);

                    GradGeomZMidEdgeInterface(MinLyrEdgeBot) =
                        GradGeomZMidEdge(MinLyrEdgeBot);
                    GradGeomZMidEdgeInterface(MaxLyrEdgeTop + 1) =
                        GradGeomZMidEdge(MaxLyrEdgeTop);

                    GradBuoyEdgeInterface(IEdge, MinLyrEdgeBot) =
                        GradBuoyEdge(MinLyrEdgeBot);
                    GradBuoyEdgeInterface(IEdge, MaxLyrEdgeTop + 1) =
                        GradBuoyEdge(MaxLyrEdgeTop);
                 });

             teamBarrier(Team);

             parallelForInner(
                 Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop + 1},
                 INNER_LAMBDA(int K) {
                    GradBuoyEdgeInterface(IEdge, K) +=
                        GradGeomZMidEdgeInterface(K) * RhoSw *
                        SpecVolEdgeInterface(K) * BVFSqEdge(K);
                 });
          }
       });
}

void SubmesoEddies::computeEddyVelocity(
    const Array1DReal &DenMixLayerDepth, const Array1DI4 &DenMixLayerIndex,
    const Array2DReal &BruntVaisalaFreqSq,
    const Array2DReal &MeanPseudoThickEdge) {

   OMEGA_SCOPE(GradBuoyEdgeInterface, this->GradBuoyEdgeInterface);
   OMEGA_SCOPE(TimeScale, this->TimeScale);
   OMEGA_SCOPE(LfMin, this->LfMin);
   OMEGA_SCOPE(DsMax, this->DsMax);
   OMEGA_SCOPE(Ce, this->Ce);
   OMEGA_SCOPE(EddyVelocity, this->EddyVelocity);

   const auto &DcEdge      = Mesh->DcEdge;
   const auto &CellsOnEdge = Mesh->CellsOnEdge;

   const auto &GeomZInterface = VCoord->GeomZInterface;

   const auto &MinLayerCell    = VCoord->MinLayerCell;
   const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   const auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   const auto NVertLayersP1    = VCoord->NVertLayersP1;

   // Replace with global constant when added
   const Real Tiny = 1e-12_Real;

   parallelForOuter(
       LaunchConfig({Mesh->NEdgesAll}, TeamScratch<Real>(NVertLayersP1)),
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int MinLyrEdgeBot = MinLayerEdgeBot(IEdge);
          const int MaxLyrEdgeTop = MaxLayerEdgeTop(IEdge);

          if (MaxLyrEdgeTop >= MinLyrEdgeBot) {

             const int JCell0 = CellsOnEdge(IEdge, 0);
             const int JCell1 = CellsOnEdge(IEdge, 1);

             const int MinLayerCell0 = MinLayerCell(JCell0);
             const int MinLayerCell1 = MinLayerCell(JCell1);

             Real PseudoThickML;
             Real GradBuoyML;
             Real BVFreqML;

             const I4 IndexMLEdge = Kokkos::min(DenMixLayerIndex(JCell0),
                                                DenMixLayerIndex(JCell1));

             // compute mixed layer averages of buoyancy gradient and
             // Brunt-Vaisala frequency
             parallelReduceInner(
                 Team, Range{MinLyrEdgeBot, IndexMLEdge},
                 INNER_LAMBDA(int K, Real &AccumThick, Real &AccumGradBuoy,
                              Real &AccumBVFreq) {
                    const Real PseudoThickKm1 =
                        ((K - 1) >= MinLyrEdgeBot
                             ? MeanPseudoThickEdge(IEdge, K - 1)
                             : 0);
                    const Real PseudoThickK = MeanPseudoThickEdge(IEdge, K);
                    const Real PseudoThickAvg =
                        0.5_Real * (PseudoThickKm1 + PseudoThickK);

                    const Real BVFSq0   = BruntVaisalaFreqSq(JCell0, K);
                    const Real BVFSq1   = BruntVaisalaFreqSq(JCell1, K);
                    const Real GradBuoy = GradBuoyEdgeInterface(IEdge, K);

                    const Real BVFEdge =
                        Kokkos::sqrt(0.5_Real * (Kokkos::max(Tiny, BVFSq0) +
                                                 Kokkos::max(Tiny, BVFSq1)));

                    AccumThick += PseudoThickAvg;
                    AccumGradBuoy += PseudoThickAvg * GradBuoy;
                    AccumBVFreq += PseudoThickAvg * BVFEdge;
                 },
                 PseudoThickML, GradBuoyML, BVFreqML);

             GradBuoyML /= PseudoThickML;
             BVFreqML /= PseudoThickML;

             // compute stream function
             ScratchArray1DReal StreamFunction(teamScratch(Team),
                                               NVertLayersP1);
             parallelForInner(
                 Team, NVertLayersP1,
                 INNER_LAMBDA(int K) { StreamFunction(K) = 0; });

             const Real MLDepthEdge = Kokkos::min(DenMixLayerDepth(JCell0),
                                                  DenMixLayerDepth(JCell1));

             const Real TScale = TimeScale(IEdge);
             const Real Ds     = Kokkos::min(DcEdge(IEdge), DsMax);

             const Real Lf1 =
                 Kokkos::abs(GradBuoyML) * MLDepthEdge / (TScale * TScale);
             const Real Lf2 = BVFreqML * MLDepthEdge / TScale;
             const Real Lf  = Kokkos::max(LfMin, Kokkos::max(Lf1, Lf2));

             const Real Factor =
                 Ce * Ds / Lf * MLDepthEdge * MLDepthEdge * GradBuoyML / TScale;

             parallelForInner(
                 Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop + 1},
                 INNER_LAMBDA(int K) {
                    const Real ZEdge =
                        0.5_Real * (GeomZInterface(JCell0, K) -
                                    GeomZInterface(JCell0, MinLayerCell0) +
                                    GeomZInterface(JCell1, K) -
                                    GeomZInterface(JCell1, MinLayerCell1));

                    const Real Mu = shapeFunction(ZEdge, MLDepthEdge);

                    StreamFunction(K) = Factor * Mu;
                 });

             teamBarrier(Team);

             // compute eddy velocity
             parallelForInner(
                 Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop},
                 INNER_LAMBDA(int K) {
                    const Real DZ = 0.5_Real * (GeomZInterface(JCell0, K) -
                                                GeomZInterface(JCell0, K + 1) +
                                                GeomZInterface(JCell1, K) -
                                                GeomZInterface(JCell1, K + 1));
                    EddyVelocity(IEdge, K) =
                        -(StreamFunction(K) - StreamFunction(K + 1)) / DZ;
                 });
          }
       });
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
