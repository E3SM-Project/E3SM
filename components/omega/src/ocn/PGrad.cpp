//===-- ocn/PGrad.cpp - Pressure Gradient Term -----------------*- C++ -*-===//
//
// Implements the PGrad manager and two discretizations: Centered and
// FiniteVolume.
//
//===----------------------------------------------------------------------===//

#include "PGrad.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

namespace OMEGA {

PressureGrad *PressureGrad::DefaultPGrad = nullptr;
std::map<std::string, std::unique_ptr<PressureGrad>> PressureGrad::AllPGrads;

//------------------------------------------------------------------------------
// Initialize the PressureGrad. Assumes that HorzMesh and VertCoord have already
// been initialized.
void PressureGrad::init() {

   // Retrieve default mesh and vertical coordinate
   HorzMesh *DefMesh = HorzMesh::getDefault();
   OMEGA_REQUIRE(DefMesh,
                 "Null default HorzMesh pointer in PressureGrad::init");
   VertCoord *DefVCoord = VertCoord::getDefault();
   OMEGA_REQUIRE(DefVCoord,
                 "Null default VertCoord pointer in PressureGrad::init");

   // Retrieve omega config
   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in PressureGrad::init");

   // Create the default PressureGrad and set pointer to it
   PressureGrad::DefaultPGrad =
       PressureGrad::create("Default", DefMesh, DefVCoord, OmegaConfig);

} // end init

//------------------------------------------------------------------------------
// Create a new PressureGrad object and add it to the map
PressureGrad *
PressureGrad::create(const std::string &Name, /// [in] Name for PressureGrad
                     const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                     const VertCoord *VCoord, ///< [in] Vertical coordinate
                     Config *Options) {       ///< [in] Configuration options

   OMEGA_REQUIRE(Mesh,
                 "Null HorzMesh pointer in PressureGrad::create with Name = {}",
                 Name);
   OMEGA_REQUIRE(
       VCoord, "Null VertCoord pointer in PressureGrad::create with Name = {}",
       Name);
   OMEGA_REQUIRE(Options,
                 "Null Config pointer in PressureGrad::create with Name = {}",
                 Name);

   // Check to see if a PressureGrad of the same name already exists and
   // if so, exit with an error
   if (AllPGrads.find(Name) != AllPGrads.end()) {
      LOG_ERROR("Attempted to create a PressureGrad with name {} but a "
                "PressureGrad of that name already exists",
                Name);
      return nullptr;
   }

   // create a new PressureGrad on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewPGrad = new PressureGrad(Mesh, VCoord, Options);
   AllPGrads.emplace(Name, NewPGrad);

   return NewPGrad;

} // end create

//------------------------------------------------------------------------------
// Get the default pressure gradient instance
PressureGrad *PressureGrad::getDefault() {

   return DefaultPGrad;

} // end get default

//------------------------------------------------------------------------------
// Constructor for PressureGrad
PressureGrad::PressureGrad(
    const HorzMesh *Mesh,    ///< [in] Horizontal mesh
    const VertCoord *VCoord, ///< [in] Vertical coordinate
    Config *Options)         ///< [in] Configuration options
    : MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop), CenteredPGrad(Mesh, VCoord),
      FiniteVolumePGrad(Mesh, VCoord), BarotropicPGrad(Mesh, VCoord) {

   // store mesh sizes
   NEdgesAll     = Mesh->NEdgesAll;
   NEdgesOwned   = Mesh->NEdgesOwned;
   NVertLayers   = VCoord->NVertLayers;
   NVertLayersP1 = NVertLayers + 1;

   // Read config options for PressureGrad type
   // and enable the appropriate functor
   Config PGradConfig("PressureGrad");
   Error Err;
   Err += Options->get(PGradConfig);
   CHECK_ERROR_ABORT(Err,
                     "PressureGrad: PressureGrad group not found in Config");
   std::string PGradTypeStr;
   Err += PGradConfig.get("PressureGradType", PGradTypeStr);

   if (PGradTypeStr == "centered" || PGradTypeStr == "Centered") {
      PressureGradChoice          = PressureGradType::Centered;
      this->CenteredPGrad.Enabled = true;
   } else if (PGradTypeStr == "finiteVolume" ||
              PGradTypeStr == "FiniteVolume") {
      PressureGradChoice              = PressureGradType::FiniteVolume;
      this->FiniteVolumePGrad.Enabled = true;
   } else {
      // Aborting rather than falling back to the centered scheme: a silent
      // fallback turns a typo, or a configuration naming a scheme that no
      // longer exists, into a run that looks like a passing Centered run
      ABORT_ERROR("PressureGrad: unknown PressureGradType '{}'; valid values "
                  "are 'Centered' and 'FiniteVolume'",
                  PGradTypeStr);
   }

   // Read the FiniteVolume sub-options. All three are optional so that a
   // configuration written before they existed still parses; where a key is
   // absent the default is used.
   I4 HorzOrder = FiniteVolumePGrad.HorzOrder;
   if (PGradConfig.existsVar("HorzOrder"))
      Err += PGradConfig.get("HorzOrder", HorzOrder);

   std::string VertReconStr = "linear";
   if (PGradConfig.existsVar("VerticalReconstruction"))
      Err += PGradConfig.get("VerticalReconstruction", VertReconStr);

   I4 QuadraturePoints = FiniteVolumePGrad.QuadraturePoints;
   if (PGradConfig.existsVar("QuadraturePoints"))
      Err += PGradConfig.get("QuadraturePoints", QuadraturePoints);

   CHECK_ERROR_ABORT(Err, "PressureGrad: error reading PressureGrad options");

   // Validate. Values that are not implemented are rejected outright rather
   // than silently falling back to an implemented one, which would make a
   // run requesting them look like a run that used them.
   if (HorzOrder != 2)
      ABORT_ERROR("PressureGrad: HorzOrder {} is not implemented for the "
                  "FiniteVolume pressure gradient; only HorzOrder 2 "
                  "(the two-cell stencil) is available",
                  HorzOrder);

   if (VertReconStr == "linear" || VertReconStr == "Linear") {
      FiniteVolumePGrad.VertRecon = PressureGradVertRecon::Linear;
   } else if (VertReconStr == "ppm" || VertReconStr == "PPM") {
      ABORT_ERROR("PressureGrad: VerticalReconstruction '{}' is not "
                  "implemented for the FiniteVolume pressure gradient; "
                  "only 'linear' is available",
                  VertReconStr);
   } else {
      ABORT_ERROR("PressureGrad: unknown VerticalReconstruction '{}'; only "
                  "'linear' is available",
                  VertReconStr);
   }

   if (QuadraturePoints < 1 || QuadraturePoints > MaxPGradQuadPoints)
      ABORT_ERROR("PressureGrad: QuadraturePoints {} is out of range; must be "
                  "between 1 and {}",
                  QuadraturePoints, MaxPGradQuadPoints);

   FiniteVolumePGrad.HorzOrder        = HorzOrder;
   FiniteVolumePGrad.QuadraturePoints = QuadraturePoints;

   // Additional mesh and coordinate data used by the FiniteVolume column scan
   CellsOnEdge  = Mesh->CellsOnEdge;
   MinLayerCell = VCoord->MinLayerCell;
   MaxLayerCell = VCoord->MaxLayerCell;
   NCellsAll    = Mesh->NCellsAll;

   // Working arrays for the FiniteVolume scheme, allocated only when that
   // scheme is selected so that a Centered run pays nothing for them
   if (PressureGradChoice == PressureGradType::FiniteVolume) {
      ReconSlopeCt =
          Array2DReal("PGradReconSlopeCt", Mesh->NCellsSize, NVertLayers);
      ReconSlopeSa =
          Array2DReal("PGradReconSlopeSa", Mesh->NCellsSize, NVertLayers);
      DeltaZIncr =
          Array2DReal("PGradDeltaZIncr", Mesh->NEdgesSize, NVertLayers);
      DeltaZMoment =
          Array2DReal("PGradDeltaZMoment", Mesh->NEdgesSize, NVertLayers);
      DeltaZFixedP =
          Array2DReal("PGradDeltaZFixedP", Mesh->NEdgesSize, NVertLayersP1);
   }

   // Temporary: initialization of tidal potential and SAL
   TidalPotential = Array1DReal("TidalPotential", Mesh->NCellsSize);
   SelfAttractionLoading =
       Array1DReal("SelfAttractionLoading", Mesh->NCellsSize);

} // end constructor

//------------------------------------------------------------------------------
// Destructor for PressureGrad
PressureGrad::~PressureGrad() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Remove PressureGrad instances before exit
void PressureGrad::clear() {

   AllPGrads.clear();

   DefaultPGrad = nullptr; // prevent dangling pointe'r

} // end clear

//------------------------------------------------------------------------------
// Remove PressureGrad from list by name
void PressureGrad::erase(const std::string &Name) {

   AllPGrads.erase(Name);

} // end erase

//------------------------------------------------------------------------------
// Get pressure gradient instance by name
PressureGrad *PressureGrad::get(const std::string &Name ///< [in] Name of
) {

   auto it = AllPGrads.find(Name);

   if (it != AllPGrads.end()) {
      return it->second.get();
   } else {
      LOG_ERROR("PressureGrad::get: Attempt to retrieve non-existent "
                "PressureGrad:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get pressure gradient

//------------------------------------------------------------------------------
// Compute the mean-preserving reconstruction slopes of temperature and
// salinity in pressure, once per cell and layer.
void PressureGrad::computeReconSlopes(
    const Array2DReal &ConservTemp, ///< [in] layer-mean temperature
    const Array2DReal &AbsSalinity, ///< [in] layer-mean salinity
    const Array2DReal &PressureMid  ///< [in] mid-layer pressure
) const {

   OMEGA_SCOPE(LocReconSlopeCt, ReconSlopeCt);
   OMEGA_SCOPE(LocReconSlopeSa, ReconSlopeSa);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);

   parallelForOuter(
       "pgrad-recon-slopes", {NCellsAll},
       KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
          const I4 KMin = LocMinLayerCell(ICell);
          const I4 KMax = LocMaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 I4 KLo, KHi;
                 linearReconStencil(K, KMin, KMax, KLo, KHi);

                 LocReconSlopeCt(ICell, K) = linearReconSlope(
                     ConservTemp(ICell, KLo), ConservTemp(ICell, KHi),
                     PressureMid(ICell, KLo), PressureMid(ICell, KHi));
                 LocReconSlopeSa(ICell, K) = linearReconSlope(
                     AbsSalinity(ICell, KLo), AbsSalinity(ICell, KHi),
                     PressureMid(ICell, KLo), PressureMid(ICell, KHi));
              });
       });

} // end computeReconSlopes

//------------------------------------------------------------------------------
// Accumulate the fixed-pressure height difference down each edge's column.
//
// The edge control volume is taken as the average of the two columns'
// interface pressures, which makes its mid-layer pressure and its pressure
// thickness exactly the edge averages of the two columns' own.
//
// The scan runs in two passes, both inside one team per edge with a barrier
// between them, so each covers only the edge's active layers. The first
// evaluates the matched-pressure integrand at the quadrature points of each
// edge layer and forms two integrals of it over the same points: the increment
// of the recurrence, and its first moment about the layer's top interface,
// which the layer mean needs. The second anchors the column and turns the
// increments into the height difference at every interface.
//
// The anchor is at the sea floor, for exactness rather than conditioning
// (design section 3.7.4). VertCoord builds geometric height upward from a
// prescribed bathymetry with a midpoint rule over each column's own layers,
// and two columns with different layer partitions accumulate different
// amounts of its error, even on profiles the reconstruction resolves exactly.
// Anchored at the surface, that difference would enter the height difference
// directly; anchored at the sea floor over a flat floor, the height difference
// is exact input and vanishes identically. Where the floor steps, the choice
// of end is a trade-off rather than a settled preference.
//
// The anchor is computed, not assumed -- it is whatever the model's geometric
// heights and interface pressures imply, evaluated at a common pressure, and
// it is the deepest instance of the same comparison the recurrence makes at
// every other interface.
//
// Where the two floors are more than a layer apart in pressure, as with
// terrain-following layers over a steep slope, the comparison is made at the
// shallower floor against the deeper column's own interface just below it, so
// that no shift to a common pressure spans more than one layer. One
// first-order expansion stretched across the whole floor offset would leave an
// error that grows with the cube of the offset and is the same at every depth.
void PressureGrad::computeColumnScan(
    const Array2DReal &PressureMid,       ///< [in] mid-layer pressure
    const Array2DReal &PressureInterface, ///< [in] interface pressure
    const Array2DReal &GeomZInterface,    ///< [in] interface geometric height
    const Array2DReal &ConservTemp,       ///< [in] layer-mean temperature
    const Array2DReal &AbsSalinity,       ///< [in] layer-mean salinity
    const Eos *EqState                    ///< [in] equation of state
) const {

   const Array2DReal SpecVol    = EqState->SpecVol;
   const Array2DReal SpecVolDCt = EqState->SpecVolDCt;
   const Array2DReal SpecVolDSa = EqState->SpecVolDSa;
   const Array2DReal SpecVolDP  = EqState->SpecVolDP;
   const I4 NQuad               = FiniteVolumePGrad.QuadraturePoints;
   const Real InvGravity        = 1.0_Real / Gravity;

   OMEGA_SCOPE(LocReconSlopeCt, ReconSlopeCt);
   OMEGA_SCOPE(LocReconSlopeSa, ReconSlopeSa);
   OMEGA_SCOPE(LocDeltaZIncr, DeltaZIncr);
   OMEGA_SCOPE(LocDeltaZMoment, DeltaZMoment);
   OMEGA_SCOPE(LocDeltaZFixedP, DeltaZFixedP);
   OMEGA_SCOPE(LocCellsOnEdge, CellsOnEdge);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);

   parallelForOuter(
       "pgrad-fv-column-scan", {NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KTop = LocMinLayerEdgeBot(IEdge);
          const I4 KBot = LocMaxLayerEdgeTop(IEdge);
          if (KBot < KTop)
             return;

          const I4 ICell0 = LocCellsOnEdge(IEdge, 0);
          const I4 ICell1 = LocCellsOnEdge(IEdge, 1);

          Real Nodes[MaxPGradQuadPoints];
          Real Weights[MaxPGradQuadPoints];
          gaussLegendreRule(NQuad, Nodes, Weights);

          // Pass one: the two integrals of the matched-pressure integrand
          // over each edge layer, from the same quadrature points.
          parallelForInner(
              Team, Range{KTop, KBot}, INNER_LAMBDA(int K) {
                 // one shared expansion for this edge layer, multiplying both
                 // columns
                 const PGradEdgeEos EdgeEos = buildEdgeEos(
                     SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, ConservTemp,
                     AbsSalinity, PressureMid, ICell0, ICell1, K);

                 const Real EdgeTop = 0.5_Real * (PressureInterface(ICell0, K) +
                                                  PressureInterface(ICell1, K));
                 const Real EdgeBot =
                     0.5_Real * (PressureInterface(ICell0, K + 1) +
                                 PressureInterface(ICell1, K + 1));
                 const Real EdgeMid  = 0.5_Real * (EdgeTop + EdgeBot);
                 const Real EdgeHalf = 0.5_Real * (EdgeBot - EdgeTop);

                 Real Incr   = 0.0_Real;
                 Real Moment = 0.0_Real;

                 for (int IQuad = 0; IQuad < NQuad; ++IQuad) {

                    const Real Press  = EdgeMid + EdgeHalf * Nodes[IQuad];
                    const Real Weight = EdgeHalf * Weights[IQuad];

                    // each column supplies its own temperature and salinity
                    // from whichever of its own layers contains this
                    // pressure, which under tilt is generally not layer K
                    const I4 KFound0 = findLayerForPress(
                        PressureInterface, ICell0, LocMinLayerCell(ICell0),
                        LocMaxLayerCell(ICell0), Press, K);
                    const I4 KFound1 = findLayerForPress(
                        PressureInterface, ICell1, LocMinLayerCell(ICell1),
                        LocMaxLayerCell(ICell1), Press, K);

                    const Real Temp0 =
                        linearReconEval(ConservTemp(ICell0, KFound0),
                                        LocReconSlopeCt(ICell0, KFound0),
                                        PressureMid(ICell0, KFound0), Press);
                    const Real Salt0 =
                        linearReconEval(AbsSalinity(ICell0, KFound0),
                                        LocReconSlopeSa(ICell0, KFound0),
                                        PressureMid(ICell0, KFound0), Press);
                    const Real Temp1 =
                        linearReconEval(ConservTemp(ICell1, KFound1),
                                        LocReconSlopeCt(ICell1, KFound1),
                                        PressureMid(ICell1, KFound1), Press);
                    const Real Salt1 =
                        linearReconEval(AbsSalinity(ICell1, KFound1),
                                        LocReconSlopeSa(ICell1, KFound1),
                                        PressureMid(ICell1, KFound1), Press);

                    const Real SpecVolDiff = matchedPressSpecVolDiff(
                        EdgeEos, Temp0, Salt0, Temp1, Salt1);

                    Incr += Weight * SpecVolDiff;
                    Moment += Weight * (Press - EdgeTop) * SpecVolDiff;
                 }

                 LocDeltaZIncr(IEdge, K)   = InvGravity * Incr;
                 LocDeltaZMoment(IEdge, K) = InvGravity * Moment;
              });

          // the scan below reads increments other threads of the team wrote
          teamBarrier(Team);

          // Pass two: the anchor at the sea floor, then the recurrence upward.
          //
          // The shallower column's floor is the deepest interface the two
          // columns share. The deeper column is represented by its own
          // interface at or just below that floor's pressure: the shared
          // interface itself wherever the floors are within one layer of each
          // other, as on a flat floor, and a shallower interface of its own
          // where they are further apart. Terrain-following layers over a
          // steep slope put the floors several layers apart, and partial
          // bottom cells put them just over a layer apart wherever a thin
          // bottom cell lies entirely above the deeper column's layer of the
          // same index. Each column's height is shifted from its own
          // interface pressure to the pressure midway between the two, by
          // integrating its own reconstruction with one expansion shared by
          // the two layers that meet there, and both short integrals vanish
          // where the two interface pressures agree.
          const Real FloorPress0 = PressureInterface(ICell0, KBot + 1);
          const Real FloorPress1 = PressureInterface(ICell1, KBot + 1);
          const bool Shallow0    = FloorPress0 <= FloorPress1;
          const I4 IDeep         = Shallow0 ? ICell1 : ICell0;
          const I4 KDeep =
              findLayerForPress(PressureInterface, IDeep,
                                LocMinLayerCell(IDeep), LocMaxLayerCell(IDeep),
                                Shallow0 ? FloorPress0 : FloorPress1, KBot);
          const I4 KAnchor0 = Shallow0 ? KBot : KDeep;
          const I4 KAnchor1 = Shallow0 ? KDeep : KBot;

          const Real AnchorPress =
              0.5_Real * (PressureInterface(ICell0, KAnchor0 + 1) +
                          PressureInterface(ICell1, KAnchor1 + 1));

          const PGradEdgeEos AnchorEos = buildEdgeEos(
              SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, ConservTemp,
              AbsSalinity, PressureMid, ICell0, KAnchor0, ICell1, KAnchor1);

          Real Anchor = GeomZInterface(ICell1, KAnchor1 + 1) -
                        GeomZInterface(ICell0, KAnchor0 + 1);

          for (int ISide = 0; ISide < 2; ++ISide) {

             const I4 ICell   = (ISide == 0) ? ICell0 : ICell1;
             const I4 KAnchor = (ISide == 0) ? KAnchor0 : KAnchor1;
             const Real Sign  = (ISide == 0) ? -1.0_Real : 1.0_Real;

             const Real ColPress = PressureInterface(ICell, KAnchor + 1);
             const Real Mid      = 0.5_Real * (AnchorPress + ColPress);
             const Real Half     = 0.5_Real * (ColPress - AnchorPress);

             Real Integral = 0.0_Real;
             for (int IQuad = 0; IQuad < NQuad; ++IQuad) {

                const Real Press  = Mid + Half * Nodes[IQuad];
                const Real Weight = Half * Weights[IQuad];

                const I4 KFound = findLayerForPress(
                    PressureInterface, ICell, LocMinLayerCell(ICell),
                    LocMaxLayerCell(ICell), Press, KAnchor);

                const Real Temp = linearReconEval(
                    ConservTemp(ICell, KFound), LocReconSlopeCt(ICell, KFound),
                    PressureMid(ICell, KFound), Press);
                const Real Salt = linearReconEval(
                    AbsSalinity(ICell, KFound), LocReconSlopeSa(ICell, KFound),
                    PressureMid(ICell, KFound), Press);

                Integral += Weight * edgeSpecVol(AnchorEos, Temp, Salt, Press);
             }

             Anchor += Sign * InvGravity * Integral;
          }

          // Where the anchor sits above the edge's deepest interface, carry it
          // down there with the matched-pressure integrand, which vanishes on
          // the exact set like every other increment of the scan.
          if (KDeep != KBot) {
             const PGradEdgeEos FloorEos = buildEdgeEos(
                 SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, ConservTemp,
                 AbsSalinity, PressureMid, ICell0, ICell1, KBot);
             const Real FloorPress = 0.5_Real * (FloorPress0 + FloorPress1);
             const Real Mid        = 0.5_Real * (AnchorPress + FloorPress);
             const Real Half       = 0.5_Real * (FloorPress - AnchorPress);

             Real Integral = 0.0_Real;
             for (int IQuad = 0; IQuad < NQuad; ++IQuad) {

                const Real Press  = Mid + Half * Nodes[IQuad];
                const Real Weight = Half * Weights[IQuad];

                const I4 KFound0 = findLayerForPress(
                    PressureInterface, ICell0, LocMinLayerCell(ICell0),
                    LocMaxLayerCell(ICell0), Press, KBot);
                const I4 KFound1 = findLayerForPress(
                    PressureInterface, ICell1, LocMinLayerCell(ICell1),
                    LocMaxLayerCell(ICell1), Press, KBot);

                const Real Temp0 =
                    linearReconEval(ConservTemp(ICell0, KFound0),
                                    LocReconSlopeCt(ICell0, KFound0),
                                    PressureMid(ICell0, KFound0), Press);
                const Real Salt0 =
                    linearReconEval(AbsSalinity(ICell0, KFound0),
                                    LocReconSlopeSa(ICell0, KFound0),
                                    PressureMid(ICell0, KFound0), Press);
                const Real Temp1 =
                    linearReconEval(ConservTemp(ICell1, KFound1),
                                    LocReconSlopeCt(ICell1, KFound1),
                                    PressureMid(ICell1, KFound1), Press);
                const Real Salt1 =
                    linearReconEval(AbsSalinity(ICell1, KFound1),
                                    LocReconSlopeSa(ICell1, KFound1),
                                    PressureMid(ICell1, KFound1), Press);

                Integral += Weight * matchedPressSpecVolDiff(
                                         FloorEos, Temp0, Salt0, Temp1, Salt1);
             }

             Anchor -= InvGravity * Integral;
          }

          LocDeltaZFixedP(IEdge, KBot + 1) = Anchor;

          // Accumulate upward, in the same shape as VertCoord's geometric
          // height. Every quantity here is a horizontal contrast, so no large
          // quantity is formed and nothing large has to cancel.
          parallelScanInner(
              Team, vertRange(KTop, KBot),
              INNER_LAMBDA(int K, Real &Accum, bool IsFinal) {
                 const I4 KLyr = KBot - K;
                 Accum += LocDeltaZIncr(IEdge, KLyr);
                 if (IsFinal)
                    LocDeltaZFixedP(IEdge, KLyr) = Anchor + Accum;
              });
       });

} // end computeColumnScan

//------------------------------------------------------------------------------
// Compute pressure gradient tendencies and add into Tend array
void PressureGrad::computePressureGrad(
    Array2DReal &Tend, const Array2DReal &PressureMid,
    const Array2DReal &PressureInterface, const Array2DReal &SpecVol,
    const Array2DReal &GeomZInterface, const Array2DReal &PseudoThick,
    const Array2DReal &ConservTemp, const Array2DReal &AbsSalinity,
    const Eos *EqState) const {

   OMEGA_SCOPE(LocCenteredPGrad, CenteredPGrad);
   OMEGA_SCOPE(LocFiniteVolumePGrad, FiniteVolumePGrad);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);
   OMEGA_SCOPE(LocTidalPotential, TidalPotential);
   OMEGA_SCOPE(LocSelfAttractionLoading, SelfAttractionLoading);

   if (PressureGradChoice == PressureGradType::Centered) {

      // computes centered geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-centered", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin   = LocMinLayerEdgeBot(IEdge);
             const int KMax   = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocCenteredPGrad(Tend, IEdge, KChunk, PressureMid,
                                     PressureInterface, GeomZInterface,
                                     LocTidalPotential,
                                     LocSelfAttractionLoading, SpecVol);
                 });
          });

   } else {

      // The per-cell reconstruction slopes are formed once and reused across
      // each cell's edges; the per-edge work is polynomial arithmetic on them
      computeReconSlopes(ConservTemp, AbsSalinity, PressureMid);

      // The column scan is a prefix sum with edge-dependent coefficients, so
      // it cannot live inside the per-vertical-chunk functor below
      computeColumnScan(PressureMid, PressureInterface, GeomZInterface,
                        ConservTemp, AbsSalinity, EqState);

      OMEGA_SCOPE(LocDeltaZFixedP, DeltaZFixedP);
      OMEGA_SCOPE(LocDeltaZMoment, DeltaZMoment);

      // computes finite-volume geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-finitevolume", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin   = LocMinLayerEdgeBot(IEdge);
             const int KMax   = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocFiniteVolumePGrad(Tend, IEdge, KChunk, PressureInterface,
                                         LocDeltaZFixedP, LocDeltaZMoment,
                                         LocTidalPotential,
                                         LocSelfAttractionLoading);
                 });
          });
   }
} // end compute pressure gradient

//------------------------------------------------------------------------------
// Compute the barotropic pressure anomaly gradient and add into Tend array
void PressureGrad::computeBarotropicPressureGrad(
    Array2DReal &Tend,                  ///< [inout] velocity tendency
    const Array1DReal &BtrPressAnomaly, ///< [in] barotropic pressure anomaly
    const Array1DReal &DepthMeanSpecVol ///< [in] depth-mean specific volume
) const {

   OMEGA_SCOPE(LocBarotropicPGrad, BarotropicPGrad);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);

   parallelForOuter(
       "pgrad-barotropic", {NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const int KMin   = LocMinLayerEdgeBot(IEdge);
          const int KMax   = LocMaxLayerEdgeTop(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocBarotropicPGrad(Tend, IEdge, KChunk, BtrPressAnomaly,
                                    DepthMeanSpecVol);
              });
       });
} // end compute barotropic pressure gradient

//------------------------------------------------------------------------------
// Constructor for centered pressure gradient functor
PressureGradCentered::PressureGradCentered(
    const HorzMesh *Mesh,   ///< [in] Horizontal mesh
    const VertCoord *VCoord ///< [in] Vertical coordinate
    )
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

//------------------------------------------------------------------------------
// Constructor for finite volume pressure gradient functor
PressureGradFiniteVolume::PressureGradFiniteVolume(
    const HorzMesh *Mesh,   ///< [in] Horizontal mesh
    const VertCoord *VCoord ///< [in] Vertical coordinate
    )
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

//------------------------------------------------------------------------------
// Constructor for barotropic pressure gradient functor
BarotropicPressureGradOnEdge::BarotropicPressureGradOnEdge(
    const HorzMesh *Mesh,   ///< [in] Horizontal mesh
    const VertCoord *VCoord ///< [in] Vertical coordinate
    )
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

} // namespace OMEGA
