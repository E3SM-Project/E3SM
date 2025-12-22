//===-- ocn/TendencyTerms.cpp - Tendency Terms ------------------*- C++ -*-===//
//
// The tendency terms that update state variables are implemented as functors,
// i.e. as classes that act like functions. This source defines the class
// constructors for these functors, which initialize the functor objects using
// the Mesh objects and info from the Config. The function call operators () are
// defined in the corresponding header file.
//
//===----------------------------------------------------------------------===//

#include "TendencyTerms.h"
#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "Tracers.h"

namespace OMEGA {

ThicknessFluxDivOnCell::ThicknessFluxDivOnCell(const HorzMesh *Mesh,
                                               const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

PotentialVortHAdvOnEdge::PotentialVortHAdvOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

KEGradOnEdge::KEGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SSHGradOnEdge::SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityDiffusionOnEdge::VelocityDiffusionOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityHyperDiffOnEdge::VelocityHyperDiffOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

WindForcingOnEdge::WindForcingOnEdge(const HorzMesh *Mesh,
                                     const VertCoord *VCoord)
    : Enabled(false), LocRhoSw(RhoSw), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot) {}

BottomDragOnEdge::BottomDragOnEdge(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : Enabled(false), Coeff(0), CellsOnEdge(Mesh->CellsOnEdge),
      NVertLayers(VCoord->NVertLayers), EdgeMask(VCoord->EdgeMask),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerHorzAdvOnCell::TracerHorzAdvOnCell(const HorzMesh *Mesh,
                                         const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeMask(VCoord->EdgeMask), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerDiffOnCell::TracerDiffOnCell(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(VCoord->EdgeMask),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerHyperDiffOnCell::TracerHyperDiffOnCell(const HorzMesh *Mesh,
                                             const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(VCoord->EdgeMask),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
