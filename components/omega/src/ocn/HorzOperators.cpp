#include "HorzOperators.h"
#include "DataTypes.h"
#include "Error.h"
#include "HorzMesh.h"
#include "VertCoord.h"

namespace OMEGA {

DivergenceOnCell::DivergenceOnCell(HorzMesh const *Mesh)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell) {}

GradientOnEdge::GradientOnEdge(HorzMesh const *Mesh)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge) {}

CurlOnVertex::CurlOnVertex(HorzMesh const *Mesh)
    : VertexDegree(Mesh->VertexDegree), EdgesOnVertex(Mesh->EdgesOnVertex),
      DcEdge(Mesh->DcEdge), AreaTriangle(Mesh->AreaTriangle),
      EdgeSignOnVertex(Mesh->EdgeSignOnVertex) {}

TangentialReconOnEdge::TangentialReconOnEdge(HorzMesh const *Mesh)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge) {}

InterpCellToEdge::InterpCellToEdge(const HorzMesh *Mesh)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      CellsOnVertex(Mesh->CellsOnVertex),
      KiteAreasOnVertex(Mesh->KiteAreasOnVertex),
      VertexDegree(Mesh->VertexDegree) {}

SecondDerivativeOnCell::SecondDerivativeOnCell(HorzMesh const *Mesh)
    : OnSphere(Mesh->OnSphere), NCellsAll(Mesh->NCellsAll),
      MaxEdges(1 + Mesh->MaxEdges), NEdgesOnCell(Mesh->NEdgesOnCell),
      EdgesOnCell(Mesh->EdgesOnCell), CellsOnCell(Mesh->CellsOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      XCell(Mesh->XCell), YCell(Mesh->YCell), ZCell(Mesh->ZCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AngleEdge(Mesh->AngleEdge),
      AreaCell(Mesh->AreaCell), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      XVertex(Mesh->XVertex), YVertex(Mesh->YVertex), ZVertex(Mesh->ZVertex),
      ThetaAbs("SphereAngle", NCellsAll), XPCell("XP", NCellsAll, MaxEdges),
      YPCell("YP", NCellsAll, MaxEdges),
      Angle2DCell("Angle2D", NCellsAll, MaxEdges),
      BCell("WorkSpaceForLeastSquares", NCellsAll, 6, MaxEdges),
      CellListCell("CellList", NCellsAll, MaxEdges) {}

MasksAndCoefficients::MasksAndCoefficients(
    HorzMesh const *Mesh, VertCoord const *VCoord, const Array3DReal DerivTwo,
    Array1DI4 NAdvCellsForEdge, Array2DI4 AdvCellsForEdge,
    Array2DI4 AdvMaskHighOrder, Array2DReal AdvCoefs, Array2DReal AdvCoefs3rd)
    : NCellsGlobal(Mesh->NCellsGlobal), NCellsAll(Mesh->NCellsAll),
      NAdvCellsMax(Mesh->MaxEdges2), // PatchCellLists("PatchCellLists",
                                     // Mesh->NEdgesOwned, Mesh->NEdgesAll+1),
      NAdvCellsForEdge(NAdvCellsForEdge), AdvCellsForEdge(AdvCellsForEdge),
      NEdgesOnEdge(Mesh->NEdgesOnEdge), NEdgesOnCell(Mesh->NEdgesOnCell),
      CellIndx("CellIndx", Mesh->NEdgesAll, Mesh->MaxEdges2 + 2),
      CellIndxSorted("CellIndxSorted", Mesh->NEdgesAll, 2, Mesh->MaxEdges2 + 2),
      CellID(Mesh->CellID), AdvMaskHighOrder(AdvMaskHighOrder),
      EdgesOnEdge(Mesh->EdgesOnEdge), CellsOnCell(Mesh->CellsOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge), AdvCoefs(AdvCoefs), AdvCoefs3rd(AdvCoefs3rd),
      BoundaryCell(VCoord->BoundaryCell), DerivTwo(DerivTwo) {}

VectorReconOnCell::VectorReconOnCell(HorzMesh const *Mesh)
    : OnSphere(Mesh->OnSphere), NEdgesReconOnCell(Mesh->NEdgesReconOnCell),
      ReconStencilCell(Mesh->ReconStencilCell),
      ReconWeightsCell(Mesh->ReconWeightsCell), LatCell(Mesh->LatCell),
      LonCell(Mesh->LonCell) {
   if (!Mesh->OnSphere)
      ABORT_ERROR("VectorReconOnCell: reconstruction stencil/weights "
                  "are only available for spherical meshes");
}

} // namespace OMEGA
