#include "HorzOperators.h"
#include "DataTypes.h"
#include "HorzMesh.h"

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
    : OnSphere(true), NCellsAll(Mesh->NCellsAll), MaxEdges(1 + Mesh->MaxEdges),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnCell(Mesh->CellsOnCell), CellsOnEdge(Mesh->CellsOnEdge),
      VerticesOnEdge(Mesh->VerticesOnEdge),

      XCell(Mesh->XCell), YCell(Mesh->YCell),
      ZCell(createDeviceMirrorCopy(Mesh->ZCellH)), DvEdge(Mesh->DvEdge),
      DcEdge(Mesh->DcEdge), AngleEdge(Mesh->AngleEdge),
      AreaCell(Mesh->AreaCell), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      XVertex(createDeviceMirrorCopy(Mesh->XVertexH)),
      YVertex(createDeviceMirrorCopy(Mesh->YVertexH)),
      ZVertex(createDeviceMirrorCopy(Mesh->ZVertexH)),

      ThetaAbs("SphereAngle", NCellsAll), XPCell("XP", NCellsAll, MaxEdges),
      YPCell("YP", NCellsAll, MaxEdges),
      Angle2DCell("Angle2D", NCellsAll, MaxEdges),
      BCell("WorkSpaceForLeastSquares", NCellsAll, 6, MaxEdges),
      CellListCell("CellList", NCellsAll, MaxEdges) {
   if (MaxMaxEdges <= Mesh->MaxEdges)
      LOG_CRITICAL(
          "SecondDerivativeOnCell::SecondDerivativeOnCell Max Edges exceeded:"
          "Max Allowed: {}  Found in Mesh:{}",
          MaxMaxEdges - 1, Mesh->MaxEdges);
}

MasksAndCoefficients::MasksAndCoefficients(
    HorzMesh const *Mesh, const Array3DReal DerivTwo,
    Array1DI4 NAdvCellsForEdge, Array2DI4 AdvCellsForEdge,
    Array1DI4 AdvMaskHighOrder, Array2DReal AdvCoefs, Array2DReal AdvCoefs3rd)
    : NCellsGlobal(Mesh->NCellsGlobal), NCellsAll(Mesh->NCellsAll),
      NAdvCellsMax(Mesh->MaxEdges2), // PatchCellLists("PatchCellLists",
                                     // Mesh->NEdgesOwned, Mesh->NEdgesAll+1),
      NAdvCellsForEdge(NAdvCellsForEdge), AdvCellsForEdge(AdvCellsForEdge),
      NEdgesOnEdge(Mesh->NEdgesOnEdge), NEdgesOnCell(Mesh->NEdgesOnCell),
      CellIndx("CellIndx", Mesh->NEdgesOwned, Mesh->MaxEdges2 + 2),
      CellIndxSorted("CellIndxSorted", Mesh->NEdgesOwned, 2,
                     Mesh->MaxEdges2 + 2),
      CellID(Mesh->CellID), AdvMaskHighOrder(AdvMaskHighOrder),
      EdgesOnEdge(Mesh->EdgesOnEdge), CellsOnCell(Mesh->CellsOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge), AdvCoefs(AdvCoefs), AdvCoefs3rd(AdvCoefs3rd),
      DerivTwo(DerivTwo) {}
} // namespace OMEGA
