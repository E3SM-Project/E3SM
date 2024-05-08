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

} // namespace OMEGA
