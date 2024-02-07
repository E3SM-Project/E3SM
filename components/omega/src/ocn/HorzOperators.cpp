#include "HorzOperators.h"
#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {

DivergenceOnCell::DivergenceOnCell(HorzMesh const *mesh)
    : NEdgesOnCell(mesh->NEdgesOnCell), EdgesOnCell(mesh->EdgesOnCell),
      DvEdge(mesh->DvEdge), AreaCell(mesh->AreaCell),
      EdgeSignOnCell(mesh->EdgeSignOnCell) {}

GradientOnEdge::GradientOnEdge(HorzMesh const *mesh)
    : CellsOnEdge(mesh->CellsOnEdge), DcEdge(mesh->DcEdge) {}

CurlOnVertex::CurlOnVertex(HorzMesh const *mesh)
    : VertexDegree(mesh->VertexDegree), EdgesOnVertex(mesh->EdgesOnVertex),
      DcEdge(mesh->DcEdge), AreaTriangle(mesh->AreaTriangle),
      EdgeSignOnVertex(mesh->EdgeSignOnVertex) {}

TangentialReconOnEdge::TangentialReconOnEdge(HorzMesh const *mesh)
    : NEdgesOnEdge(mesh->NEdgesOnEdge), EdgesOnEdge(mesh->EdgesOnEdge),
      WeightsOnEdge(mesh->WeightsOnEdge) {}

} // namespace OMEGA
