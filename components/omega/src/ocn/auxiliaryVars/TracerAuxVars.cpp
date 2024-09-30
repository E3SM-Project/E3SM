#include "TracerAuxVars.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

TracerAuxVars::TracerAuxVars(const std::string &AuxStateSuffix,
                             const HorzMesh *Mesh, const I4 NVertLevels,
                             const I4 NTracers)
    : HTracersOnEdge("ThickTracersOnEdge" + AuxStateSuffix, NTracers,
                     Mesh->NEdgesSize, NVertLevels),
      Del2TracersOnCell("Del2TracerOnCell" + AuxStateSuffix, NTracers,
                        Mesh->NCellsSize, NVertLevels),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell) {}

} // namespace OMEGA
