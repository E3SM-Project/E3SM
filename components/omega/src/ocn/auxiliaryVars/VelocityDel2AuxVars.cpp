#include "VelocityDel2AuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

const std::string VelocityDel2AuxVars::Del2EdgeName    = "VelDel2Edge";
const std::string VelocityDel2AuxVars::Del2DivCellName = "VelDel2DivCell";
const std::string VelocityDel2AuxVars::Del2RelVortVertexName =
    "VelDel2RelVortVertex";

VelocityDel2AuxVars::VelocityDel2AuxVars(const HorzMesh *mesh, int NVertLevels)
    : Del2Edge("VelDel2Edge", mesh->NEdgesSize, NVertLevels),
      Del2DivCell("VelDel2DivCell", mesh->NCellsSize, NVertLevels),
      Del2RelVortVertex("VelDel2RelVortVertex", mesh->NVerticesSize,
                        NVertLevels),
      NEdgesOnCell(mesh->NEdgesOnCell), EdgesOnCell(mesh->EdgesOnCell),
      EdgeSignOnCell(mesh->EdgeSignOnCell), DcEdge(mesh->DcEdge),
      DvEdge(mesh->DvEdge), AreaCell(mesh->AreaCell),
      EdgesOnVertex(mesh->EdgesOnVertex), CellsOnEdge(mesh->CellsOnEdge),
      VerticesOnEdge(mesh->VerticesOnEdge),
      EdgeSignOnVertex(mesh->EdgeSignOnVertex),
      AreaTriangle(mesh->AreaTriangle), VertexDegree(mesh->VertexDegree) {
   addMetaData();
   defineIOFields();
}

void VelocityDel2AuxVars::addMetaData() const {
   auto EdgeDim      = MetaDim::get("NEdges");
   auto CellDim      = MetaDim::get("NCells");
   auto VertexDim    = MetaDim::get("NVertices");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get("auxiliaryVars");

   const Real FillValue = -9.99e30;

   auto Del2EdgeMeta = ArrayMetaData::create(
       Del2EdgeName,
       "laplacian of horizontal velocity on edges", /// long Name or
                                                    /// description
       "m^-1 s^-1",                                 /// units
       "",                                          /// CF standard Name
       std::numeric_limits<Real>::min(),            /// min valid value
       std::numeric_limits<Real>::max(),            /// max valid value
       FillValue,         /// scalar used for undefined entries
       2,                 /// number of dimensions
       {EdgeDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(Del2EdgeName);

   auto Del2DivCellMeta =
       ArrayMetaData::create(Del2DivCellName,
                             "divergence of laplacian of horizontal velocity "
                             "on cells",  /// long Name or description
                             "m^-2 s^-1", /// units
                             "",          /// CF standard Name
                             std::numeric_limits<Real>::min(), /// min valid
                                                               /// value
                             std::numeric_limits<Real>::max(), /// max valid
                                                               /// value
                             FillValue, /// scalar used for undefined entries
                             2,         /// number of dimensions
                             {CellDim, VertDim} /// dim pointers
       );
   AuxMetaGroup->addField(Del2DivCellName);

   auto Del2RelVortVertexMeta = ArrayMetaData::create(
       Del2RelVortVertexName,
       "curl of laplacian of horizontal velocity on cells", /// long Name or
                                                            /// description
       "m^-2 s^-1",                                         /// units
       "",                                                  /// CF standard Name
       std::numeric_limits<Real>::min(),                    /// min valid value
       std::numeric_limits<Real>::max(),                    /// max valid value
       FillValue,           /// scalar used for undefined entries
       2,                   /// number of dimensions
       {VertexDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(Del2RelVortVertexName);
}

void VelocityDel2AuxVars::defineIOFields() const {
   int Err;

   Err = IOField::define(Del2EdgeName);
   Err = IOField::attachData(Del2EdgeName, Del2Edge);

   Err = IOField::define(Del2DivCellName);
   Err = IOField::attachData(Del2DivCellName, Del2DivCell);

   Err = IOField::define(Del2RelVortVertexName);
   Err = IOField::attachData(Del2RelVortVertexName, Del2RelVortVertex);
}

} // namespace OMEGA
