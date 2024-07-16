#include "VelocityDel2AuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

VelocityDel2AuxVars::VelocityDel2AuxVars(const std::string &AuxStateSuffix,
                                         const HorzMesh *Mesh, int NVertLevels)
    : Del2Edge("VelDel2Edge" + AuxStateSuffix, Mesh->NEdgesSize, NVertLevels),
      Del2DivCell("VelDel2DivCell" + AuxStateSuffix, Mesh->NCellsSize,
                  NVertLevels),
      Del2RelVortVertex("VelDel2RelVortVertex" + AuxStateSuffix,
                        Mesh->NVerticesSize, NVertLevels),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgesOnVertex(Mesh->EdgesOnVertex), CellsOnEdge(Mesh->CellsOnEdge),
      VerticesOnEdge(Mesh->VerticesOnEdge),
      EdgeSignOnVertex(Mesh->EdgeSignOnVertex),
      AreaTriangle(Mesh->AreaTriangle), VertexDegree(Mesh->VertexDegree) {}

void VelocityDel2AuxVars::registerFields(
    const std::string &AuxGroupName) const {
   addMetaData(AuxGroupName);
   defineIOFields();
}

void VelocityDel2AuxVars::unregisterFields() const {
   IOField::erase(Del2Edge.label());
   IOField::erase(Del2DivCell.label());
   IOField::erase(Del2RelVortVertex.label());
   MetaData::destroy(Del2Edge.label());
   MetaData::destroy(Del2DivCell.label());
   MetaData::destroy(Del2RelVortVertex.label());
}

void VelocityDel2AuxVars::addMetaData(const std::string &AuxGroupName) const {
   auto EdgeDim      = MetaDim::get("NEdges");
   auto CellDim      = MetaDim::get("NCells");
   auto VertexDim    = MetaDim::get("NVertices");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get(AuxGroupName);

   const Real FillValue = -9.99e30;

   auto Del2EdgeMeta = ArrayMetaData::create(
       Del2Edge.label(),
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
   AuxMetaGroup->addField(Del2Edge.label());

   auto Del2DivCellMeta =
       ArrayMetaData::create(Del2DivCell.label(),
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
   AuxMetaGroup->addField(Del2DivCell.label());

   auto Del2RelVortVertexMeta = ArrayMetaData::create(
       Del2RelVortVertex.label(),
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
   AuxMetaGroup->addField(Del2RelVortVertex.label());
}

void VelocityDel2AuxVars::defineIOFields() const {
   int Err;

   Err = IOField::define(Del2Edge.label());
   Err = IOField::attachData(Del2Edge.label(), Del2Edge);

   Err = IOField::define(Del2DivCell.label());
   Err = IOField::attachData(Del2DivCell.label(), Del2DivCell);

   Err = IOField::define(Del2RelVortVertex.label());
   Err = IOField::attachData(Del2RelVortVertex.label(), Del2RelVortVertex);
}

} // namespace OMEGA
