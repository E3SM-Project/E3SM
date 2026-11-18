#include "TracerForcingVars.h"
#include "Eos.h"
#include "Field.h"
#include "Tracers.h"
#include "VertCoord.h"

#include <limits>

namespace OMEGA {

TracerForcingVars::TracerForcingVars(const std::string &Suffix,
                                     const HorzMesh *Mesh)
    : SnowFluxCell("SnowFlux" + Suffix, Mesh->NCellsSize),
      RainFluxCell("RainFlux" + Suffix, Mesh->NCellsSize),
      EvaporationFluxCell("EvaporationFlux" + Suffix, Mesh->NCellsSize),
      SeaIceFreshWaterFluxCell("SeaIceFreshWaterFlux" + Suffix,
                               Mesh->NCellsSize),
      IceRunoffFluxCell("IceRunoffFlux" + Suffix, Mesh->NCellsSize),
      RiverRunoffFluxCell("RiverRunoffFlux" + Suffix, Mesh->NCellsSize),
      LatentHeatFluxEvapCell("LatentHeatFluxEvap" + Suffix, Mesh->NCellsSize),
      SensibleHeatFluxCell("SensibleHeatFlux" + Suffix, Mesh->NCellsSize),
      LongWaveHeatFluxUpCell("LongWaveHeatFluxUp" + Suffix, Mesh->NCellsSize),
      LongWaveHeatFluxDownCell("LongWaveHeatFluxDown" + Suffix,
                               Mesh->NCellsSize),
      SeaIceHeatFluxCell("SeaIceHeatFlux" + Suffix, Mesh->NCellsSize),
      ShortWaveHeatFluxCell("ShortWaveHeatFlux" + Suffix, Mesh->NCellsSize),
      SeaIceSaltFluxCell("SeaIceSaltFlux" + Suffix, Mesh->NCellsSize),
      ExtraEnergyCell("ExtraEnergy" + Suffix, Mesh->NCellsSize),
      SurfaceTracerFluxCell("SurfaceTracerFlux" + Suffix,
                            Tracers::getNumTracers(), Mesh->NCellsSize) {
   deepCopy(SnowFluxCell, 0.0_Real);
   deepCopy(RainFluxCell, 0.0_Real);
   deepCopy(EvaporationFluxCell, 0.0_Real);
   deepCopy(SeaIceFreshWaterFluxCell, 0.0_Real);
   deepCopy(IceRunoffFluxCell, 0.0_Real);
   deepCopy(RiverRunoffFluxCell, 0.0_Real);
   deepCopy(LatentHeatFluxEvapCell, 0.0_Real);
   deepCopy(SensibleHeatFluxCell, 0.0_Real);
   deepCopy(LongWaveHeatFluxUpCell, 0.0_Real);
   deepCopy(LongWaveHeatFluxDownCell, 0.0_Real);
   deepCopy(SeaIceHeatFluxCell, 0.0_Real);
   deepCopy(ShortWaveHeatFluxCell, 0.0_Real);
   deepCopy(SeaIceSaltFluxCell, 0.0_Real);
   deepCopy(ExtraEnergyCell, 0.0_Real);
}

void TracerForcingVars::registerFields(const std::string &MeshName) const {
   const int NDims = 1;
   std::vector<std::string> DimNames(NDims);

   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NCells" + DimSuffix;

   auto SnowFluxField =
       Field::create(SnowFluxCell.label(), "snow freshwater flux", "kg m-2 s-1",
                     "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);
   auto RainFluxField =
       Field::create(RainFluxCell.label(), "rain freshwater flux", "kg m-2 s-1",
                     "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);
   auto EvaporationFluxField =
       Field::create(EvaporationFluxCell.label(), "evaporation freshwater flux",
                     "kg m-2 s-1", "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);
   auto SeaIceFreshWaterFluxField = Field::create(
       SeaIceFreshWaterFluxCell.label(), "sea-ice freshwater flux",
       "kg m-2 s-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto IceRunoffFluxField =
       Field::create(IceRunoffFluxCell.label(), "ice runoff freshwater flux",
                     "kg m-2 s-1", "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);
   auto RiverRunoffFluxField = Field::create(
       RiverRunoffFluxCell.label(), "river runoff freshwater flux",
       "kg m-2 s-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), NDims, DimNames);

   auto LatentHeatFluxEvapField =
       Field::create(LatentHeatFluxEvapCell.label(), "latent heat flux",
                     "W m-2", "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);
   auto SensibleHeatFluxField =
       Field::create(SensibleHeatFluxCell.label(), "sensible heat flux",
                     "W m-2", "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);
   auto LongWaveHeatFluxUpField = Field::create(
       LongWaveHeatFluxUpCell.label(), "upward longwave heat flux", "W m-2", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       NDims, DimNames);
   auto LongWaveHeatFluxDownField = Field::create(
       LongWaveHeatFluxDownCell.label(), "downward longwave heat flux", "W m-2",
       "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto SeaIceHeatFluxField =
       Field::create(SeaIceHeatFluxCell.label(), "sea-ice heat flux", "W m-2",
                     "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);
   auto ShortWaveHeatFluxField =
       Field::create(ShortWaveHeatFluxCell.label(), "shortwave heat flux",
                     "W m-2", "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);

   auto SeaIceSaltFluxField =
       Field::create(SeaIceSaltFluxCell.label(), "sea-ice salt flux",
                     "kg m-2 s-1", "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);

   const int NTracerDims = 2;
   std::vector<std::string> TracerFluxDims(NTracerDims);
   TracerFluxDims[0]           = "NTracers";
   TracerFluxDims[1]           = "NCells" + DimSuffix;
   auto SurfaceTracerFluxField = Field::create(
       SurfaceTracerFluxCell.label(),
       "surface tracer flux for KPP non-local mixing", "", "",
       std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
       NTracerDims, TracerFluxDims);
   SurfaceTracerFluxField->setOptionalRead(true);

   FieldGroup::addFieldToGroup(SnowFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(RainFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(EvaporationFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SeaIceFreshWaterFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(IceRunoffFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(RiverRunoffFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(LatentHeatFluxEvapCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SensibleHeatFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(LongWaveHeatFluxUpCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(LongWaveHeatFluxDownCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SeaIceHeatFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(ShortWaveHeatFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SeaIceSaltFluxCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(SurfaceTracerFluxCell.label(), "Forcing");

   SnowFluxField->attachData<Array1DReal>(SnowFluxCell);
   RainFluxField->attachData<Array1DReal>(RainFluxCell);
   EvaporationFluxField->attachData<Array1DReal>(EvaporationFluxCell);
   SeaIceFreshWaterFluxField->attachData<Array1DReal>(SeaIceFreshWaterFluxCell);
   IceRunoffFluxField->attachData<Array1DReal>(IceRunoffFluxCell);
   RiverRunoffFluxField->attachData<Array1DReal>(RiverRunoffFluxCell);
   LatentHeatFluxEvapField->attachData<Array1DReal>(LatentHeatFluxEvapCell);
   SensibleHeatFluxField->attachData<Array1DReal>(SensibleHeatFluxCell);
   LongWaveHeatFluxUpField->attachData<Array1DReal>(LongWaveHeatFluxUpCell);
   LongWaveHeatFluxDownField->attachData<Array1DReal>(LongWaveHeatFluxDownCell);
   SeaIceHeatFluxField->attachData<Array1DReal>(SeaIceHeatFluxCell);
   ShortWaveHeatFluxField->attachData<Array1DReal>(ShortWaveHeatFluxCell);
   SeaIceSaltFluxField->attachData<Array1DReal>(SeaIceSaltFluxCell);
   SurfaceTracerFluxField->attachData<Array2DReal>(SurfaceTracerFluxCell);
}

void TracerForcingVars::unregisterFields() const {
   Field::destroy(SnowFluxCell.label());
   Field::destroy(RainFluxCell.label());
   Field::destroy(EvaporationFluxCell.label());
   Field::destroy(SeaIceFreshWaterFluxCell.label());
   Field::destroy(IceRunoffFluxCell.label());
   Field::destroy(RiverRunoffFluxCell.label());
   Field::destroy(LatentHeatFluxEvapCell.label());
   Field::destroy(SensibleHeatFluxCell.label());
   Field::destroy(LongWaveHeatFluxUpCell.label());
   Field::destroy(LongWaveHeatFluxDownCell.label());
   Field::destroy(SeaIceHeatFluxCell.label());
   Field::destroy(ShortWaveHeatFluxCell.label());
   Field::destroy(SeaIceSaltFluxCell.label());
   Field::destroy(SurfaceTracerFluxCell.label());
}
} // namespace OMEGA
