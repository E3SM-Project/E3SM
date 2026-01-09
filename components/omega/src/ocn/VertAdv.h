#ifndef OMEGA_VERTADV_H
#define OMEGA_VERTADV_H
//===-- ocn/VertAdv.h - vertical advection definitions ----------*- C++ -*-===//
//
/// \file
/// \brief Contains methods and variables for vertical advection
///
/// This header defines the VertAdv class which contains methods and variables
/// for calculating the vertical velocity and for vertical transport of
/// thickness, velocity, and tracers.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Error.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

// Enum for flux types
enum class VertFluxOption { Second, Third, Fourth };
// Enum for advection scheme type
enum class VertAdvOption { Standard, FCT };

class VertAdv {

 public:
   // Logicals to enable/disable advection of specific fields
   bool ThickVertAdvEnabled  = false;
   bool VelVertAdvEnabled    = false;
   bool TracerVertAdvEnabled = false;

   // enums storing options chosen by user
   VertFluxOption VertFluxChoice;
   VertAdvOption VertAdvChoice;

   // Mesh dimensions
   I4 NVertLayers;
   I4 NVertLayersP1;
   I4 NCellsOwned;
   I4 NCellsAll;
   I4 NCellsSize;
   I4 NEdgesOwned;
   I4 NEdgesAll;
   I4 NEdgesSize;

   // Number of tracers
   I4 NTracers;

   // Coefficient for blending 3rd-order and 4th-order reconstruction of
   // tracers with VertFluxOption::Third. Coef3rdOrder = 1 gives purely
   // 3rd-order, Coef3rdOrder = 0 gives purely 4th-order.
   Real Coef3rdOrder;

   Array2DReal VerticalVelocity; ///< pseudovelocity through top of cell
   Array2DReal
       TotalVerticalVelocity;    ///< transport velocity through top of Cell
   Array3DReal VertFlux;         ///< fluxes at vertical interfaces
   Array3DReal LowOrderVertFlux; ///< low-order fluxes for FCT

   // VertAdv instance name and group name
   std::string Name;
   std::string GroupName;

   // Field names
   std::string VerticalVelocityFldName;
   std::string TotalVertVelocityFldName;
   std::string VertFluxFldName;
   std::string LowOrderVertFluxFldName;

   // public methods

   // Initialize the default VertAdv instance
   static void init();

   /// Creates a new vertical advection object by calling the constructor and
   /// puts it in the AllVertAdvs map.
   static VertAdv *
   create(const std::string &Name, ///< [in] name for new VertAdv
          const HorzMesh *Mesh,    ///< [in] associated HorzMesh
          const VertCoord *VCoord, ///< [in] associated VertCoord
          Config *Options          ///< [in] configuration options
   );

   /// Destructor - deallocates all memory and deletes a VertAdv
   ~VertAdv();

   /// Deallocates arrays
   static void clear();

   /// Remove a VertAdv by name
   static void erase(std::string InName);

   /// Retrieve the default VertAdv
   static VertAdv *getDefault();

   /// Retreive a VertAdv by name
   static VertAdv *get(std::string name);

   /// Read and set config options
   void readConfigOptions(Config *Options);

   /// Determine transport due to vertical advection from divergence of
   /// horizontal advection.
   void computeVerticalVelocity(
       const Array2DReal &NormalVelocity,    ///< [in] horizontal velocity
       const Array2DReal &FluxLayerThickEdge ///< [in] layer thickness at edges
   );

   /// Compute pseudo thickness tendency due to vertical advection
   void computeThicknessVAdvTend(
       const Array2DReal &ThickTend ///< [inout] thickness tendency
   );

   /// Compute velocity tendency due to vertical advection
   void computeVelocityVAdvTend(
       const Array2DReal &VelTend, ///< [inout] horizontal velocity tendency
       const Array2DReal &NormalVelocity,    ///< [in] horizontal velocity
       const Array2DReal &FluxLayerThickEdge ///< [in] layer thickness at edges
   );

   /// Compute tracer tendency due to vertical advection, The LayerThickness at
   /// the beginning of the time step is passed for standard advection, while
   /// the provisional thickness including horizontal thickness flux is passed
   /// for flux-corrected transport. The TimeStep is only needed as an argument
   /// for flux-corrected transport
   void computeTracerVAdvTend(
       const Array3DReal &TracerTend,     ///< [inout] tracer tendencies
       const Array3DReal &Tracers,        ///< [in] tracer array
       const Array2DReal &LayerThickness, ///< [in] layer thickness
       const TimeInterval TimeStep =
           TimeInterval() ///< [in] (optional) time step
   );

   /// Compute tracer fluxes due to vertical advection
   void computeVerticalFluxes(
       const Array3DReal &Tracers,       ///< [in] tracer array
       const Array2DReal &LayerThickness ///< [in] layer thickness
   );

   /// Compute tracer tendencies due to vertical advection using standard
   /// advection scheme
   void computeStdVAdvTend(
       const Array3DReal &TracerTend ///< [inout] tracer tendencies
   );

   /// Compute tracer tendencies due to vertical advection using flux-corrected
   /// transport scheme
   void computeFCTVAdvTend(
       const Array3DReal &TracerTend,    ///< [inout] tracer tendencies
       const Array3DReal &Tracers,       ///< [in] tracer array
       const Array2DReal &ProvThickness, ///< [in] provisional layer thickness
       const Real Dt                     ///< [in] time step
   );

 private:
   // Vertical loop bounds from VertCoord
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;

   // Arrays from HorzMesh
   Array2DI4 CellsOnEdge;
   Array2DI4 EdgesOnCell;
   Array1DI4 NEdgesOnCell;
   Array1DReal AreaCell;
   Array1DReal DvEdge;
   Array2DReal EdgeSignOnCell;

   // small number to avoid numerical difficulties in computing
   // scale factors for FCT.
   Real Eps = 1e-10;

   const HorzMesh *Mesh;
   const VertCoord *VCoord;

   static VertAdv *DefaultVertAdv;
   static std::map<std::string, std::unique_ptr<VertAdv>> AllVertAdvs;

   // private methods

   /// construct a new vertical coordinate object
   VertAdv(const std::string &Name, ///< [in] Name for new VertAdv
           const HorzMesh *Mesh,    ///< [in] associated HorzMesh
           const VertCoord *VCoord, ///< [in] associated VertCoord
           Config *Options          ///< [in] configuration options
   );

   /// define field metadata
   void defineFields();

   // Forbid copy and move construction
   VertAdv(const VertAdv &) = delete;
   VertAdv(VertAdv &&)      = delete;

}; // end class VertAdv

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_VERTICALADV_H
