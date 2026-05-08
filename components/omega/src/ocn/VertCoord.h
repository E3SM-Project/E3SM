#ifndef OMEGA_VERTCOORD_H
#define OMEGA_VERTCOORD_H
//===-- base/VertCoord.h - vertical coordinate definitions ------*- C++ -*-===//
//
/// \file
/// \brief Contains the vertical mesh variables and methods for Omega
///
/// This header defines the VertCoord class which contains information related
/// to the vertical coordinate, the vertical extent of the mesh, and the active
/// vertical layers for each ocean column.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Error.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

#include <memory>
#include <string>

namespace OMEGA {

enum class MovementWeightType {
   Fixed,  /// Distribute perturbations in top level
   Uniform /// Uniform stretching
};

KOKKOS_INLINE_FUNCTION int vertRange(int KMin, int KMax) {
   return KMax - KMin + 1;
}

KOKKOS_INLINE_FUNCTION int vertRangeChunked(int KMin, int KMax) {
   return (vertRange(KMin, KMax) + VecLength - 1) / VecLength;
}

KOKKOS_INLINE_FUNCTION int chunkStart(int KChunk, int KMin) {
   return KMin + KChunk * VecLength;
}

KOKKOS_INLINE_FUNCTION int chunkLength(int KChunk, int KStart, int KMax) {
   if constexpr (VecLength == 1) {
      return 1;
   } else {
      return (KStart + VecLength - 1) > KMax ? (KMax - KStart + 1) : VecLength;
   }
}

class VertCoord {

 private:
   // Variables from Decomp
   I4 NCellsOwned;
   I4 NCellsAll;
   I4 NCellsSize;
   I4 NEdgesOwned;
   I4 NEdgesHalo0;
   I4 NEdgesAll;
   I4 NEdgesSize;
   I4 NVerticesOwned;
   I4 NVerticesHalo0;
   I4 NVerticesAll;
   I4 NVerticesSize;
   I4 VertexDegree;
   Array2DI4 CellsOnEdge;
   Array2DI4 CellsOnVertex;

   // Choice of VertCoorMovementWeight type
   MovementWeightType MvmtWgtChoice;

   static VertCoord *DefaultVertCoord;
   static std::map<std::string, std::unique_ptr<VertCoord>> AllVertCoords;

   // methods

   /// construct a new vertical coordinate object
   VertCoord(const std::string &Name,  ///< [in] Name for new VertCoord
             const Decomp *MeshDecomp, ///< [in] associated Decomp
             Halo *MeshHalo,           ///< [in] mesh halo exchanger
             Config *Options,          ///< [in] configuration options
             const bool ReadStream,    ///< [in] logical to read stream
             const int NVertLayers     ///< [in] int to set vertical dim
   );

   /// define field metadata
   void defineFields();

   // Forbid copy and move construction
   VertCoord(const VertCoord &) = delete;
   VertCoord(VertCoord &&)      = delete;

 public:
   // Vertical dimension
   I4 NVertLayers;
   I4 NVertLayersP1;

   // Variables computed
   Array2DReal PressureInterface;
   Array2DReal PressureMid;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array2DReal GeopotentialMid;
   Array2DReal LayerThicknessTarget;

   HostArray2DReal PressureInterfaceH;
   HostArray2DReal PressureMidH;
   HostArray2DReal ZInterfaceH;
   HostArray2DReal ZMidH;
   HostArray2DReal GeopotentialMidH;
   HostArray2DReal LayerThicknessTargetH;

   // Vertical loop bounds
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeTop;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeBot;
   Array1DI4 MinLayerVertexTop;
   Array1DI4 MaxLayerVertexTop;
   Array1DI4 MinLayerVertexBot;
   Array1DI4 MaxLayerVertexBot;

   HostArray1DI4 MinLayerCellH;
   HostArray1DI4 MaxLayerCellH;
   HostArray1DI4 MinLayerEdgeTopH;
   HostArray1DI4 MaxLayerEdgeTopH;
   HostArray1DI4 MinLayerEdgeBotH;
   HostArray1DI4 MaxLayerEdgeBotH;
   HostArray1DI4 MinLayerVertexTopH;
   HostArray1DI4 MaxLayerVertexTopH;
   HostArray1DI4 MinLayerVertexBotH;
   HostArray1DI4 MaxLayerVertexBotH;

   // p star coordinate variables
   Array1DReal VertCoordMovementWeights;
   Array2DReal RefPseudoThickness;

   HostArray1DReal VertCoordMovementWeightsH;
   HostArray2DReal RefPseudoThicknessH;

   // Masks
   Array2DReal EdgeMask;        ///< Mask to determine if computations should be
                                ///  done on edge
   HostArray2DReal EdgeMaskH;   ///< Mask to determine if computations should be
                                ///  done on edge
   Array2DReal CellMask;        ///< Mask to determine if computations should be
                                ///  done on cell
   HostArray2DReal CellMaskH;   ///< Mask to determine if computations should be
                                ///  done on cell
   Array2DReal VertexMask;      ///< Mask to determine if computations should be
                                ///  done on vertex
   HostArray2DReal VertexMaskH; ///< Mask to determine if computations should be
                                ///  done on vertex

   // BottomDepth read from mesh file
   Array1DReal BottomDepth;

   HostArray1DReal BottomDepthH;

   // TODO: Temporary handling of SurfacePressure
   Array1DReal SurfacePressure;

   HostArray1DReal SurfacePressureH;

   // VertCoord instance name and FieldGroup names
   std::string Name;
   std::string InitGroupName;
   std::string GroupName;

   // Field names
   std::string MinLayerCellFldName;   ///< Field name for MinLayerCell
   std::string MaxLayerCellFldName;   ///< Field name for MaxLayerCell
   std::string BottomDepthFldName;    ///< Field name for BottomDepth
   std::string RefPseudoThickFldName; ///< Field name for RefPseudoThickness
   std::string
       VCoordMvmtWgtsFldName;      ///< Field name for VertCoordMovementWeights
   std::string PressInterfFldName; ///< Field name for interface pressure
   std::string PressMidFldName;    ///< Field name for midpoint pressure
   std::string ZInterfFldName;     ///< Field name for interface Z height
   std::string ZMidFldName;        ///< Field name for midpoint Z height
   std::string GeopotFldName;      ///< Field name for geopotential
   std::string LyrThickTargetFldName; ///< Field name for target thickness

   // methods

   /// Initialization of default vertical coordinate
   static void init(const bool ReadStream = true, const int NVertLayers = 0);

   /// Creates a new vertical coordinate object by calling the constructor and
   /// puts it in the AllVertCoords map.
   static VertCoord *
   create(const std::string &Name,      /// [in] name for new VertCoord
          const Decomp *MeshDecomp,     /// [in] associated Decomp
          Halo *MeshHalo,               /// [in] mesh halo exchanger
          Config *Options,              /// [in] configuration options
          const bool ReadStream = true, /// [in] optional logical to read stream
          const int NVertLayers = 0 /// [in] optional int to set vertical dim
   );

   /// Copy member arrays from device to host
   void copyToHost();

   /// Copy member arrays from host to device
   void copyToDevice();

   /// Destructor - deallocates all memory and deletes a VertCoord
   ~VertCoord();

   /// Deallocates arrays
   static void clear();

   /// Remove a VertCoord by name
   static void erase(std::string InName);

   /// Retrieve the default VertCoord
   static VertCoord *getDefault();

   /// Retreive a VertCoord by name
   static VertCoord *get(std::string name);

   // Array initialization methods
   void setStreamArrays(const bool ReadStream, Halo *MeshHalo);
   void minMaxLayerEdge(Halo *MeshHalo);
   void minMaxLayerVertex(Halo *MeshHalo);
   void initMovementWeights();

   /// Initialize computational masks
   void setMasks();

   /// Sums the mass thickness times g from the top layer down, starting with
   /// the surface pressure
   void
   computePressure(const Array2DReal &LayerThickness, ///< [in] pseudo thickness
                   const Array1DReal &SurfacePressure ///< [in] surface pressure
   );

   /// Sum the mass thickness times specific volume from the bottom layer up,
   /// starting with the bottom elevation
   void
   computeZHeight(const Array2DReal &LayerThickness, ///< [in] pseudo thickness
                  const Array2DReal &SpecVol);       ///< [in] specific volume

   /// Sum the z height times g, the tidal potential, and self attraction and
   /// loading
   void computeGeopotential(
       const Array1DReal &TidalPotential, ///< [in] tidal potential
       const Array1DReal
           &SelfAttractionLoading ///< [in] self attraction and loading
   );

   /// Determine mass thickness used for the target vertical coordinate
   void computeTargetThickness();

}; // end class VertCoord

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_VERTCOORD_H
