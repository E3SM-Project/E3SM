#ifndef OMEGA_MOC_H
#define OMEGA_MOC_H

//===-- analysis/analysisGroups/MOC.h - MOC AnalysisGroup -------*- C++ -*-===//
//
/// \file
/// \brief Defines the MOC analysis group for computing Meridional Overturning
/// Circulation
///
/// MOC is a bundled AnalysisGroup that computes the Meridional Overturning
/// Circulation (MOC) streamfunction using a latitude-binned method. Each named
/// region produces a streamfunction output field. An optional list of transects
/// (index-paired with Regions) provides boundary conditions that anchor the
/// streamfunction to the transport through a southern boundary transect.
///
/// Latitude-binned MOC (for regions):
/// 1. Converting vertical pseudo-velocity to geometric coordinates
/// 2. Computing vertical flux (VerticalVelocity × AreaCell)
/// 3. Optionally applying a regional mask (ExtractRegion)
/// 4. Accumulating flux into latitude bins
/// 5. Performing horizontal integration (cumulative sum from south to north),
///    optionally seeded with a boundary transect transport (BC)
/// 6. Converting to Sverdrups (Sv = 1e6 m³/s)
///
/// Transect BC chain (internal, not written to output):
/// 1. Converting MeanPseudoThickEdge to geometric layer thickness
/// 2. Computing edge transport (LayerThickness × NormalVelocity × DvEdge)
/// 3. Accumulating transport across transect edges
/// 4. Performing vertical integration (cumulative sum from bottom to top)
///
/// Configuration options:
/// - NumBins: Number of latitude bins (default: 180)
/// - MinLat: Minimum latitude in degrees (default: -90.0)
/// - MaxLat: Maximum latitude in degrees (default: 90.0)
/// - Regions: List of region names for regional MOC (optional)
/// - Transects: List of transect names (index-paired with non-Global Regions;
///              required when any non-Global region is specified)
/// - ReductionPeriod: Time-averaging period (e.g., "1month")
/// - SnapshotPeriod: Instantaneous output period (e.g., "1month")
///
/// Example usage in config:
/// \code
/// MOC:
///   NumBins: 180
///   MinLat: -90.0
///   MaxLat: 90.0
///   Regions: [Global, REGION1, REGION2]
///   Transects: [TRANSECT1, TRANSECT2]
///   ReductionPeriod: [1Month]
///   SnapshotPeriod: [1Day]
/// \endcode
///
/// The MOC streamfunction represents the zonally integrated meridional mass
/// transport as a function of latitude and depth. Units are Sverdrups (Sv),
/// where 1 Sv = 10^6 m³/s.
///
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AnalysisGroup.h"
#include "Config.h"
#include "operators/Ops.h"
#include <string>

namespace OMEGA {

/// MOC is a bundled AnalysisGroup that computes the Meridional Overturning
/// Circulation (MOC) streamfunction using a latitude-binned method. Converts
/// vertical pseudo-velocity to geometric coordinates, computes vertical flux
/// (velocity × area), accumulates flux into latitude bins, and integrates
/// horizontally (south to north) to produce the streamfunction, optionally
/// seeded by a transect boundary condition. Converts results to Sverdrups and
/// supports optional temporal averaging.
class MOC : public AnalysisGroup {
 public:
   /// Constructs a MOC analysis group. Reads configuration options (number
   /// of bins, latitude range, regions, output frequency), builds operator
   /// chains for MOC computation, creates IOStreams for output, and
   /// associates operators with streams.
   MOC(const std::string &GroupName, ///< [in] name of this group
       Config &AnalysisGroupOptions, ///< [in] group configuration
       Analysis *AnalysisManager     ///< [in] analysis manager
   );

   /// Default destructor
   ~MOC() = default;

 private:
   /// Builds the complete operator chain string for a single MOC region.
   /// Constructs the chain string encoding all transformation steps:
   /// 1. PseudoToGeometric: Convert vertical pseudo-velocity to geometric
   /// 2. BinaryMultiply: Compute vertical flux (velocity × area)
   /// 3. ExtractRegion: Apply regional mask (optional)
   /// 4. BinnedAccumulator: Accumulate flux into latitude bins
   /// 5. PrefixSum: Horizontal integration (south to north across bins)
   /// 6. ScalarMultiply: Convert to Sverdrups
   ///
   /// The chain string is later parsed by buildTemporalChains() which
   /// appends temporal operators and calls parseChainAndBuildOps().
   ///
   /// Returns the complete chain string (stem) for temporal operator appending.
   std::string buildMOCChain(
       const std::string &RegionName, ///< [in] region name (e.g., "Global")
       const std::string
           &RegionMaskName, ///< [in] mask field name (empty = global)
       I4 NumBins,          ///< [in] number of latitude bins
       Real MinLat,         ///< [in] minimum latitude in degrees
       Real MaxLat,         ///< [in] maximum latitude in degrees
       const std::string &BCFieldName, ///< [in] BC field name (empty = no BC)
       Analysis *AnalysisManager       ///< [in] analysis manager
   );

   /// Builds an internal transect BC chain. Steps:
   /// 1. PseudoToGeometric: Convert MeanPseudoThickEdge to geometric layer
   ///    thickness
   /// 2. BinaryMultiply: Multiply by NormalVelocity (thickness × velocity)
   /// 3. BinaryMultiply: Multiply by DvEdge (transport × edge width)
   /// 4. TransectAccumulator: Accumulate across transect edges
   /// 5. PrefixSum: Vertical integration (bottom to top)
   /// Registers ops immediately via parseChainAndBuildOps.
   /// Returns the PrefixSum output field name for use as BC.
   std::string
   buildTransectBCChain(const std::string &TransectName, ///< [in] transect name
                        Analysis *AnalysisManager ///< [in] analysis manager
   );

   /// Configuration storage for each MOC chain. Vector of Config objects
   /// containing operator-specific parameters (NumBins, Dimension, Reverse,
   /// TransectName). Indexed in the same order as the combined region and
   /// transect chain stems, these configs are passed to buildTemporalChains().
   std::vector<Config> ChainConfigs;

   /// Name of the static latitude bin boundary field (size NumBins+1).
   /// Created once and added to all MOC output streams.
   std::string BinBoundaryFieldName;

   /// Names of the per-region area-weighted representative interface depth
   /// coordinate fields (size NVertLayersP1). One entry per region; each is
   /// added to the MOC output streams. Regional entries use the region's mask
   /// via the upstream ExtractRegion.
   std::vector<std::string> DepthFieldNames;

   /// Builds the per-region depth coordinate chain (GeomZInterface averaged
   /// horizontally, area-weighted, optionally masked to the region) and returns
   /// the resulting output field name.
   std::string buildDepthCoordChain(
       const std::string &RegionName, ///< [in] region name (e.g., "Global")
       const std::string
           &RegionMaskName,      ///< [in] mask field name (empty = global)
       Analysis *AnalysisManager ///< [in] analysis manager
   );

}; // end class MOC

} // end namespace OMEGA

#endif
