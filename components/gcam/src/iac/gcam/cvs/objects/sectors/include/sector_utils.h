/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 * 
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 * 
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community 
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 * 
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 * 
 */

#ifndef _SECTOR_UTILS_H_
#define _SECTOR_UTILS_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file sector_utils.h
 * \ingroup Objects
 * \brief The SectorUtils class header file.
 * \author Josh Lurz, Sonny Kim
 */

#include <string>
#include <vector>
#include "util/base/include/hash_map.h"
#include "util/base/include/value.h"
#include "util/base/include/time_vector.h"

class IInfo;

/*! 
 * \ingroup Objects
 * \brief This class contains a set of static helper methods which the various
 *        sectors and subsectors use.
 * \author Josh Lurz
 */
class SectorUtils {
public:
    static bool createTrialSupplyMarket( const std::string& aRegionName,
                                         const std::string& aSectorName,
                                         const IInfo* aTechnologyInfo);

    static void addToTrialDemand( const std::string& aRegionName,
                                  const std::string& aSectorName,
                                  const double aSupply,
                                  const int aPeriod );

    static double getTrialSupply( const std::string& aRegionName,
                                  const std::string& aSectorName,
                                  const int aPeriod );

    static void askToCreateTrialSupply( const std::string& aRegionName,
                                        const std::string& aSectorName );

    static double calcFixedOutputScaleFactor( const double aMarketDemand,
                                              const double aFixedOutput );

    static double normalizeShares( std::vector<double>& aShares );

    static double calcPriceRatio( const std::string& aRegionName,
                                  const std::string& aSectorName,
                                  const int aBasePeriod,
                                  const int aCurrentPeriod );

    static const std::string createTFEMarketName( const std::string& aSectorName );

    static void setFinalEnergyFlag( const std::string& aRegionName,
                                    const std::string& aSectorName );

    static bool isFinalEnergySector( const std::string& aRegionName,
                                     const std::string& aSectorName );

    static double getVariance( const std::string& aResourceName,
                               const std::string& aRegionName,
                               const int aPeriod );

    static int getDemandNormPeriod( const int aPeriod );

    static double getCapacityFactor( const std::string& aResourceName,
                                     const std::string& aRegionName,
                                     const int aPeriod );
    
    static double convertEnergyToCapacity( const double aCapacityFactor,
                                           const double aEnergy );

    static double convertCapacityToEnergy( const double aCapacityFactor,
                                           const double aCapacity );

    static const std::string getTrialMarketName( const std::string& aSectorName );

     /*! \brief Fills missing period elements in a Value vector with values linearly 
     *        interpolated from initialized or read-in values.
     * \detail This method is intended for enabling variable time-step capability
     *         and filling in values that have not been read-in or initialized
     *         in a PeriodVector.
     * \param aValueVector a period vector of Values.
     */
    static void fillMissingPeriodVectorInterpolated( objects::PeriodVector<Value>& aPeriodVector );

     /*! \brief Fills missing period elements in a Value vector with values 
     *         available from the next initialized or read-in values.
     * \detail This method is intended for enabling variable time-step capability
     *         and filling in values that have not been read-in or initialized
     *         in a PeriodVector.
     * \param aValueVector a period vector of Values.
     */
    static void fillMissingPeriodVectorNextAvailable( objects::PeriodVector<Value>& aPeriodVector );

protected:

    static HashMap<std::string, std::string> sTrialMarketNames;
};

#endif // _SECTOR_UTILS_H_
