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

/*!
 * \file wind_backup_calculator.cpp
 * \ingroup Objects
 * \brief WindBackupCalculator class source file.
 * \author Marshall Wise, Josh Lurz, Sonny Kim
 */

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>

#include "sectors/include/wind_backup_calculator.h"
#include "util/base/include/util.h"
#include "util/base/include/xml_helper.h"
#include "sectors/include/sector_utils.h"
#include "marketplace/include/marketplace.h"

using namespace std;
using namespace xercesc;

/*!
 * \brief Constructor.
 */
WindBackupCalculator::WindBackupCalculator()
{}

WindBackupCalculator* WindBackupCalculator::clone() const {
    return new WindBackupCalculator( *this );
}

bool WindBackupCalculator::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

const string& WindBackupCalculator::getName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& WindBackupCalculator::getXMLNameStatic() {
    const static string XML_NAME = "wind-backup-calculator";
    return XML_NAME;
}

bool WindBackupCalculator::XMLParse( const xercesc::DOMNode* node ){
    // This backup calculator does not need to parse any data.
    return true;
}

void WindBackupCalculator::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void WindBackupCalculator::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

// Documentation is inherited.
void WindBackupCalculator::initCalc( const IInfo* aTechInfo ) {
    // No information needs to be passed in
}

double WindBackupCalculator::getAverageBackupCapacity( const string& aSector,
                                                       const string& aElectricSector,
                                                       const string& aResource,
                                                       const string& aRegion,
                                                       const double aReserveMargin,
                                                       const double aAverageGridCapacityFactor,
                                                       const int aPeriod ) const
{
    // Preconditions
    assert( !aSector.empty() );
    assert( !aElectricSector.empty() );
    assert( !aResource.empty() );
    assert( !aRegion.empty() );
    assert( aReserveMargin >= 0 );
    assert( aAverageGridCapacityFactor > 0 );
    
    // This can be called from initCalc and the value may not be initialized
    // yet. Should be correct in equilibrium.
    double resourceCapacityFactor = SectorUtils::getCapacityFactor( aResource, aRegion, aPeriod );

    if( resourceCapacityFactor < util::getSmallNumber() ){
        return 0;
    }

    double backupFraction = getBackupCapacityFraction( aSector,
                                                       aElectricSector,
                                                       aResource,
                                                       aRegion,
                                                       aReserveMargin,
                                                       aAverageGridCapacityFactor,
                                                       aPeriod );
    
    // This is confusing but mathematically correct.  The backupCapacityFraction is in units of 
    // GW per GW.  The denominator (intermittent sector capacity GW)needs to be converted to energy, 
    // therefore the quotient is divided by the conversion from capacity to energy, 
    // using the resource capacity factor, which is the same as multiplying by
    // the conversion from energy to capacity as done here. Units returned here are GW/EJ.

    return backupFraction;
}

double WindBackupCalculator::getMarginalBackupCapacity( const string& aSector,
                                                        const string& aElectricSector,
                                                        const string& aResource,
                                                        const string& aRegion,
                                                        const double aReserveMargin,
                                                        const double aAverageGridCapacityFactor,
                                                        const int aPeriod ) const
{
    // Preconditions
    assert( !aSector.empty() );
    assert( !aElectricSector.empty() );
    assert( !aResource.empty() );
    assert( !aRegion.empty() );
    assert( aReserveMargin >= 0 );
    assert( aAverageGridCapacityFactor > 0 );
    
    const double EJ_PER_GWH = 3.6E-6;
    const double HOURS_PER_YEAR = 8760;
    const double UC = 1 / (EJ_PER_GWH * HOURS_PER_YEAR); // [GWe/EJ]

    double variance = SectorUtils::getVariance( aResource, aRegion, aPeriod );
    double resourceCapacityFactor = SectorUtils::getCapacityFactor( aResource, aRegion, aPeriod );
    double trialCapacityShare = SectorUtils::getTrialSupply( aRegion, aSector, aPeriod )
                              * ( aAverageGridCapacityFactor / resourceCapacityFactor );
    
    // Compute terms for Winds operating reserve due to intermittency formula
    // This is the derivative of the total backup capacity equation.
    //TODO  We need to get the documentation on this derivative, which was done
    // by Stephen Herwig for Josh Lurz. Otherwise, might simplify this and use average
    // for cost as an approximation.
    double backupCapacity = variance * trialCapacityShare * UC / aReserveMargin / resourceCapacityFactor
                          / sqrt( 1.0 + variance / pow( aReserveMargin, 2 ) * pow( trialCapacityShare, 2 ) );
    
    assert( backupCapacity >= 0 && util::isValidNumber( backupCapacity ) );
    // return value has units of GW/EJ
    return backupCapacity;
}

/*!
 * \brief Compute the backup capacity fraction.
 * \details Compute operating reserve capacity based on formula from NREL WINDS
 *          model. First compute in terms of backup capacity per total
 *          intermittent capacity, then convert to backup capacity as fraction
 *          of wind resource output in energy terms, since that is what the
 *          model and market are based on.
 * \param aSector The name of the sector which requires backup capacity.
 * \param aElectricSector The name of the electricity sector into which the
 *        sector having a backup amount calculated for will feed.
 * \param aResource The name of the resource the sector consumes.
 * \param aRegion Name of the containing region.
 * \param aReserveMargin Reserve margin for the electricity sector.
 * \param aAverageGridCapacityFactor The average electricity grid capacity
 *        factor.
 * \param aPeriod Model period.
 * \return Percent of reserve capacity per unit of intermittent capacity (e.g.,
 *         GW/GW).
 */
double WindBackupCalculator::getBackupCapacityFraction( const string& aSector,
                                                        const string& aElectricSector,
                                                        const string& aResource,
                                                        const string& aRegion,
                                                        const double aReserveMargin,
                                                        const double aAverageGridCapacityFactor,
                                                        const int aPeriod ) const
{
    // Preconditions
    assert( !aSector.empty() );
    assert( !aElectricSector.empty() );
    assert( !aResource.empty() );
    assert( !aRegion.empty() );
    assert( aReserveMargin >= 0 );
    assert( aAverageGridCapacityFactor > 0 );

    const double EJ_PER_GWH = 3.6E-6;
    const double HOURS_PER_YEAR = 8760;
    const double UC = 1 / (EJ_PER_GWH * HOURS_PER_YEAR); // [GWe/EJ]

    double variance = SectorUtils::getVariance( aResource, aRegion, aPeriod );
    double resourceCapacityFactor = SectorUtils::getCapacityFactor( aResource, aRegion, aPeriod );
    double trialCapacityShare = SectorUtils::getTrialSupply( aRegion, aSector, aPeriod )
                              * ( aAverageGridCapacityFactor / resourceCapacityFactor );

    // Compute terms for Winds operating reserve due to intermittency formula
    double backupCapacityFraction = aReserveMargin * UC / resourceCapacityFactor / trialCapacityShare *
        ( pow( ( 1.0 + variance / pow( aReserveMargin, 2 ) * pow( trialCapacityShare, 2) ), 0.5 ) - 1.0 );

    // Backup capacity fraction is positive.
    assert( util::isValidNumber( backupCapacityFraction ) && backupCapacityFraction >= 0 );
    return backupCapacityFraction;
}

double WindBackupCalculator::getReserveTotal( const string& aElectricSector,
                                              const string& aRegion,
                                              const double aReserveMargin,
                                              const double aAverageGridCapacityFactor,
                                              const int aPeriod ) const
{
    // Get current regional electricity supplied.
    const Marketplace* marketplace = scenario->getMarketplace();
    double elecSupply = marketplace->getDemand( aElectricSector, aRegion, aPeriod );

    // Electricity supply must be positive unless it is period 0 where the trial
    // supply market has not been setup. 
    assert( elecSupply >= 0 || aPeriod == 0 );

    // If there is no electricity supply there is no backup requirement.
    if( elecSupply < util::getSmallNumber() ) {
        return 0;
    }

    // Compute reserve capacity as the product of electricity reserve margin
    // and the total regional electric capacity (which is converted from
    // electric supply in Energy units to capacity units using its average
    // capacity factor).
    return aReserveMargin * SectorUtils::convertEnergyToCapacity( aAverageGridCapacityFactor, elecSupply );
}

double WindBackupCalculator::getSectorCapacity( const string& aRegionName,
                                                const string& aSectorName,
                                                const string& aDependentSectorName,
                                                const string& aResourceName,
                                                const int aPeriod ) const
{
    // Get trial amount of energy produced by the sector.
    double sectorOutputRatio = SectorUtils::getTrialSupply( aRegionName, aSectorName, aPeriod );
    double sectorSupply = sectorOutputRatio * scenario->getMarketplace()->getDemand( aDependentSectorName,
                                              aRegionName, aPeriod );

    if( sectorSupply < util::getSmallNumber() ){
        return 0;
    }

    double resourceCapacityFactor = SectorUtils::getCapacityFactor( aResourceName, aRegionName, aPeriod );
    if( resourceCapacityFactor < util::getSmallNumber() ){
        return 0;
    }

    // Using average capacity factor for resource, translate sector supply
    // in electricity or energy terms in EJ to equivalent capacity terms in
    // GW or gigawatts. Note that model's energy units need to be converted
    // to capacity units for use in the NREL WINDS operating reserve
    // computation.
    return SectorUtils::convertEnergyToCapacity( resourceCapacityFactor, sectorSupply );
}
