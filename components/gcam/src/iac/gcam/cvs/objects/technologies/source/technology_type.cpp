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
* \file technology_type.cpp
* \ingroup Objects
* \brief TechnologyType class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <string>
#include <cassert>

#include "util/base/include/util.h"
#include "technologies/include/technology_type.h"
#include "technologies/include/base_technology.h"
#include "investment/include/investment_utils.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"

extern Scenario* scenario;

using namespace std;

//! Constructor
TechnologyType::TechnologyType(){
}

/*! \brief Add a vintage to the technology type.
* \param aTech A pointer to a BaseTechnology.
* \return Whether the BaseTechnology was successfully added.
* \note The TechnologyType object does NOT assume ownership over the
*       technologies.
* \author Josh Lurz
*/
bool TechnologyType::addVintage( BaseTechnology* aTech ){
    // Check if a vintage already exists for the year.
    if( util::searchForValue( mVintages, aTech->getYear() ) != 0 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "A vintage already exists with year: " << aTech->getYear() 
                << " of the Technology Type: " << aTech->getName() << "." << endl;
        return false;
    }

    // Add the vintage.
    mVintages[ aTech->getYear() ] = aTech;
    return true;
}

/*! \brief Return the total capital stock up to and including a given year.
* \param aUpToYear Last year to include in capital stock summation.
* \return Total capital stock
* \author Josh Lurz
*/
double TechnologyType::getTotalCapitalStock( const int aUpToYear ) const {
    double totalStock = 0;
    for( CVintageIterator cVintage = mVintages.begin(); cVintage != mVintages.end(); ++cVintage ){
        // Check if the year of the technology is less than or equal to the top
        // year.
        if( cVintage->first <= aUpToYear ){
            totalStock += cVintage->second->getCapitalStock();
        }
    }
    return totalStock;
}

/*! \brief Initialize technologies earlier than the given year with the
*          technology from the given year.
* \details This method will attempt to merge the technology at aBaseYear
*          with any earlier technologies, it will NOT attempt to copy it
*          backwards if there were missing technologies.  This method will
*          also call initCalc on those technologies to ensure that they
*          have a chance to initialize themselves once paramaters are copied.
* \param aBaseYear Year to use as the base.
* \param aMoreSectorInfo A more sector info which will be necessary for initCalc
*                        of a technology.
* \param aRegionName The containing region's name.
* \param aSectorName The containing sector's name.
* \param aNationalAccount The national accounts object necessary for initCalc of
*                         a technology.
* \param aDemographics The demographics object necessary for initCalc of a technology.
* \param aCapitalStock The total capital stock for the sector.  I don't think this
*                      is still utilized and can probably be removed.
* \author Josh Lurz
*/
void TechnologyType::initializeTechsFromBase( const int aBaseYear, const MoreSectorInfo* aMoreSectorInfo, const string& aRegionName,
                                  const string& aSectorName, NationalAccount& aNationalAccount,
                                  const Demographic* aDemographics, const double aCapitalStock )
{
    // Find the base technology.
    const BaseTechnology* baseTech = util::searchForValue( mVintages, aBaseYear );
    if( !baseTech ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid base year. Cannot initialize previous base technologies." << endl;
        return;
    }

    // Initialize any previous vintages. 
    const int basePeriod = scenario->getModeltime()->getyr_to_per( aBaseYear );
    for( VintageIterator cVintage = mVintages.begin(); cVintage != mVintages.end(); ++cVintage ){
        // Check if the year of the technology is less than or equal to the base
        // year
        if( cVintage->first < aBaseYear ){
            cVintage->second->copyParam( baseTech, basePeriod );
            cVintage->second->initCalc( aMoreSectorInfo, aRegionName, aSectorName, aNationalAccount,
                aDemographics, aCapitalStock, basePeriod );
        }
    }
}

/*! \brief Either initialize or create a technology for the given year with an
*          existing base technology.
* \param aNewTechYear The year of the technology to initialize or create.
* \param aCurrTechYear The year of a technology to use to initialize the new
*        technology.
* \return A pointer to the newly created technology, null otherwise. 
* \note The callee must take care of the dynamically allocated memory possibly
*       returned by this function.
* \author Josh Lurz
*/
BaseTechnology* TechnologyType::initOrCreateTech( const int aNewTechYear, const int aCurrTechYear ){
    assert( mVintages.size() > 0 );

    // Find the base technology.
    const BaseTechnology* baseTech = util::searchForValue( mVintages, aCurrTechYear );
    if( !baseTech ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid current technology year. Cannot initialize current technology." << endl;
        return 0;
    }

    // Check if the technology exists. 
    BaseTechnology* newTech = util::searchForValue( mVintages, aNewTechYear );
    if( newTech != 0 ){
        const int newPeriod =
            scenario->getModeltime()->getyr_to_per( aNewTechYear );
        newTech->copyParam( baseTech, newPeriod );
        newTech = 0;
    }
    else { // Need to create the technology.
        newTech = baseTech->clone();
        newTech->setYear( aNewTechYear );
        addVintage( newTech );
    }
    return newTech;
}

/*! \brief Set total investment for the given year and a level of annual
*          investment at that year.
* \param aRegionName Region name.
* \param aPrevYear The year of the previous investment to interpolate from.
* \param aCurrentYear The year of the technology to set the investment level
*        for. 
* \param aAnnualInvestment The annual investment in this technology type at the
*        current year.
* \param aPeriod The period in which to set the investment.
* \return The total amount of annual investment that occurred.
* \author Josh Lurz
* \todo period and current year are redundant.
*/
double TechnologyType::setTotalInvestment( const string& aRegionName, const int aPrevYear, const int aCurrentYear,
                                           const double aAnnualInvestment, const int aPeriod )
{
    /*! \pre The year of the previous technology to use for interpolation is
    *        less than the new technology.
    */
    assert( aPrevYear < aCurrentYear );
    /*! \pre The annual investment passed to the technology type is a valid
    *        number and greater than zero. 
    */
    assert( util::isValidNumber( aAnnualInvestment ) );
    assert( aAnnualInvestment > 0 );

    // Find the previous technology level investment.
    /*const BaseTechnology* prevTech = util::searchForValue( mVintages, aPrevYear );
    double prevAnnualInvestment = prevTech ? prevTech->getAnnualInvestment( -1 ) : 0;*/
    
    // Find the new technology.
    BaseTechnology* currTech = util::searchForValue( mVintages, aCurrentYear );
    assert( currTech );
    
    /*
    // Set the total investment level and annual investment.
    double totalInvestment = InvestmentUtils::interpolateAndSumFlows( prevAnnualInvestment, aAnnualInvestment,
                                                                      aCurrentYear - aPrevYear );

    // Override total investment to annual times the time step.
    totalInvestment = aAnnualInvestment * scenario->getModeltime()->gettimestep( aPeriod );

    // Make sure total investment is valid. 
    assert( util::isValidNumber( totalInvestment ) );
    assert( totalInvestment > 0 );*/

    double totalInvestment = aAnnualInvestment;
    double annualInvestment = aAnnualInvestment;

    // Set the new technologies investment level and return the amount actually
    // invested.
    return currTech->setInvestment( aRegionName, annualInvestment, totalInvestment, aPeriod );
}
