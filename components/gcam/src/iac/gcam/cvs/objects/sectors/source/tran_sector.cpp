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
* \file tran_sector.cpp
* \ingroup Objects
* \brief transporation technology class source file.
* \author Marshall Wise, Sonny Kim, Josh Lurz
* \date $Date: 2006/02/21 15:32:03 $
* \version $Revision: 1.20.2.1 $
*/

#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <cassert>

#include "marketplace/include/marketplace.h"
#include "util/base/include/model_time.h"
#include "containers/include/scenario.h"
#include "sectors/include/tran_sector.h"
#include "containers/include/gdp.h"

// xml headers
#include "util/base/include/xml_helper.h"
#include <xercesc/dom/DOMNode.hpp>

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default constructor
TranSector::TranSector() {
    mPercentLicensed.resize( scenario->getModeltime()->getmaxper() );
}

void TranSector::completeInit( const string& aRegionName, const IInfo* aRegionInfo ) {
    EnergyFinalDemand::completeInit( aRegionName, aRegionInfo );
}

/*! \brief Initialize the TranSector.
* \param aRegionName Region name.
* \param aPeriod Period for which to initialize the TranSector.
*/
void TranSector::initCalc( const string& aRegionName,
                           const GDP* aGDP,
                           const Demographic* aDemographics,
                           const int aPeriod )
{

    // Set the energy demand for period 1 which is fixed.
    if( aPeriod == 0 || aPeriod == 1 ){
        mServiceDemands[ aPeriod ] = mBaseService[ 0 ];

        // Calculate the base scalers.
        calcBaseScalers( aGDP, aPeriod );
    }

    EnergyFinalDemand::initCalc( aRegionName, aGDP, aDemographics, aPeriod );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& TranSector::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& TranSector::getXMLNameStatic() {
    const static string XML_NAME = "tranSector";
    return XML_NAME;
}

//! Parses any input variables specific to derived classes
bool TranSector::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    if( nodeName == "percentLicensed" ) {
        XMLHelper<double>::insertValueIntoVector( curr, mPercentLicensed, scenario->getModeltime() );
    }
    // TODO: Change the name of this variable.
    else if( nodeName == "aeei" ){
        // TODO: Fix this with merge.
        XMLHelper<double>::insertValueIntoVector( curr, mTechChange,
                                                  scenario->getModeltime() );
    }
    else {
        return false;
    }
    return true;
}

/*! \brief XML output stream for derived classes
*
* Function writes output due to any variables specific to derived classes to XML
* \author Josh Lurz
* \param out reference to the output stream
* \param tabs A tabs object responsible for printing the correct number of tabs. 
*/
void TranSector::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    XMLWriteVector( mPercentLicensed, "percentLicensed", out, tabs, modeltime, 0.0 );
    // TODO: Fix this with merge.
    XMLWriteVector( mTechChange, "aeei", out, tabs, modeltime );
}

//! Write object to debugging xml output stream.
void TranSector::toDebugXMLDerived( const int aPeriod, ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mPercentLicensed[ aPeriod ], "percentLicensed", out, tabs );
    XMLWriteElement( mTechChange[ aPeriod ], "aeei", out, tabs );
}


/*!
* \brief Calculate the base scalers for the demand function.
* \details
* \param 
*/
void TranSector::calcBaseScalers( const GDP* aGDP, const int aPeriod ){
        
    // Calculate base year scalers. The GDP is known in the base period and not
    // adjusted for energy prices, so the GDP ratio cannot change throughout the
    // iteration. The price ratio is assumed to be 1 and so does not effect the
    // equation.
    // TODO: Convert to demand function.
    if ( mDemandFunction->isPerCapitaBased() ) {
        double scaledGdpPerCapita = aGDP->getBestScaledGDPperCap( aPeriod ); 
        mBaseScaler = mServiceDemands[ 0 ] * mPercentLicensed[ aPeriod ] 
        * pow( scaledGdpPerCapita, -1 * mIncomeElasticity[ aPeriod ] );

        mBaseScalerNotLic = mServiceDemands[ 0 ] * ( 1 - mPercentLicensed[ aPeriod ] )
            * pow( scaledGdpPerCapita,-1 * mIncomeElasticity[ aPeriod ] );
    }
    else {
        // Note that the percent licensed is not used here. This is a standard
        // demand function.
        double gdpRatio = aGDP->getApproxScaledGDP( aPeriod );
        mBaseScaler = mServiceDemands[ 0 ] * pow( gdpRatio, -1 *
                      mIncomeElasticity[ aPeriod ] );
    }
}

//! Aggrgate sector energy service demand function.
void TranSector::setFinalDemand( const string& aRegionName,
                                 const Demographic* aDemographics,
                                 const GDP* gdp,
                                 const int aPeriod )
{ 
    double scaledGdpPerCapita = gdp->getBestScaledGDPperCap(aPeriod); 
             
    double gdp1 = gdp->getApproxScaledGDP(aPeriod);
    
    // Data for periods 0 and 1 is read-in, so do not recalculate. The known
    // demands will be added to the marketplace.
    if ( aPeriod > 1 ){
        // note normalized to previous year not base year
        // has implications for how technical change is applied.

        // TODO: Sector utils function.
        const Marketplace* marketplace = scenario->getMarketplace();
        double sectorPrice = marketplace->getPrice( mName, aRegionName, aPeriod );
        double prevSectorPrice = marketplace->getPrice( mName, aRegionName, aPeriod - 1 );
        assert( prevSectorPrice > 0 );
        double priceRatio  = sectorPrice / prevSectorPrice;

        double priceRatioNotLic = priceRatio;
        
        double serviceDemand;
        // TODO: Convert to demand function.
        if( mDemandFunction->isPerCapitaBased() ){ // demand based on per capita GDP
            serviceDemand = mBaseScaler * pow( priceRatio, mPriceElasticity[ aPeriod ] )
                * pow( scaledGdpPerCapita, mIncomeElasticity[ aPeriod ] )
                + mBaseScalerNotLic * pow( priceRatioNotLic, mPriceElasticity[ aPeriod ] )
                * pow( scaledGdpPerCapita, mIncomeElasticity[ aPeriod ] );

            // need to multiply above by population ratio (current population/base year
            // population).  The gdp ratio provides the population ratio.
            // TODO: Fix this on GDP merge.
            serviceDemand *= gdp1 / scaledGdpPerCapita;
        }
        else { // demand based on scale of GDP
            serviceDemand = mBaseScaler * pow( priceRatio, mPriceElasticity[ aPeriod ] ) 
                * pow( gdp1, mIncomeElasticity[ aPeriod ] );
        }

        // adjust demand for AEEI, autonomous end-use energy intensity
        // note: not using cummulative technical change
        mServiceDemands[ aPeriod ] = serviceDemand;
    }

    // Add the demand to the marketplace.
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->addToDemand( mName, aRegionName, mServiceDemands[ aPeriod ], aPeriod );
}

