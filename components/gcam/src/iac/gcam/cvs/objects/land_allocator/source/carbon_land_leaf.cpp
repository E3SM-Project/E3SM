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
 * \file carbon_land_leaf.cpp
 * \ingroup Objects
 * \brief CarbonLandLeaf class source file.
 * \author Kate Calvin
 */

#include "util/base/include/definitions.h"
#include "land_allocator/include/carbon_land_leaf.h"
#include "land_allocator/include/land_use_history.h"
#include "util/base/include/xml_helper.h"
#include "ccarbon_model/include/land_carbon_densities.h"
#include "util/base/include/summary.h"
#include "util/base/include/ivisitor.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/iinfo.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

/*!
 * \brief Default constructor.
 * \param aParent Pointer to this leafs's parent.
 * \author James Blackwood
*/
CarbonLandLeaf::CarbonLandLeaf( const ALandAllocatorItem* aParent ):
// Default the name to the empty string. It will be read in during XML parsing.
LandLeaf( aParent, "" )
{
}

//! Destructor
CarbonLandLeaf::~CarbonLandLeaf() {
}

bool CarbonLandLeaf::XMLDerivedClassParse( const std::string& aNodeName,
                                              const xercesc::DOMNode* aCurr )
{
    return true;
}

void CarbonLandLeaf::toInputXML( ostream& aOut, Tabs* aTabs ) const {

}

const string& CarbonLandLeaf::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \return The constant XML_NAME as a static.
*/
const string& CarbonLandLeaf::getXMLNameStatic() {
    const static string XML_NAME = "CarbonLandLeaf";
    return XML_NAME;
}

/*!
* \brief Sets a the profit rate of a land leaf
* \details This method adjusts the profit rate of an unmanaged land leaf
*          to account for the carbon value of land if the ag subsidy is
*          is active and a carbon price exists. 
* \param aRegionName Region.
* \param aPeriod Period.
*/
void CarbonLandLeaf::setUnmanagedLandProfitRate( const string& aRegionName,
                                                    double aAverageProfitRate,
                                                    const int aPeriod )
{
    bool agSubsidy = Configuration::getInstance()->getBool( "agSubsidy", true );
    if ( agSubsidy ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Ag subsidy should not be used in combination with ";
        mainLog << "carbon plantations." << endl;
    }

    double profitRate = 0.0;

    // The base profit rate is based on the carbon density and the carbon price
    
    // Check if a carbon market exists and has a non-zero price.
    const Marketplace* marketplace = scenario->getMarketplace();
    double carbonPrice = marketplace->getPrice( "CO2", aRegionName, aPeriod, false );

    // Retrieve proportional tax rate.
    const IInfo* marketInfo = marketplace->getMarketInfo( "CO2", aRegionName, aPeriod, false );
    // Note: the key includes the region name.
    const double proportionalTaxRate = 
            ( marketInfo && marketInfo->hasValue( "proportional-tax-rate" + aRegionName ) ) 
            ? marketInfo->getDouble( "proportional-tax-rate" + aRegionName, true )
            : 1.0;

    // If a carbon price exists, calculate the subsidy
    if( carbonPrice != Marketplace::NO_MARKET_PRICE && carbonPrice > 0.0 ){
        // Carbon price is in 1990$, but land value is in 1975$ so we need to convert
        const double dollar_conversion_75_90 = 2.212;
        carbonPrice /= dollar_conversion_75_90;

        // Adjust carbon price with the proportional tax rate.
        carbonPrice *= proportionalTaxRate;

        // With carbon content in Tg C/KHa, convert to $/KHa.
        const double tC_in_TgC = 1000000.0;

        // We are only subsidizing for carbon contents above the read in minimum
        const int year = scenario->getModeltime()->getper_to_yr( aPeriod );
        double incrementalAboveCDensity = mCarbonContentCalc->getActualAboveGroundCarbonDensity( year )
            - mMinAboveGroundCDensity;
        double incrementalBelowCDensity = mCarbonContentCalc->getActualBelowGroundCarbonDensity( year )
            - mMinBelowGroundCDensity;

        // Calculate the carbon value as the total carbon content of the land
        // multiplied by the carbon price and the interest rate.
        profitRate = ( incrementalAboveCDensity * mCarbonContentCalc->getAboveGroundCarbonSubsidyDiscountFactor()
            + incrementalBelowCDensity * mCarbonContentCalc->getBelowGroundCarbonSubsidyDiscountFactor() )
            * carbonPrice * ( mInterestRate - mCarbonPriceIncreaseRate[ aPeriod ] )* tC_in_TgC;
    }

    mProfitRate[ aPeriod ] = profitRate;
}  


/*!
* \brief Sets a the profit rate of a land leaf
* \param aRegionName Region.
* \param aProductName Name of land leaf
* \param aProfitRate Profit rate
* \param aPeriod Period.
*/
void CarbonLandLeaf::setProfitRate( const string& aRegionName,
                                 const string& aProductName,
                                 const double aProfitRate,
                                 const int aPeriod )
{
    // This shouldn't do anything for unmanaged land leafs
}
