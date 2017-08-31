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
* \file consumer.cpp
* \ingroup Objects-SGM
* \brief The Consumer class source file.
*
*  Detailed Description.
*
* \author Pralit Patel
* \author Sonny Kim
*/
#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>

#include "consumers/include/consumer.h"
#include "util/base/include/xml_helper.h"
#include "functions/include/iinput.h"
#include "functions/include/ifunction.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "functions/include/function_manager.h"
#include "util/base/include/ivisitor.h"
#include "technologies/include/primary_output.h"
#include "emissions/include/aghg.h"
#include "functions/include/node_input.h"
#include "functions/include/function_utils.h" // TODO: can remove once initCalc stuff is sorted out
#include "containers/include/national_account.h"
#include "containers/include/iinfo.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

typedef vector<AGHG*>::const_iterator CGHGIterator;
typedef vector<AGHG*>::iterator GHGIterator;

//! Default Constructor
Consumer::Consumer():
mUtilityParameterA( 1.0 ) {
}

void Consumer::initCalc( const MoreSectorInfo* aMoreSectorInfo,
                         const string& aRegionName, 
                         const string& aSectorName,
                         NationalAccount& nationalAccount,
                         const Demographic* aDemographics,
                         const double aCapitalStock,
                         const int aPeriod )
{
    const Modeltime* modeltime = scenario->getModeltime();
    const bool isInitialYear = modeltime->getper_to_yr( aPeriod ) == year;
    if( isInitialYear ) {
        mNestedInputRoot->initCalc( aRegionName, aSectorName, isNewInvestment( aPeriod ), isTrade(), aPeriod );
        if( aPeriod == 0 ) {
            mNestedInputRoot->initialize();
            // consumers do not use alpha zero
            mNestedInputRoot->setCoefficient( 0, aPeriod );
        }
    }
    BaseTechnology::initCalc( aMoreSectorInfo, aRegionName, aSectorName,
                              nationalAccount, aDemographics, aCapitalStock,
                              aPeriod );
    mLeafInputs = FunctionUtils::getLeafInputs( mNestedInputRoot );
}

//! Calculate Demand
void Consumer::calcInputDemand( double aConsumption, const string& aRegionName, 
                                const string& aSectorName, int aPeriod ) 
{
    double output = mNestedInputRoot->calcInputDemand( aRegionName, aSectorName, aPeriod, 
                aConsumption, mUtilityParameterA, mNestedInputRoot->getCoefficient( aPeriod ) ); // alpha zero is the root's alpha
    // this is likely a currency output but for consumers it does not matter
    mOutputs[ 0 ]->setPhysicalOutput( output, aRegionName, 0, aPeriod );

}

void Consumer::updateMarketplace( const string& aSectorName, const string& aRegionName, const int aPeriod ) {
    // need to create the list here so that marketplaces get set up correctly
    // TODO: I could just create an updateMarketplace in the node input
    Marketplace* marketplace = scenario->getMarketplace();
    for( unsigned int i = 0; i < mLeafInputs.size(); i++ ) {
        // don't add govement deficit to marketplace demand
        // TODO: it would be better to check the type but that will not be set until
        //  initCalc
        if( mLeafInputs[i]->getName() != "Capital" ) {
            // really currency
            double tempDemand = mLeafInputs[ i ]->getPhysicalDemand( aPeriod );
            assert( util::isValidNumber( tempDemand ) );

            marketplace->addToDemand( mLeafInputs[ i ]->getName(), aRegionName, tempDemand, aPeriod );
        }
    }
    mLeafInputs.clear();
}

void Consumer::calcEmissions( const string& aGoodName, const string& aRegionName, const int aPeriod ) {
    // Loop over GHGs and calculate emissions.
    for( GHGIterator ghg = mGhgs.begin(); ghg != mGhgs.end(); ++ghg ){
        assert( *ghg );
        // Consumers have no physical output so pass in zero.
        (*ghg)->calcEmission( aRegionName, mLeafInputs, mOutputs, 0, 0, aPeriod );
    }
}

void Consumer::postCalc( const string& aRegionName, const string& aSectorName, const int aPeriod ) {
    const Modeltime* modeltime = scenario->getModeltime();
    if( year == modeltime->getper_to_yr( aPeriod ) ){
        // Account for exports and import.  We do this in postCalc because it is inconsequential to
        // operation however the national account is not available here.  As a temporary solution we added it
        // into the market info for Capital.
        Marketplace* marketplace = scenario->getMarketplace();
        IInfo* capitalMarketInfo = marketplace->getMarketInfo( "Capital", aRegionName, aPeriod, true );

        // nominal is the current year quantity * the current year prices, real is the current year
        // quantity * base year prices
        // consumers do not import anything so just account for domestic consumption
        double currExportsNominal = capitalMarketInfo->getDouble( "export-nominal", false );
        double currExportsReal = capitalMarketInfo->getDouble( "export-real", false );

        for( vector<IInput*>::const_iterator it = mLeafInputs.begin(); it != mLeafInputs.end(); ++it ) {
            if( !(*it)->hasTypeFlag( IInput::FACTOR ) ) {
                // subtract from exports because it was consumed domestically
                currExportsNominal -= (*it)->getPhysicalDemand( aPeriod ) *
                    marketplace->getPrice( (*it)->getName(), aRegionName, aPeriod );
                currExportsReal -= (*it)->getPhysicalDemand( aPeriod ) *
                    marketplace->getPrice( (*it)->getName(), aRegionName, 0 );
            }
        }
        capitalMarketInfo->setDouble( "export-nominal", currExportsNominal );
        capitalMarketInfo->setDouble( "export-real", currExportsReal );
    }

    // make sure we clear the base year quantities that we have kept around for
    // the calculation of our price index before they get written to the output
    // TODO: maybe this should just be in the household
    if( aPeriod == modeltime->getmaxper() - 1 && year != modeltime->getStartYear() ) {
        for( vector<IInput*>::iterator it = mLeafInputs.begin(); it != mLeafInputs.end(); ++it ) {
            (*it)->setPhysicalDemand( 0, aRegionName, 0 );
        }
    }
    mLeafInputs.clear();
}

/*!
 * \breif Calculates and adds to the real GNP national account passed in.
 * \details Sums the physical demand times base year price for that good.  That
 *          amount is then added to the GNP_REAL account and also returned so
 *          that it may be added to the account the the correct consumer.
 * \param aNationalAccount The national account of which the GNP_REAL will be
 *          added to.
 * \param aRegionName The name of the region the national account belongs to.
 * \param aPeriod The current model period.
 * \return The amount added to the GNP_REAL account.
 * \pre Demands must be calculated first.
 * \note Since the GNP is not currently used during model operation it may be
 *          a good idea to call this during post calc to cut down run time.
 */
double Consumer::calcRealGNP( NationalAccount& aNationalAccount, const string& aRegionName, int aPeriod) const {
    const int basePeriod = 0;
    double realTotal = 0;
    const Marketplace* marketplace = scenario->getMarketplace();

    for( unsigned int i = 0; i < mLeafInputs.size(); i++ ) {
        if( !mLeafInputs[ i ]->hasTypeFlag( IInput::CAPITAL ) ) {
        realTotal += mLeafInputs[ i ]->getPhysicalDemand( aPeriod ) *
            marketplace->getPrice( mLeafInputs[ i ]->getName(), aRegionName, basePeriod );
            //mLeafInputs[ i ]->getPricePaid( aRegionName, basePeriod );
        }
    }
    aNationalAccount.addToAccount( NationalAccount::GNP_REAL, realTotal );
    return realTotal;
}

void Consumer::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitConsumer( this, aPeriod );
    BaseTechnology::accept( aVisitor, aPeriod );
    aVisitor->endVisitConsumer( this, aPeriod );
}
