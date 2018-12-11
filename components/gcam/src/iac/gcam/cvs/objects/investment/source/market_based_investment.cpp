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
 * \file market_based_investment.cpp
 * \ingroup Objects
 * \brief MarketBasedInvestor class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/xml_helper.h"
#include "investment/include/market_based_investment.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/scenario.h"
#include "marketplace/include/imarket_type.h"
#include "investment/include/investment_utils.h"
#include "investment/include/levelized_cost_calculator.h"
#include "investment/include/output_share_levelized_cost_calculator.h"
#include "investment/include/rate_logit_distributor.h"
#include "util/logger/include/ilogger.h"
#include "functions/include/function_utils.h"
#include "util/base/include/configuration.h"
#include "investment/include/investable_counter_visitor.h"
#include "investment/include/set_share_weight_visitor.h"
#include "investment/include/get_distributed_investment_visitor.h"
#include "investment/include/iinvestable.h"

using namespace std;
extern Scenario* scenario;

//! Constructor
MarketBasedInvestor::MarketBasedInvestor():
mInvestmentLogitExp( 1 ),
mInvestments( scenario->getModeltime()->getmaxper() ),
mFixedInvestments( scenario->getModeltime()->getmaxper(), -1.0 ){
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way the tag is always consistent for both read-in and output and
*          can be easily changed. The "==" operator that is used when parsing,
*          required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz
* \return The constant XML_NAME as a static.
*/
const std::string& MarketBasedInvestor::getXMLNameStatic() {
    const static string XML_NAME = "market-based-investor";
    return XML_NAME;
}

/*! \brief Parses any data from XML.
*
* \author Josh Lurz
* \param aCurr pointer to the current node in the XML input tree
*/
void MarketBasedInvestor::XMLParse( const xercesc::DOMNode* aCurr ) {
    /*! \pre make sure we were passed a valid node. */
    assert( aCurr );

    // get all child nodes.
    const xercesc::DOMNodeList* nodeList = aCurr->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        const xercesc::DOMNode* curr = nodeList->item( i );
        // Skip any text nodes.
        if( curr->getNodeType() == xercesc::DOMNode::TEXT_NODE ){
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "FixedInvestment" ){
            XMLHelper<double>::insertValueIntoVector( curr, mFixedInvestments, scenario->getModeltime() );
        }
        else if( nodeName == "InvestmentLogitExp" ){
            mInvestmentLogitExp = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING ); 
            mainLog << "Unrecognized node " << nodeName << " found while parsing " << getXMLNameStatic() 
                    << "." << endl;
        }
    }
}

/*! \brief Write the object to an XML output stream for debugging.
* \param aPeriod Period to write debugging information for.
* \param aOut The output stream to write to.
* \param aTabs The object which tracks the number of tabs to write.
*/
void MarketBasedInvestor::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    
    XMLWriteElement( mInvestments[ aPeriod ], "investment", aOut, aTabs,
        scenario->getModeltime()->getper_to_yr( aPeriod ) );
    
    XMLWriteElement( mFixedInvestments[ aPeriod ], "FixedInvestment", aOut, aTabs,
        scenario->getModeltime()->getper_to_yr( aPeriod ) );
    
    XMLWriteElement( mInvestmentLogitExp, "InvestmentLogitExp", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Write the object to an XML output stream.
* \param aOut The output stream to write to.
* \param aTabs The object which tracks the number of tabs to write.
*/
void MarketBasedInvestor::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteVector( mFixedInvestments, "FixedInvestment", aOut, aTabs, scenario->getModeltime(), -1.0 );

    XMLWriteElementCheckDefault( mInvestmentLogitExp, "InvestmentLogitExp", aOut, aTabs, 1.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Complete the initialization of the market based investor.
* \details This function performs a series of tasks needed before it can be
*          used. It stores the region and sector name, and then sets up the
*          necessary market. This is done by constructing a unique market name,
*          creating the market, and setting the market to solve for all periods
*          where investment is not fixed at the sector level.
* \param aRegionName Name of the region containing this investor.
* \param aSectorName Name of the sector the investor is investing in.
* \author Josh Lurz
*/
void MarketBasedInvestor::completeInit( const string& aRegionName, const string& aSectorName ){
    // Store the region and sector name
    mRegionName = aRegionName;
    mSectorName = aSectorName;
    // Investment market must be for current region only.
    string marketRegionName( mRegionName );
    // Set the name of the "good" for this market. This must be unique.
    mMarketName = mSectorName + "-" + getXMLNameStatic();


    const int START_PERIOD = scenario->getModeltime()->getBasePeriod();

    // Create the trial market. 
    Marketplace* marketplace = scenario->getMarketplace();
    if ( marketplace->createMarket( mRegionName, marketRegionName, mMarketName, IMarketType::NORMAL ) ) {
        // Set the market to solve if the investment in the period is not fixed.
        if( !Configuration::getInstance()->getBool( "CalibrationActive" ) ){
            for( int per = START_PERIOD; per < scenario->getModeltime()->getmaxper(); ++per ){
                if( mFixedInvestments[ per ] == -1 ){
                    marketplace->setMarketToSolve( mMarketName, mRegionName, per );
                }
            }
        }
    }
    else { // This should not occur unless there is invalid data read in.
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Multiple sectors cannot use the same trial market for investment." << endl;
    }

    // Warn if fixed investment for the period 0 was read in, as it will be ignored.
    if( mFixedInvestments[ 0 ] != -1 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING ); 
        mainLog << "Ignoring fixed investment in the base period for sector "
                << mSectorName << " in region " << mRegionName << endl;
    }
}


/*! \brief Initialization of the market based investor for each period.
* \param aRegionName Name of the region containing this investor.
* \param aSectorName Name of the sector the investor is investing in.
* \author Josh Lurz
*/
void MarketBasedInvestor::initCalc( vector<IInvestable*>& aInvestables,
                                    NationalAccount& aNationalAccount, 
                                    const Demographic* aDemographic,
                                    const int aPeriod )
{

    Marketplace* marketplace = scenario->getMarketplace();
    /*! \pre Check that the period is not nonsensical */
    assert( aPeriod >= 0 );
    // In period 0 calculate the initial trial investment and set it in the marketplace
    if( aPeriod == 0 ){
        // Subtract fixed investments so that trial amount is the variable portion only.
        // If all fixed, then null variable or trial amount.
        double totalInvestments = InvestmentUtils::sumInvestment( aInvestables, aPeriod );
        double fixedInvestments = InvestmentUtils::sumFixedInvestment( aInvestables, aPeriod );
        mInvestments[ aPeriod ] = totalInvestments - fixedInvestments;
        assert( mInvestments[ aPeriod ] >= 0 );
        if( mInvestments[ aPeriod ] > 0 ) {
            // if we have variable investments set the initial trail value
            marketplace->setPrice( mMarketName, mRegionName, mInvestments[ aPeriod ], aPeriod, true );
        }
        else {
            // if all of the investment is fixed we need to make sure that this market does not
            // try to solve as well as set the fixed investment in this investor so that it knows
            // not to try to set the efficiency conditions
            marketplace->unsetMarketToSolve( mMarketName, mRegionName, aPeriod );
            mFixedInvestments[ aPeriod ] = fixedInvestments;
        }

        if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
            InvestableCounterVisitor investableCounterVisitor( mRegionName, aInvestables.size() );
            visitInvestables( aInvestables, &investableCounterVisitor, aPeriod );
        }

    }
}

/*! \brief Distribute the total new investment as solved in the marketplace.
* \todo Reword this comment.
* \details In period 0, investment is summed not calculated. This function
*          determines the level of variable trial investment from the solved
*          parameter in the market, or the price unless the sector investment is
*          fixed. If the sector is fixed, the market should not have been set to
*          solve, and the trial value is ignored. This trial investment does not
*          include subsector and technology fixed investment, it is only the
*          variable investment. The trial investment plus the sum of all fixed
*          investment at the subsector level and below is then distributed to
*          the children of the investor using the LevelizedCostCalculator and
*          RateLogitDistributor, stored internally, and returned.


* \param aInvestables The vector of children which will receive investment.
* \param aNationalAccount The national accounts container.
* \param aDemographic A pointer to the Demographics object.
* \param aPeriod The period in which investment is calculated and distributed.
* \return The total investment which occurred.
* \author Josh Lurz
*/
double MarketBasedInvestor::calcAndDistributeInvestment( vector<IInvestable*>& aInvestables,
                                                         NationalAccount& aNationalAccount, 
                                                         const Demographic* aDemographic,
                                                         const int aPeriod )
{   
    // For the base period use read-in investments. Do not calculate new investments.
    // Revise to calculate investment is base period is to be solved.

    double totalInvestment = 0;    
    // Create a standard profit rate calculator. This will be used to calculate the average
    // levelized cost and to distribute investment. 
    LevelizedCostCalculator levelizedCostCalculator;
    // Check if parent level investment is fixed.
    if( mFixedInvestments[ aPeriod ] != -1 ){
        // Use sector level fixed investment instead of the trial value. If all
        // children are fixed and that sum does not equal this value, there will
        // be a warning later when the investment is distributed.
        totalInvestment = mFixedInvestments[ aPeriod ];
        // Check if total child investment is greater than this amount, as that
        // will always override this amount.
        double childSumFixed = InvestmentUtils::sumFixedInvestment( aInvestables, aPeriod );
        if( childSumFixed > mFixedInvestments[ aPeriod ] ){
            totalInvestment = childSumFixed;
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Overriding parent level investment with child level investment sum. " << endl;
        }
    }
    // Otherwise use the solved value.
    else {
        // Determine the amount of fixed investment in the sector.
        const double fixedInvestment = InvestmentUtils::sumFixedInvestment( aInvestables, aPeriod );
        /*! \invariant Fixed investment is positive. */
        assert( fixedInvestment >= 0 );
        
        Marketplace* marketplace = scenario->getMarketplace();
        double trialInvestment = marketplace->getPrice( mMarketName, mRegionName, aPeriod, true );
        // Warn if trial investment reaches zero.
  /*      if( trialInvestment < util::getSmallNumber() ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Zero trial investment for with average levelized cost of " 
                    << sectorLevelizedCost << " in sector " << mSectorName 
                    << " in " << mRegionName << "." << endl;
        }
 */
        /*! \invariant The trial level of investment should always be zero or
        *              greater. 
        */
        assert( trialInvestment >= 0 );
        // Probably should warn on zero, not entirely sure if it is possible(all
        // fixed subsectors?)

        // Calculate the total investment for the sector as all fixed
        // investment, which cannot be avoided, plus the trial quantity from the
        // market.
        totalInvestment = fixedInvestment + trialInvestment;
    }
    if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
        SetShareWeightVisitor setShareWeightVisitor( mRegionName );
        visitInvestables( aInvestables, &setShareWeightVisitor, aPeriod );
    }

    // Create a logit based investment distributor.
    RateLogitDistributor invDistributor;
    
    // Use the investment distributor to distribute the trial investment plus
    // the fixed investment.
    mInvestments[ aPeriod ] = invDistributor.distribute( &levelizedCostCalculator,
                                                         aInvestables,
                                                         aNationalAccount,
                                                         mInvestmentLogitExp,
                                                         mRegionName,
                                                         mSectorName,
                                                         totalInvestment,
                                                         aPeriod );

    // Check that total investment and distributed investment are equal.
 /*   if( !util::isEqual( totalInvestment, mInvestments[ aPeriod ] ) ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << totalInvestment - mInvestments[ aPeriod ]
                << " difference between desired investment and distributed investment at the sector level in "
                << mSectorName << " in " << mRegionName << endl;
    }*/
    /*! \post A positive amount of investment occurred. */
    assert( mInvestments[ aPeriod ] >= 0 );

    if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
        GetDistributedInvestmentVisitor getDistributedInvestmentVisitor( mRegionName );
        visitInvestables( aInvestables, &getDistributedInvestmentVisitor, aPeriod );
    }

    // Return the total amount of investment actually distributed
    return mInvestments[ aPeriod ];
}

/*! \brief Set the efficienty conditions for the investment market for the sector.
* \details This method sets one side of the equation, supply, to the price received
*          for the good. It sets the other side of the equation, demand, to the
*          average levelized cost for the sector. This is calculated using a
*          OutputShareLevelizedCostCalculator object so that the price also reflects
*          fixed investment technologies. The equation will solve when price received
*          is equal to the levelized cost, or profit equals 0. 
* \param aInvestables The vector of children which will receive investment.
* \param aNationalAccount The national accounts container.
* \param aDemographic A pointer to the Demographics object.
* \param aPeriod The period in which investment is calculated and distributed.
* \author Sonny Kim
*/
void MarketBasedInvestor::setEfficiencyConditions( vector<IInvestable*>& aInvestables,
                                                   NationalAccount& aNationalAccount, 
                                                   const Demographic* aDemographic,
                                                   const int aPeriod ) const
{
    // do not calculate efficiency conditions for fixed investment
    // otherwise they will never solve
    if( mFixedInvestments[ aPeriod ] != -1 ) {
        return;
    }
    /*
    if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
        SetShareWeightVisitor setShareWeightVisitor( mRegionName );
        visitInvestables( aInvestables, &setShareWeightVisitor, aPeriod );
    }
    */

    Marketplace* marketplace = scenario->getMarketplace();
    const double currInvestment = marketplace->getPrice( mMarketName, mRegionName, aPeriod );
    // Get the price received for the good from the marketInfo for the
    // sector market. We need to check to make sure this is up to date.
    const double priceReceived = FunctionUtils::getPriceReceived( mRegionName, mSectorName, aPeriod );
    /*! \invariant Price received should always be positive and non-zero.*/
    assert( priceReceived > 0 );
    /*! \pre The market supply is zero as this is the only place it is added
    */
    assert( marketplace->getDemand( mMarketName, mRegionName, aPeriod ) 
        < util::getVerySmallNumber() );
    marketplace->addToDemand( mMarketName, mRegionName, currInvestment * priceReceived, aPeriod, true );

    // Create a levelized cost rate calculator to calculate the sector average
    // levelized cost.
    // TODO: make sure investment has already been distributed
    OutputShareLevelizedCostCalculator levelizedCostCalculator;
    const double sectorLevelizedCost = levelizedCostCalculator.calcSectorExpectedProfitRate( aInvestables,
        aNationalAccount, mRegionName, mSectorName, mInvestmentLogitExp, false, false, aPeriod );

    // Set the sector levelized cost as the right hand side. It will have
    // been cleared at the end of the last iteration, so there should not be
    // a residual.
    /*! \pre The market demand is zero as this is the only place it is added
    *        to. 
    */
    // Set the left hand side of the equation to the price received for the
    // good.
    assert( marketplace->getSupply( mMarketName, mRegionName, aPeriod ) 
        < util::getVerySmallNumber() );

    // The base year investment numbers may not reflect optimum investment patterns
    // so we allow the levelized cost and price recieved to be off in the base year.
    // Do write out a warning so that the user knows in case they were intended to be
    // equal.
    if( aPeriod == 0 && !util::isEqual(sectorLevelizedCost,  priceReceived, .001 ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Reseting LC for " << mMarketName << " in " << mRegionName << " from " << sectorLevelizedCost
            << " to " << priceReceived << endl;
        marketplace->addToSupply( mMarketName, mRegionName, currInvestment * priceReceived, aPeriod, true );
        return;
    }
    // to here
    marketplace->addToSupply( mMarketName, mRegionName, currInvestment * sectorLevelizedCost, aPeriod, true );
}
/*! \brief Set the efficienty conditions for the investment market for the sector.
* \details In period 0, investment is summed not calculated. This function
*          determines the level of variable trial investment from the solved
* \param aInvestables The vector of children which will receive investment.
* \param aNationalAccount The national accounts container.
* \param aDemographic A pointer to the Demographics object.
* \param aPeriod The period in which investment is calculated and distributed.
* \author Sonny Kim
*/
void MarketBasedInvestor::visitInvestables( vector<IInvestable*>& aInvestables,
                                            IVisitor* aVisitor,
                                            const int aPeriod ) const
{
    for( unsigned int i = 0; i < aInvestables.size(); ++i ){
        aInvestables[ i ]->accept( aVisitor, aPeriod );
    }
}
