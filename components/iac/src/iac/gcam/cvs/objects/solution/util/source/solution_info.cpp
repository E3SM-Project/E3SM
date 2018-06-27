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
* \file solution_info.cpp
* \ingroup Solution
* \brief SolutionInfo class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include "util/base/include/util.h"
#include <string>
#include <iostream>
#include <cassert>
#include "solution/util/include/solution_info.h"
#include "solution/util/include/solution_info_set.h"
#include "solution/util/include/solver_library.h"
#include "marketplace/include/market.h"
#include "util/base/include/supply_demand_curve.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/info.h"

using namespace std;

//! Constructor
SolutionInfo::SolutionInfo( Market* aLinkedMarket )
:linkedMarket( aLinkedMarket ),
X( 0 ),
storedX( 0 ),
demand( 0 ),
storedDemand( 0 ),
supply( 0 ),
storedSupply( 0 ),
XL( 0 ),
XR( 0 ),
EDL( 0 ),
EDR( 0 ),
bracketed( false ),
mBisected( false ),
mSolutionTolerance( 0 ),
mSolutionFloor( 0 ),
mBracketInterval( 0 ),
mMaxNRPriceJump( 0 )
{
    assert( aLinkedMarket );
}

//! Equals operator based on the equality of the linked market.
bool SolutionInfo::operator==( const SolutionInfo& rhs ) const {
    return ( linkedMarket == rhs.linkedMarket );
}

//! Not equals operator based on the equals operator.
bool SolutionInfo::operator!=( const SolutionInfo& rhs ) const {
    return !( *this == rhs );
}

/*!
 * \brief Initialize the SolutionInfo.
 * \details Initializes X, supply, and demand from the linked market.  We will use the solution info
 *          values from the parser to set the market specific solution tolerance, solution floor,
 *          bracket interval, max nr price jump, and delta price if they are set.  If they do not have
 *          values for the solution tolerance, or floor the passed in defaults will be used.  Note that
 *          default values for bracket interval, max nr price jump, and delta price can not be checked here
 *          as those may be SolverComponenet specific and so defaults will be checked when a user
 *          attempts to get those values.
 * \param aDefaultSolutionTolerance The default solution tolerance to use if a market specific value
 *                                  was not provided.
 * \param aDefaultSolutionFloor The default solution floor to use if a market specific value
 *                                  was not provided.
 * \param aSolutionInfoValues The object which will contain any parsed solution parameters for this
 *                            solution info object.
 */
void SolutionInfo::init( const double aDefaultSolutionTolerance, const double aDefaultSolutionFloor,
                         const SolutionInfoParamParser::SolutionInfoValues& aSolutionInfoValues ) {
    // Update price, supply and demand from the market.
    X = linkedMarket->getRawPrice();
    supply = linkedMarket->getRawSupply();
    demand = linkedMarket->getRawDemand();
    resetBrackets();
    
    // Initialize parameters from the solution info values if they are set
    mSolutionTolerance = aSolutionInfoValues.mSolutionTolerance == 0 ?
        aDefaultSolutionTolerance : aSolutionInfoValues.mSolutionTolerance;
    mSolutionFloor = aSolutionInfoValues.mSolutionFloor == 0 ?
        aDefaultSolutionFloor : aSolutionInfoValues.mSolutionFloor;
    mBracketInterval = aSolutionInfoValues.mBracketInterval;
    mMaxNRPriceJump = aSolutionInfoValues.mMaxNRPriceJump;
    mDeltaPrice = aSolutionInfoValues.mDeltaPrice;
}

/*! \brief Return whether the SolutionInfo has had its bracketed flag set. 
* \note this does not mean it is still within the bracket.
* \return Whether the bracketed flag has been set.
*/
bool SolutionInfo::isBracketed() const {
    return bracketed;
}

/*! \brief Set the SolutionInfos bracketed flag.
*/
void SolutionInfo::setBracketed(){
    bracketed = true;
    // Fixes policy market erroneously set to true
    // without actual price brackets.
    if( XL == XR ){
        bracketed = false;
    }
}

//! Return the size of the bracket.
double SolutionInfo::getBracketSize() const {
    return fabs( XR - XL );
}

//! Return dynamically calculated bracket interval.
double SolutionInfo::getCurrentBracketInterval() const {

    return fabs( XR - XL );
}

/*!
 * \brief Gets the bracket interval which should be used when calculating
 *        brackets.
 * \details If a bracket interval has been set for this solution info that will
 *          be returned otherwise the given default will be used.
 * \param aDefaultBracketInverval The default bracket interval to use.
 * \return The correct bracket interval to use when finding brackets.
 */
double SolutionInfo::getBracketInterval( const double aDefaultBracketInverval ) const {
    return mBracketInterval == 0 ? aDefaultBracketInverval : mBracketInterval;
}

/*!
 * \brief Gets the max price change which should be used when calculating
 *        new prices using the newton raphson method.
 * \details If a max price jump has been set for this solution info that will
 *          be returned otherwise the given default will be used.
 * \param aDefaultMaxPriceJump The default max price jump to use.
 * \return The correct max price jump to use calculating new prices via NR.
 */
double SolutionInfo::getMaxNRPriceJump( const double aDefaultMaxPriceJump ) const {
    return mMaxNRPriceJump == 0 ? aDefaultMaxPriceJump : mMaxNRPriceJump;
}

/*!
 * \brief Gets the delta price change which should be used when calculating
 *        derivatives using the newton raphson method.
 * \details If a delta price has been set for this solution info that will
 *          be returned otherwise the given default will be used.
 * \param aDefaultDeltaPrice The default delta price to use.
 * \return The correct delta price to use for calculating derivatives in NR.
 */
double SolutionInfo::getDeltaPrice( const double aDefaultDeltaPrice ) const {
    return mDeltaPrice == 0 ? aDefaultDeltaPrice : mDeltaPrice;
}

/*! \brief Get the price */
double SolutionInfo::getPrice() const {
    return X;
}

//! Set the price
void SolutionInfo::setPrice( const double aPrice ){
   X = aPrice;
}

//! Set the price to the median value between the left and right bracket.
void SolutionInfo::setPriceToCenter(){
    X = ( XL + XR ) / 2;
}

/*! \brief Get Demand */
double SolutionInfo::getDemand() const {
    return demand;
}

//! Remove a given amount from the raw demand for the linked market.
void SolutionInfo::removeFromRawDemand( const double aDemand ){
    linkedMarket->removeFromRawDemand( aDemand );
}


/*! \brief Get the Supply */
double SolutionInfo::getSupply() const {
    return supply;
}

//! Remove a given amount from the raw supply for the linked market.
void SolutionInfo::removeFromRawSupply( const double aSupply ){
    linkedMarket->removeFromRawSupply( aSupply );
}

//! Get the Excess Demand, calculated dynamically.
double SolutionInfo::getED() const {
    return demand - supply;
}

//! Get the ED at the left bracket.
double SolutionInfo::getEDLeft() const {
    return EDL;
}

//! Get the ED at the right bracket.
double SolutionInfo::getEDRight() const {
    return EDR;
}

//! Store X, demand and supply so that they can later be used to calculate derivatives.
void SolutionInfo::storeValues() {
    storedX = X;
    storedDemand = demand;
    storedSupply = supply;
}

//! Restore X, demand and supply.
void SolutionInfo::restoreValues() {
    X = storedX;
    demand = storedDemand;
    supply = storedSupply;
    linkedMarket->setRawDemand( demand );
    linkedMarket->setRawSupply( supply );
    linkedMarket->setRawPrice( X );
}

/*! \brief Return the name of the SolutionInfo object.
* \author Josh Lurz
* \return The name of the market the SolutionInfo is connected to.
*/
const string& SolutionInfo::getName() const {
    return linkedMarket->getName();
}

/*! \brief Return the market type of the linked market object.
* \author Sonny Kim
* \return The type of the market the SolutionInfo is connected to.
*/
const IMarketType::Type SolutionInfo::getType() const {
    return linkedMarket->getType();
}

/*! \brief Return the name of market type of the linked market object.
* \author Sonny Kim
* \return The name of market type the SolutionInfo is connected to.
*/
string SolutionInfo::getTypeName() const {
    return linkedMarket->convert_type_to_string( linkedMarket->getType() );
}

/*! \brief Sets the price contained in the Solver into its corresponding market.
* \author Josh Lurz
* \details This function sets the market prices from the prices in the SolutionInfo
*/
void SolutionInfo::updateToMarket() {
    linkedMarket->setRawPrice( X );
}

/*! \brief Get the demands, supplies, prices and excess demands from the market and set them into their corresponding places in the SolutionInfo.
* \author Josh Lurz
* \details This function gets the price, supply and demand out the market and sets those values into the SolutionInfo
* object. This function also updates the ED value to one based on the retrieved supply and demand.
*/
void SolutionInfo::updateFromMarket() {
    X = linkedMarket->getRawPrice();
    supply = linkedMarket->getRawSupply();
    demand = linkedMarket->getRawDemand();
}

/*! \brief Determines if a price or demand market has been unbracketed and attempts to restore it back to a bracketed state. 
* \details This function checks if the linked market is a PriceMarket or DemandMarket and is unbracketed. It then adjust the brackets
* to attempt to bring the market into a bracketed state.
* \note This was originally a kludge and has never been replaced.
* \author Josh Lurz
*/
void SolutionInfo::adjustBracket() {
    // How much to adjust brackets by.
    const static double ADJUSTMENT_FACTOR = 1.5;

    // Adjust the brackets if it is a price or demand market. 
    if( linkedMarket->getType() == IMarketType::PRICE
        || linkedMarket->getType() == IMarketType::DEMAND )
    {
        double rawDemand = linkedMarket->getRawDemand();

        if( XL < rawDemand ) {
            XL = rawDemand * ADJUSTMENT_FACTOR;
        }
        if( XR > rawDemand ) {
            XR = rawDemand / ADJUSTMENT_FACTOR;
        }
    }
}

// Expand the bracket slightly.
void SolutionInfo::expandBracket( const double aAdjFactor ) {
    XL *= aAdjFactor;
    XR /= aAdjFactor;
}

/*! \brief Get the relativeED
* \return The relative excess demand. 
*/
double SolutionInfo::getRelativeED() const {
    return SolverLibrary::getRelativeED( getED(), demand, mSolutionFloor );
}

/*! \brief Determine whether a market is within the solution tolerance. 
* \author Josh Lurz
* \details This function determines if a market is solved to within the solution tolerance. It does this by checking if the 
* relative excess demand is less than the solution tolerance.
* \return Whether the market is within the solution tolerance. 
*/
bool SolutionInfo::isWithinTolerance() const {
   return ( getRelativeED() < mSolutionTolerance );
}

//! Determine whether a SolutionInfo is solvable for the current method.
bool SolutionInfo::shouldSolve( const bool isNR ) const {
    return isNR ? linkedMarket->shouldSolveNR() : linkedMarket->shouldSolve();
}

/*! \brief Calculate the log change in raw price.
* \author Josh Lurz
* \details This function determines the change in the log of the price and the stored price.
* It will return a very small number if either price is zero.
* \return The change in the logs of the price and stored price.
*/
double SolutionInfo::getLogChangeInRawPrice() const {
   double change = 0;
   // Case 1: price or Previous price is zero.
   if( storedX == 0 || X == 0 ) {
      change = util::getVerySmallNumber();
   }
   // Case 2: price and Previous price are both positive.
   else if( ( X > 0 ) && ( storedX > 0 ) ) {
      change = log( X ) - log( storedX );
   }
   // Case 3: price or Previous price is negative. This should not occur.
   else {
      assert( false );
   }
   // Check for a change of 0. Should this use a very very small number?
   if( change == 0 ){
       change = util::getVerySmallNumber();
   }
   return change;
}

/*! \brief Calculate the log change in raw demand.
* \author Josh Lurz
* \details This function determines the change in the log of the demand and the stored demand.
* It will return a very small number if either demand is zero.
* \return The change in the logs of the demand and stored demand.
*/
double SolutionInfo::getLogChangeInRawDemand() const {
   double change = 0;

   // Case 1: Demand or Previous Demand is zero.
   if( storedDemand == 0 || demand == 0 ) {
      change = util::getVerySmallNumber();
   }
   // Case 2: Demand and Previous Demand are both positive.
   else if( ( demand > 0 ) && ( storedDemand > 0 ) ) {
      change = log( demand ) - log( storedDemand );
   }
   // Case 3: Demand or Previous Demand is negative. This should not occur.
   else {
      assert( false );
   }
   return change;
}

/*! \brief Calculate the log change in raw supply.
* \author Josh Lurz
* \details This function determines the change in the log of the supply and the stored supply.
* It will return a very small number if either supply is zero.
* \return The change in the logs of the supply and stored supply.
*/
double SolutionInfo::getLogChangeInRawSupply() const {
   double change = 0;

   // Case 1: supply or Previous supply is zero.
   if( storedSupply == 0 || supply == 0 ) {
      change = util::getVerySmallNumber();
   }
   // Case 2: supply and Previous supply are both positive.
   else if( ( supply > 0 ) && ( storedSupply > 0 ) ) {
      change = log( supply ) - log( storedSupply );
   }
   // Case 3: supply or Previous supply is negative. This should not occur.
   else {
      assert( false );
   }
   return change;
}

/*! \brief Calculate demand elasticities relative to the current SolutionInfo
* \author Sonny Kim
* \details Calculate the elasticities of demand for all markets relative to the linked market. 
* If the price is 0 the cross-elasticity will be set to a very small number.
* \param solvableMarkets Vector of SolutionInfo objects to use to calculate demand elasticities.
*/
void SolutionInfo::calcDemandElas( const SolutionInfoSet& solvableMarkets ) {
    // TODO: I think a map would make it easier to only update a single cross-derivative. 
    const double dprice = getLogChangeInRawPrice();
    demandElasticities.resize( solvableMarkets.getNumSolvable() );
    
    for ( unsigned int i = 0; i < solvableMarkets.getNumSolvable(); ++i ) {
        double ddemand = solvableMarkets.getSolvable( i ).getLogChangeInRawDemand();
        demandElasticities[ i ] = ddemand / dprice;
        assert( util::isValidNumber( demandElasticities[ i ] ) );
    }
}

/*! \brief Calculate supply elasticities.
* \author Sonny Kim
* \details Calculate the elasticities of supply for all markets relative to the current SolutionInfo. 
* If the price is 0 the elasticity will be set to a very small number.
* \param solvableMarkets Vector of SolutionInfo objects to use to calculate supply elasticities.
*/
void SolutionInfo::calcSupplyElas( const SolutionInfoSet& solvableMarkets ) {
   const double dprice = getLogChangeInRawPrice();
   supplyElasticities.resize( solvableMarkets.getNumSolvable() );

   for ( unsigned int i = 0; i < solvableMarkets.getNumSolvable(); ++i ) {
      double dsupply = solvableMarkets.getSolvable( i ).getLogChangeInRawSupply();
      supplyElasticities[ i ] = dsupply / dprice;
      assert( util::isValidNumber( supplyElasticities[ i ] ) );
   }
}

/*! \brief This function checks for any SolutionInfo objects that are converging too slowly
* and resets the bracket.
* \return Whether the market was re-bracketed 
*/
bool SolutionInfo::checkAndResetBrackets(){
    // try re-bracketing if the current bracket is empty.
    if( !isCurrentlyBracketed() ){
        resetBrackets();
        return true;
    }
    return false;
}

//! Increase price by 1+multiplier, or set it to lowerBound if it is below the lower bound.
void SolutionInfo::increaseX( const double multiplier, const double lowerBound ){
    if( X >= lowerBound ) {
        X *= ( 1 + multiplier );
    }
    else if ( fabs( X ) < lowerBound ) {
        X = lowerBound;
    }
    else { // X is negative.
        assert( false );
    }
}

//! Decrease price by 1/(1+multiplier), or set it to 0 if it is below the lower bound.
void SolutionInfo::decreaseX( const double multiplier, const double lowerBound ){
    if( X >= lowerBound ) {
        X /= ( 1 + multiplier );
    }
    else if ( fabs( X ) < lowerBound ) {
        X = 0;
    }
    else {
        assert( false );
    }
}
//! Set the right bracket to X and right bracket ED to the current value of ED.
void SolutionInfo::moveRightBracketToX(){
    XR = X;
    EDR = getED();
}

//! Set the left bracket to X and left bracket ED to the current value of ED.
void SolutionInfo::moveLeftBracketToX(){
    XL = X; 
    EDL = getED();
}

//! Reset left and right bracket to X.
void SolutionInfo::resetBrackets(){
    bracketed = false;
    moveRightBracketToX();
    moveLeftBracketToX();
}

//! Dynamic check to see if the brackets span X.
bool SolutionInfo::isCurrentlyBracketed() const {
    // Market is bracketed if EDL and EDR have different signs
    // and left and right prices are different
    if ( util::sign( EDR ) != util::sign( EDL ) && XL != XR ) {
        return true;
    }
    return false;
}

//! Check if the market has correctly solved.
bool SolutionInfo::isSolved() const {
    // Check if either the market is within tolerance or for a special market
    // type dependent solution condition.
    return ( isWithinTolerance()
             || linkedMarket->meetsSpecialSolutionCriteria() );
}

//! Get the demand elasticity of this market with respect to another market.
// This should be replaced with a map and market name.
double SolutionInfo::getDemandElasWithRespectTo( const unsigned int aMarketNumber ) const {
    if( aMarketNumber > demandElasticities.size() ){
        ILogger& mainLog = ILogger::getLogger( "solver_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid market number in SolutionInfo::getDemandElasWithRespectTo of : " << aMarketNumber << endl;
        return 0;
    }
    return demandElasticities.at( aMarketNumber );
}

//! Return the list of regions contained in the linked market for this SolutionInfo.
// Probably too much linkages here, not sure how to fix.
const vector<const objects::Atom*>& SolutionInfo::getContainedRegions() const {
    return linkedMarket->getContainedRegions();
}

//! Get the supply elasticity of this market with respect to another market.
// This should be replaced with a map and market name.
double SolutionInfo::getSupplyElasWithRespectTo( const unsigned int aMarketNumber ) const {
    if( aMarketNumber > supplyElasticities.size() ){
        ILogger& solverLog = ILogger::getLogger( "solver_log" );
        solverLog.setLevel( ILogger::ERROR );
        solverLog << "Invalid market number in SolutionInfo::getSupplyElasWithRespectTo of : " << aMarketNumber << endl;
        return 0;
    }
    return supplyElasticities.at( aMarketNumber );
}

//! Return whether this market is currently unsolved and has a singularity. 
bool SolutionInfo::isUnsolvedAndSingular(){
    // Check if the market is not solved, is solved under bisection and is not solved under NR.
    return( !isSolved() && shouldSolve( false ) && !shouldSolve( true ) );
}

void SolutionInfo::setBisectedFlag(){
    mBisected = true;
}

void SolutionInfo::unsetBisectedFlag(){
    mBisected = false;
}

bool SolutionInfo::hasBisected() const {
    return mBisected;
}
//! Print out information from the SolutionInfo to an output stream.
void SolutionInfo::print( ostream& aOut ) const {

    aOut.setf(ios_base::left,ios_base::adjustfield); // left alignment
    aOut.precision(6); // for floating-point
    aOut.width(36); aOut << getName() << ", ";
    aOut.width(10); aOut << X << ", ";
    aOut.width(10); aOut << XL << ", ";
    aOut.width(10); aOut << XR << ", ";
    aOut.width(10); aOut << getED() << ", ";
    aOut.width(10); aOut << EDL << ", ";
    aOut.width(10); aOut << EDR << ", ";
    aOut.width(10); aOut << getRelativeED() << ", ";
    aOut.width(3); aOut << isBracketed() << ", ";
    aOut.width(10); aOut << supply << ", ";
    aOut.width(10); aOut << demand << ", ";
    aOut.width(10); aOut << getTypeName() << ", ";
    aOut.setf(ios_base::fmtflags(0),ios_base::floatfield); //reset to default
}

SupplyDemandCurve SolutionInfo::createSDCurve(){
    return SupplyDemandCurve( linkedMarket );
}

void SolutionInfo::printDerivatives( ostream& aOut ) const {
    aOut << getName() << " SupplyDer: ";
    for( unsigned int i = 0; i < supplyElasticities.size(); ++i ){
        aOut << "," << supplyElasticities[ i ];
    }
    aOut << endl;
        aOut << getName() << " DemandDer: ";
    for( unsigned int i = 0; i < demandElasticities.size(); ++i ){
        aOut << "," << demandElasticities[ i ];
    }
    aOut << endl;
}
