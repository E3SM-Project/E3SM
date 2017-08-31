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
 * \file asimple_carbon_calc.cpp
 * \ingroup Objects
 * \brief ASimpleCarbonCalc class source file.
 * \author James Blackwood
 */

#include "util/base/include/definitions.h"
#include <cassert>

#include "ccarbon_model/include/asimple_carbon_calc.h"
#include "ccarbon_model/include/carbon_model_utils.h"
#include "util/base/include/ivisitor.h"
#include "util/base/include/util.h"
#include "land_allocator/include/land_use_history.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;
using namespace objects;

extern Scenario* scenario;

ASimpleCarbonCalc::ASimpleCarbonCalc():
mTotalEmissions( CarbonModelUtils::getStartYear(), CarbonModelUtils::getEndYear() ),
mLandUseHistory( 0 ),
mSoilTimeScale( CarbonModelUtils::getSoilTimeScale() ),
precalc_sigmoid_diff(  /* only allocate space if necessary i.e. mature age > 1 */ ),
mStoredEmissions( 0 ),
mHasCalculatedHistoricEmiss( false )
{
    int endYear = CarbonModelUtils::getEndYear();
    const Modeltime* modeltime = scenario->getModeltime();
    
    // Note we are not allocating space for period zero since that is historical
    // and can never be calculated more than once.
    mStoredEmissions[ 0 ] = 0;
    for( int period = 1; period < mStoredEmissions.size(); ++period ){
        int currYear = modeltime->getper_to_yr( period ) - modeltime->gettimestep( period ) + 1;
        mStoredEmissions[ period ] = new YearVector<double>( currYear, endYear, 0.0 );
    }
}

//! Default destructor
ASimpleCarbonCalc::~ASimpleCarbonCalc() {
    for( int period = 0; period < mStoredEmissions.size(); ++period ){
        delete mStoredEmissions[ period ];
    }
}

void ASimpleCarbonCalc::initLandUseHistory( const LandUseHistory* aHistory )
{
    mLandUseHistory = aHistory;
}

void ASimpleCarbonCalc::calc( const int aPeriod, const int aEndYear ) {
    const Modeltime* modeltime = scenario->getModeltime();
    
    // If this is a land-use history year...
    if( aPeriod == 0 ) {
        /*!
         * \warning Land-use history emissions can only be calculated once regardless
         *          of how many times the model will be run.
         */
        if( !mHasCalculatedHistoricEmiss && aEndYear == CarbonModelUtils::getEndYear() ) {
            // This code requires our land use history to be accurate.
            // AboveGroundCarbon is overwritten in these years
            // BelowGroundCarbon affects future model periods that are not overwritten
            const double aboveGroundCarbonDensity = mLandUseHistory->getHistoricAboveGroundCarbonDensity();
            const double belowGroundCarbonDensity = mLandUseHistory->getHistoricBelowGroundCarbonDensity();
            
            double prevLand = mLandUseHistory->getAllocation( CarbonModelUtils::getStartYear() - 1 );
            for( int year = CarbonModelUtils::getStartYear(); year <= mLandUseHistory->getMaxYear(); ++year ) {
                double currLand = mLandUseHistory->getAllocation( year );
                double landDifference = prevLand - currLand;
                calcAboveGroundCarbonEmission( landDifference * aboveGroundCarbonDensity, year, aEndYear, mTotalEmissions );
                calcBelowGroundCarbonEmission( landDifference * belowGroundCarbonDensity, year, aEndYear, mTotalEmissions );
                prevLand = currLand;
            }
            mHasCalculatedHistoricEmiss = true;
        }
        
        // Set the historical 1975 land into the land use vector so that it is
        // available for the next calc.
        mLandUse[ aPeriod ] = mLandUseHistory->getAllocation( modeltime->getper_to_yr( aPeriod ) );
    }
    else {
        // using model calculated allocations
        const int modelYear = modeltime->getper_to_yr( aPeriod );
        const int prevModelYear = modelYear - modeltime->gettimestep( aPeriod );
        int year = prevModelYear + 1;
        YearVector<double>& currEmissions = *mStoredEmissions[ aPeriod ];
        
        // clear the previously calculated emissions first
        for( ; year <= aEndYear; ++year ) {
            mTotalEmissions[ year ] -= currEmissions[ year ];
            currEmissions[ year ] = 0.0;
        }
        
        year = prevModelYear;
        double currLand = mLandUse[ aPeriod - 1 ];
        const double avgAnnualChangeInLand = ( mLandUse[ aPeriod ] - currLand )
            / modeltime->gettimestep( aPeriod );
        double prevCarbonAbove = currLand * getActualAboveGroundCarbonDensity( year );
        double prevCarbonBelow = currLand * getActualBelowGroundCarbonDensity( year );
        for( ++year; year <= modelYear; ++year ) {
            currLand += avgAnnualChangeInLand;
            double currCarbonAbove = currLand * getActualAboveGroundCarbonDensity( year );
            double currCarbonBelow = currLand * getActualBelowGroundCarbonDensity( year );
            calcAboveGroundCarbonEmission( prevCarbonAbove - currCarbonAbove, year, aEndYear, currEmissions );
            calcBelowGroundCarbonEmission( prevCarbonBelow - currCarbonBelow, year, aEndYear, currEmissions );
            prevCarbonAbove = currCarbonAbove;
            prevCarbonBelow = currCarbonBelow;
        }
        
        // add current emissions to the total
        for( year = prevModelYear + 1; year <= aEndYear; ++year ) {
            mTotalEmissions[ year ] += currEmissions[ year ];
        }
    }
}

double ASimpleCarbonCalc::getNetLandUseChangeEmission( const int aYear ) const {
    return mTotalEmissions[ aYear ];
}

void ASimpleCarbonCalc::setTotalLandUse( const double aLandUse, const int aPeriod ) {
    mLandUse[ aPeriod ] = aLandUse;
}

void ASimpleCarbonCalc::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitCarbonCalc( this, aPeriod );
    
    acceptDerived( aVisitor, aPeriod );
    
    aVisitor->endVisitCarbonCalc( this, aPeriod );
}

void ASimpleCarbonCalc::acceptDerived( IVisitor* aVisitor, const int aPeriod ) const {
    // do nothing currently
}

/*!
 * \brief Calculate the emission from above ground carbon for a given year.
 * \details Above ground carbon is emitted as a pulse.
 * \param aYear Year.
 * \param aEndYear The last future year to calculate to.
 * \param aEmissVector A vector to accumulate emissions into.
 */
void ASimpleCarbonCalc::calcAboveGroundCarbonEmission( const double aCarbonDiff,
                                                       const int aYear,
                                                       const int aEndYear,
                                                       YearVector<double>& aEmissVector )
{
    // If no emissions or sequestration occurred, then exit.
    if( util::isEqual( aCarbonDiff, 0.0 ) ){
        return;
    }
    
    // Finally, calculate net land use change emissions from changes in
    // above ground carbon.
    if ( getMatureAge() > 1 && aCarbonDiff < 0.0 ){  
        // If carbon content increases, then carbon was sequestered.
        // Carbon sequestration is stretched out in time, based on mMatureAge, because some
        // land cover types (notably forests) don't mature instantly.
        calcSigmoidCurve( aCarbonDiff, aYear, aEndYear, aEmissVector );
    }
    else { 
        // If carbon content decreases, then emissions have occurred.
        // We are assuming that all emissions happen in a single year.
        // Note if the mature age is just one year then sequestration (negative
        // emission) can just be added here as well.
        aEmissVector[ aYear ] += aCarbonDiff;
    }   
}

/*!
 * \brief Calculate the emission from below ground carbon for the given year.
 * \details Below ground, or soil carbon, is not emitted as a pulse but at a
 *          constant rate such that all carbon is emitted at the soil time scale.
 * \param aYear Year.
 * \param aEndYear The last future year to calculate to.
 * \param aEmissVector A vector to accumulate emissions into.
 */
void ASimpleCarbonCalc::calcBelowGroundCarbonEmission( const double aCarbonDiff,
                                                       const int aYear,
                                                       const int aEndYear,
                                                       YearVector<double>& aEmissVector )
{
    // If no emissions or sequestration occurred, then exit.
    if( util::isEqual( aCarbonDiff, 0.0 ) ){
        return;
    }
    
    // Set emissions (or, if negative, uptake) from now until the end year.
    const int endYear = std::min( static_cast<int> ( aYear + mSoilTimeScale ), aEndYear + 1);
    const double annualEmissions = ( 1.0 / mSoilTimeScale ) * aCarbonDiff;
    for( int currYear = aYear; currYear < endYear; ++currYear ){
        aEmissVector[ currYear ] += annualEmissions;
    }
}

/*!
 * \brief    Calculate the sigmoidal sequestration curve.
 * \details  Called by calcAboveGroundCarbonEmission.
 * \param    carbonDifference Annual change in carbon for aYear 
 * \param    aYear Year.
 * \param    aEndYear The last future year to calculate to.
 * \param    aEmissVector A vector to accumulate emissions into.
 */
void ASimpleCarbonCalc::calcSigmoidCurve( const double aCarbonDiff,
                                          const int aYear,
                                          const int aEndYear,
                                          YearVector<double>& aEmissVector )
{
    /*!
     * \pre This calculation will not be correct for a mature age of a single
     *      year.
     */
    assert( getMatureAge() > 1 );
    
    for( int currYear = aYear; currYear <= aEndYear; ++currYear ){
        // To avoid expensive calculations the difference in the sigmoid curve
        // has already been precomputed.
        aEmissVector[ currYear ] += precalc_sigmoid_diff[ currYear - aYear ] * aCarbonDiff;
    }
}

/*!
* \brief Returns a discount factor for the carbon subsidy for soil carbon.
* \details This method approximates a discount factor to adjust the carbon subsidy
*          to reflect the slow uptake of soil carbons and forest vegetation
*          growth. The carbon subsidy is based on a constant carbon tax; thus,
*          carbon uptake in the future should be valued less than carbon uptake
*          in the initial period
* \return soil carbon subsidy discount factor
*/
double ASimpleCarbonCalc::getBelowGroundCarbonSubsidyDiscountFactor( ){
    // If carbon uptake occurs in the first year, we do not discount it.
    if ( mSoilTimeScale == 1 ) {
        return 1.0;
    }

    // We are approximating this curve as alpha / (SoilTimeScale - beta)
    // Alpha and beta are chosen by minimizing least squared error 
    // between actual carbon subsidy discount and functional estimate
    // Note: these parameters assume a discount rate of 0.05
    const double ALPHA = 28.27;
    const double BETA = -25.39;
    return ALPHA / ( mSoilTimeScale - BETA ); 
	
	// We are approximating this curve as alpha / (SoilTimeScale - beta)
    // Alpha and beta are chosen by minimizing least squared error 
    // between actual carbon subsidy discount and functional estimate
    // Note: these parameters assume a discount rate of 0.025
    /* const double ALPHA = 59.9;
    const double BETA = -56.19;
    return ALPHA / ( mSoilTimeScale - BETA );*/
}

/*!
* \brief Returns a discount factor for the carbon subsidy for vegetation carbon.
* \details This method approximates a discount factor to adjust the carbon subsidy
*          to reflect the slow uptake of soil carbons and forest vegetation
*          growth. The carbon subsidy is based on a constant carbon tax; thus,
*          carbon uptake in the future should be valued less than carbon uptake
*          in the initial period
* \return above ground carbon subsidy discount factor
*/
double ASimpleCarbonCalc::getAboveGroundCarbonSubsidyDiscountFactor( ){
    // If carbon uptake occurs in the first year, we do not discount it.
    if ( getMatureAge() == 1 ) {
        return 1.0;
    }

    // We are approximating this curve as a quadratic with an offset of
    // 250 (If the mature age is 250, all carbon uptake occurs far enough 
    // in the future that you wouldn't base decisions on it. So, for a 
    // mature age of 250 the multiplier is zero
    // quadCoef is chosen by minimizing least squared error 
    // between actual carbon subsidy discount and functional estimate
    // Note: these parameters assume a discount rate of 0.05
    /* const double QUADCOEF = 2.7e-10;
    const int MAXMATUREAGE = 250; // Mature age where carbon subsidy is zero
    return QUADCOEF * pow( double(getMatureAge() - MAXMATUREAGE), 4); */
	
	// We are approximating this curve as a quadratic with an offset of
    // 250 (If the mature age is 250, all carbon uptake occurs far enough 
    // in the future that you wouldn't base decisions on it. So, for a 
    // mature age of 250 the multiplier is zero
    // quadCoef is chosen by minimizing least squared error 
    // between actual carbon subsidy discount and functional estimate
    // Note: these parameters assume a discount rate of 0.025
    const double QUADCOEF = 1.53e-05;
    const int MAXMATUREAGE = 250; // Mature age where carbon subsidy is zero
    return QUADCOEF * pow( double(getMatureAge() - MAXMATUREAGE), 2);
}

void ASimpleCarbonCalc::setSoilTimeScale( const int aTimeScale ) {
    mSoilTimeScale = aTimeScale;
}

double ASimpleCarbonCalc::getAboveGroundCarbonStock( const int aYear ) const {
    // TODO: decide what to do for carbon stock
	return 0;//mAboveGroundCarbonStock[ aYear ];
}

double ASimpleCarbonCalc::getBelowGroundCarbonStock( const int aYear ) const {
    // TODO: decide what to do for carbon stock
	return 0;//mBelowGroundCarbonStock[ aYear ];
}
