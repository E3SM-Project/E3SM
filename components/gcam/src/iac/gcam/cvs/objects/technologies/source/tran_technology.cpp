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
* \file tran_technology.cpp
* \ingroup Objects
* \brief transporation technology class source file.
* \author Sonny Kim, Josh Lurz, Steve Smith
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <cmath>
#include "technologies/include/tran_technology.h"
#include "emissions/include/aghg.h"
#include "containers/include/scenario.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "util/logger/include/ilogger.h"
#include "functions/include/ifunction.h"
#include "containers/include/iinfo.h"
#include "technologies/include/ical_data.h"
#include "technologies/include/ioutput.h"
#include "technologies/include/iproduction_state.h"
#include "technologies/include/marginal_profit_calculator.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

const string TranTechnology::XML_NAME = "tranTechnology";

//! Constructor.
TranTechnology::TranTechnology( const string& aName, const int aYear ): Technology( aName, aYear ) {
    mLoadFactor = 1;
	mServiceOutput = 0;
}

TranTechnology* TranTechnology::clone() const {
    return new TranTechnology( *this );
}

const std::string& TranTechnology::getXMLName() const {
    return XML_NAME;
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
const std::string& TranTechnology::getXMLNameStatic() {
    return XML_NAME;
}

bool TranTechnology::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    if( nodeName == "loadFactor" ){
        mLoadFactor = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "serviceOutput" ){
        mServiceOutput = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

void TranTechnology::toInputXMLDerived( ostream& out, Tabs* tabs ) const {  
    XMLWriteElementCheckDefault( mLoadFactor, "loadFactor", out, tabs, 1.0 );
    XMLWriteElementCheckDefault( mServiceOutput, "serviceOutput", out, tabs, 0.0 );
}

void TranTechnology::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const { 
    XMLWriteElement( mLoadFactor, "loadFactor", out, tabs );
    XMLWriteElement( mServiceOutput, "serviceOutput", out, tabs );
    XMLWriteElement( getOutput( period ) / mLoadFactor, "vehicleOutput", out, tabs );
    XMLWriteElement( getOutput( period ), "serviceOutput", out, tabs );
}   

void TranTechnology::completeInit( const string& aRegionName,
                                   const string& aSectorName,
                                   const string& aSubsectorName,
                                   DependencyFinder* aDepFinder,
                                   const IInfo* aSubsectorInfo,
                                   ILandAllocator* aLandAllocator )
{
    Technology::completeInit( aRegionName, aSectorName, aSubsectorName, aDepFinder,
                              aSubsectorInfo, aLandAllocator );
}

void TranTechnology::initCalc( const string& aRegionName,
                               const string& aSectorName, 
							   const IInfo* aSubsectorInfo,
                               const Demographic* aDemographics,
                               PreviousPeriodInfo& aPrevPeriodInfo,
							   const int aPeriod )   
{
	// initialize mOutput to read-in service output
    // TODO: This is not correct because this will add to the market. This will
    //       only happen in the first iteration however, so should be okay most
    //       of the time.
	if( aPeriod <= 1 ) {
        // Primary output is at location 0.
        mOutputs[ 0 ]->setPhysicalOutput( mServiceOutput,
                                          aRegionName,
                                          mCaptureComponent.get(),
                                          aPeriod );
    }

    Technology::initCalc( aRegionName, aSectorName, aSubsectorInfo,
                          aDemographics, aPrevPeriodInfo, aPeriod );

    // Check if illegal values have been read in
    if ( mLoadFactor == 0 ) {
        mLoadFactor = 1;
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "LoadFactor was zero in technology: " << mName << ". Reset to 1." << endl;
    }
}

void TranTechnology::postCalc( const string& aRegionName, const int aPeriod ) {
    Technology::postCalc( aRegionName, aPeriod );
}

double TranTechnology::getTotalInputCost( const string& aRegionName,
                                          const string& aSectorName,
                                          const int aPeriod ) const
{
    return getEnergyCost( aRegionName, aSectorName, aPeriod) + 
           getNonEnergyCost( aRegionName, aSectorName, aPeriod);
}

/*! \brief This function calculates the sum of the Carbon Values for all GHG's
*          in this Technology.
* \details The function first checks if a carbon tax exists for the Technology,
*          and if it does loops through all GHGs to calculate a sum carbon
*          value. The GHG function which it calls, getGHGValue() calculates the
*          carbon equivalent of all GHG's contained in this Technology. The
*          totalGHGCost attribute of the Technology is then set to this new
*          value. This function is slightly modified from the base version
*          to account for the vehicle intensity unit.
* \author Sonny Kim, Josh Lurz
* \param aRegionName The region containing this Technology.
* \param aSectorName The sector containing this Technology.
* \param aPeriod The period to calculate this value for.
* \return The total emissions and storage cost of all ghgs.
*/
double TranTechnology::getTotalGHGCost( const string& aRegionName,
                                    const string& aSectorName,
                                    const int aPeriod ) const
{
    const double CVRT90 = 2.212; // 1975 $ to 1990 $
    const double JPERBTU = 1055; // 1055 Joules per BTU
	const double GIGA = 1.0E9; // for getting price per GJ
    double totalGHGCost = 0;
    // totalGHGCost and carbontax must be in same unit as fuel price
    for( unsigned int i = 0; i < ghg.size(); i++ ) {
        totalGHGCost += ghg[ i ]->getGHGValue( aRegionName, mInputs, mOutputs, mCaptureComponent.get(), aPeriod );
    }
    // totalGHGCost is in 1975$/GJ(Btu/veh-mi) due to the vehicle intensity.
    return totalGHGCost * JPERBTU / GIGA * CVRT90;
}

double TranTechnology::calcSecondaryValue( const string& aRegionName,
                                           const int aPeriod ) const
{
    // NOTE: This entire function is copied so the units can be adjusted.
    const double CVRT90 = 2.212; // 1975 $ to 1990 $
    const double JPERBTU = 1055; // 1055 Joules per BTU
	const double GIGA = 1.0E9; // for getting price per GJ

    double totalValue = 0;
    // Add all costs from the GHGs.
    // NOTE: Negative value for GHG.
    for( unsigned int i = 0; i < ghg.size(); ++i ) {
        totalValue -= ghg[ i ]->getGHGValue( aRegionName, mInputs, mOutputs,
                                             mCaptureComponent.get(), aPeriod );
    }

    // Add all values from the outputs. The primary output is included in this
    // loop but will have a value of zero.
    for( unsigned int i = 0; i < mOutputs.size(); ++i ) {
        totalValue += mOutputs[ i ]->getValue( aRegionName,
                      mCaptureComponent.get(), aPeriod );
    }

    // TODO: Remove this function once units framework code is added.
    // totalValue is in 1975$/GJ(Btu/veh-mi) due to the vehicle intensity.
    return totalValue * JPERBTU / GIGA * CVRT90;
}

double TranTechnology::getEnergyCost( const string& aRegionName,
                                      const string& aSectorName,
                                      const int aPeriod ) const
{
    // NOTE: This entire function is copied so the units can be adjusted.
    const double CVRT90 = 2.212; // 1975 $ to 1990 $
    const double JPERBTU = 1055; // 1055 Joules per BTU
	const double GIGA = 1.0E9; // for getting price per GJ

    // Initialize energy cost.
    double cost = 0;
    // Calculate energy costs.
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        if( mInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ) {
            // TODO: Leontief assumption.
            // Assumes prices in 1975$ and vehicle intensity in BTU/veh-mi.
            // Converts cost to 1990$.
            cost += mInputs[ i ]->getPrice( aRegionName, aPeriod )
                    * mInputs[ i ]->getCoefficient( aPeriod )
                    / mAlphaZero * JPERBTU / GIGA * CVRT90;
        }
    }
    return cost;
}

double TranTechnology::getNonEnergyCost( const string& aRegionName,
                                      const string& aSectorName,
                                      const int aPeriod ) const
{
    // NOTE: This entire function is copied so the units can be adjusted.

    // Initialize non-energy cost.
    double cost = 0;
    // Calculate non-energy costs.
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        if( !mInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ) {
            // TODO: Leontief assumption.
            // Assumes prices in 1990$ per vehicle mile.
            cost += mInputs[ i ]->getPrice( aRegionName, aPeriod )
                    * mInputs[ i ]->getCoefficient( aPeriod )
                    / mAlphaZero;
        }
    }
    return cost;
}


void TranTechnology::calcCost( const string& aRegionName,
                               const string& aSectorName,
                               const int aPeriod )
{
    double techCost = getTotalInputCost( aRegionName, aSectorName, aPeriod )
                      * mPMultiplier - calcSecondaryValue( aRegionName, aPeriod );
    
    // Convert cost to cost per service instead of cost per vehicle.
    // For example,  convert $/vehicle-mi into $/pass-mi or $/ton-mi 

    mCosts[ aPeriod ] = max( techCost / mLoadFactor, util::getSmallNumber() );
}

void TranTechnology::production( const string& aRegionName, const string& aSectorName,
                                 double aVariableDemand, double aFixedOutputScaleFactor,
                                 const GDP* aGDP, const int aPeriod )
{
    // Can't have a scale factor and positive demand.
    assert( aFixedOutputScaleFactor == 1 || aVariableDemand == 0 );

    // Can't have negative variable demand.
    assert( aVariableDemand >= 0 && util::isValidNumber( aVariableDemand ) );

    // Check for positive variable demand and positive fixed output.
    assert( mFixedOutput == getFixedOutputDefault() || util::isEqual( aVariableDemand, 0.0 ) );

    // Check that a state has been created for the period.
    assert( mProductionState[ aPeriod ] );

    // Early exit optimization to avoid running through the demand function and
    // emissions calculations for non-operating technologies.
    if( !mProductionState[ aPeriod ]->isOperating() ) {
        return;
    }

    // Construct a marginal profit calculator. This allows the calculation of 
    // marginal profits to be lazy.
    MarginalProfitCalculator marginalProfitCalc( this );

    // Use the production state to determine output. This ensures the correct action
    // is taken when the technology is retired.
    double primaryOutput =
        mProductionState[ aPeriod ]->calcProduction( aRegionName,
                                                     aSectorName,
                                                     aVariableDemand,
                                                     &marginalProfitCalc,
                                                     aFixedOutputScaleFactor,
                                                     mShutdownDeciders,
                                                     aPeriod );

    // Convert from service demand (pass-km) to vehicle demand (vehicle-km)
    double vehicleOutput = primaryOutput / mLoadFactor;
        
    // for transportation technology use intensity instead of efficiency
    // convert from million Btu to EJ, (mInput in EJ)
    const double ECONV = 1.055e-9;

    // TODO: Improve this by adjusting the inputs for technical change and
    // possibly creating a new production function.
    // Using ECONV here is confusing since vehicle ouput is not in energy units.
    // This conversion is necessary to account for the intensity unit.
    const double fuelUsage = vehicleOutput * ECONV;

    // TODO: Would need to calculate the shutdown coef here if transportation
    //       stocks were vintaged.
    mProductionFunction->calcDemand( mInputs, fuelUsage, aRegionName,
                                     aSectorName, 1, aPeriod, 0, mAlphaZero );

    calcEmissionsAndOutputs( aRegionName, primaryOutput, aGDP, aPeriod );  
}

double TranTechnology::getCalibrationOutput( const bool aHasRequiredInput,
                                             const string& aRequiredInput,
                                             const int aPeriod ) const
{
    // TODO: Remove function and use the base class when units framework is
    //       available.

    /*! \pre If the caller requests only output for a specific fixed input, the
       *        specific input name must be passed. 
       */
    assert( !aHasRequiredInput || ( !aRequiredInput.empty() && aRequiredInput != "allInputs" ) );

    // Check if this is an existing vintage which cannot have a calibration value.
    if( !mProductionState[ aPeriod ]->isNewInvestment() ){
        return -1;
    }

    // If an input is required and the technology does not have it return early.
    if( aHasRequiredInput && !hasInput( aRequiredInput ) ) {
        return -1;
    }

    // Check if the technology has a calibrated output value.
    if( mCalValue.get() ) {
        return mCalValue->getCalOutput() * mLoadFactor;
    }

    // Conversion from input to output.
    const double ECONV = 1.055e-9;
    double totalCalOutput = -1;
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        // Check if either the caller does not care whether this technology uses a
        // certain input, or it is used.
        if( !aHasRequiredInput || mInputs[ i ]->getName() == aRequiredInput ) {
            // Calibrated output uses the first calibrated coefficient found.
            // All coefficients are checked for consistency, so the input used
            // is arbitrary.
            double calInput = mInputs[ i ]->getCalibrationQuantity( aPeriod );
            if( calInput >= 0 ) {
                // TODO: Remove leontief assumption.
                totalCalOutput = calInput / mInputs[ i ]->getCoefficient( aPeriod )
                                 * mAlphaZero / ECONV * mLoadFactor;
                break;
            }
        }
    }
    return totalCalOutput;
}

double TranTechnology::calcShare( const string& aRegionName,
                                  const string& aSectorName,
                                  const GDP* aGDP,
                                  const double aLogitExp,
                                  const int aPeriod ) const
{
    return Technology::calcShare( aRegionName, aSectorName, aGDP, aLogitExp, aPeriod );
}

void TranTechnology::acceptDerived( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitTranTechnology( this, aPeriod );
    aVisitor->endVisitTranTechnology( this, aPeriod );
}

void TranTechnology::doInterpolations( const Technology* aPrevTech, const Technology* aNextTech ) {
    Technology::doInterpolations( aPrevTech, aNextTech );
    
    const TranTechnology* prevTranTech = static_cast<const TranTechnology*> ( aPrevTech );
    const TranTechnology* nextTranTech = static_cast<const TranTechnology*> ( aNextTech );
    /*!
     * \pre We are given a valid TranTechnology for the previous tech.
     */
    assert( prevTranTech );
    
    /*!
     * \pre We are given a valid TranTechnology for the next tech.
     */
    assert( nextTranTech );
    
    // Interpolate load factors
    mLoadFactor = util::linearInterpolateY( year, prevTranTech->year, nextTranTech->year,
                                            prevTranTech->mLoadFactor, nextTranTech->mLoadFactor );
}
