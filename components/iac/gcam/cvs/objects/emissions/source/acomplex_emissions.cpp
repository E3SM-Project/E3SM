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
 * \file acomplex_emissions.cpp
 * \ingroup Objects
 * \brief AComplexEmissions class header file.
 * \author Jim Naslund
 */

#include "util/base/include/definitions.h"

#include <cmath>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "emissions/include/acomplex_emissions.h"
#include "containers/include/scenario.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"
#include "emissions/include/ghg_mac.h"
#include "util/base/include/model_time.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/iinfo.h"
#include "technologies/include/ioutput.h"
#include "containers/include/gdp.h"
#include "emissions/include/aemissions_coef.h"
#include "emissions/include/input_emissions_coef.h"
#include "emissions/include/read_emissions_coef.h"
#include "emissions/include/aggr_emissions_coef.h"
#include "functions/include/function_utils.h"
#include "marketplace/include/cached_market.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default constructor.
AComplexEmissions::AComplexEmissions():
AGHG(),
emAdjust( 0 ),
gdpCap( 0 ),
maxCntrl( -1000 ),
techDiff( 0 ),
mGDPMidControl( 0 ),
finalEmissCoef( -1 ),
mGDPControlRange( 0 ),
adjMaxCntrl( 1 ),
multMaxCntrl( 1 ),
gwp( 1 )
{
    // default unit for emissions
    mEmissionsUnit = "Tg";
}

//! Default destructor.
AComplexEmissions::~AComplexEmissions(){
}

//! Copy constructor.
AComplexEmissions::AComplexEmissions( const AComplexEmissions& other )
: AGHG( other ){
    copy( other );
}

//! Assignment operator.
AComplexEmissions& AComplexEmissions::operator=( const AComplexEmissions& other ){
    if( this != &other ){
        // If there was a destructor it would need to be called here.
        AGHG::operator=( other );
        copy( other );
    }
    return *this;
}

//! Copy helper function.
void AComplexEmissions::copy( const AComplexEmissions& other ){
    gwp = other.gwp;
    maxCntrl = other.maxCntrl;
    mGDPMidControl = other.mGDPMidControl;
    mGDPControlRange = other.mGDPControlRange;
    gdpCap = other.gdpCap;
    techDiff = other.techDiff;
    adjMaxCntrl = other.adjMaxCntrl;
    multMaxCntrl = other.multMaxCntrl;
    emAdjust = other.emAdjust;
    finalEmissCoef = other.finalEmissCoef;
    gwp = other.gwp;
    mEmissionsUnit = other.mEmissionsUnit;

    // Deep copy the auto_ptrs
    if( other.ghgMac.get() ){
        ghgMac.reset( other.ghgMac->clone() );
    }
    if( other.mEmissionsCoef.get() ){
        mEmissionsCoef.reset( other.mEmissionsCoef->clone() );
    }
}

void AComplexEmissions::copyGHGParameters( const AGHG* prevGHG ){
    assert( prevGHG ); // Make sure valid pointer was passed

    // Copy values that always need to be the same for all periods.
    // Note that finalEmissCoef is not copied, since maxCntrl has already been set appropriately

    // Ensure that prevGHG can be cast to AComplexEmissions* otherwise return early
    // TODO: Fix this, could use a double dispatch approach to avoid the cast. See 
    // the copyParam/copyParamsInto solution in IInput.
    const AComplexEmissions* prevComplexGHG = static_cast<const AComplexEmissions*>( prevGHG );
    if( !prevComplexGHG ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Bad dynamic cast occurred in copyGHGParameters." << endl;
        return;
    }


    // Check the consistency of GHG driver from input file
    // this will call the getGHGDriverName from base class (AGHG) to retrieve the GHG driver tag
    if ( getGHGDriverName() != prevGHG->getGHGDriverName() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Warning, the driver has been changed from "<< prevGHG->getGHGDriverName() << " to "
                << getGHGDriverName() << ". This may cause a calculation error." << endl;
    }

    if ( !gwp ) { 
        gwp = prevComplexGHG->gwp; // only copy if GWP has not changed
    }
    
    if ( util::isValidNumber( prevComplexGHG->maxCntrl ) ) {
      maxCntrl = prevComplexGHG->maxCntrl;
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Invalid calcuation of maxCntrl in GHG " << getName() << endl;    
    }
    
    mGDPMidControl = prevComplexGHG->mGDPMidControl;
    mGDPControlRange = prevComplexGHG->mGDPControlRange;

    adjMaxCntrl = prevComplexGHG->adjMaxCntrl;

    // Adjust for maximum control level once GDP per capita is determined
    // This could better be put into a post-calculation processing function if we implemented that in general
    adjustMaxCntrl( prevComplexGHG->gdpCap );

    // Copy values that could change, so only copy if these are still zero (and, thus, were never read-in)
    if ( !techDiff ) { 
        techDiff = prevComplexGHG->techDiff; // only copy if techDiff has not changed
    }

    // Default value for emissCoef is zero, so only copy if a new value was not read in
    // OR
    // If an emissions value was input in a previous period and none was input this period
    // the copy from the previous period
    if( !mEmissionsCoef->getCoef() 
        || prevComplexGHG->mEmissionsCoef->getOverride() ){
        // Note that the overrideCoef method does nothing in the AEmissionsCoef class.
        // A subclass would have to override it if they want to allow it to be overridden.
        mEmissionsCoef->overrideCoef( prevComplexGHG->mEmissionsCoef->getCoef() );
    }

    // If Mac curve was input then copy it, as long as one was not read in for this period
    if ( !ghgMac.get() && prevComplexGHG->ghgMac.get() ) {
        ghgMac.reset( prevComplexGHG->ghgMac->clone() );
    }

}

double AComplexEmissions::getGHGValue( const string& aRegionName,
                                       const vector<IInput*>& aInputs,
                                       const vector<IOutput*>& aOutputs,
                                       const ICaptureComponent* aSequestrationDevice,
                                        const int aPeriod ) const 
{   
    // Constants
    const double CVRT90 = 2.212; // 1975 $ to 1990 $
    
    double GHGTax = mCachedMarket->getPrice( getName(), aRegionName, aPeriod, false );
    if( GHGTax == Marketplace::NO_MARKET_PRICE ){
        GHGTax = 0;
    }

    // apply carbon equivalent to emiss coefficient
    double generalizedCost = GHGTax * gwp * mEmissionsCoef->getCoef() / CVRT90;

    // The generalized cost returned by the GHG may be negative if
    // emissions crediting is occuring.
    return generalizedCost;
}

void AComplexEmissions::calcEmission( const string& aRegionName, 
                                      const vector<IInput*>& aInputs,
                                      const vector<IOutput*>& aOutputs,
                                      const GDP* aGDP,
                                      ICaptureComponent* aSequestrationDevice,
                                      const int aPeriod )
{
    // Primary output is always stored at position zero and used to drive
    // emissions.
    assert( aOutputs[ 0 ] );
    double primaryOutput = aOutputs[ 0 ]->getPhysicalOutput( aPeriod );

    double macReduction = 0;
    gdpCap = aGDP->getPPPGDPperCap( aPeriod );

/*
* This is a crude way to determine the input driver. This is problematic since different input units are 
* not accounted for (although as a driver, this is less relevant because all inputs currently scale to each other)
* but also does not allow an unambiguous definition of an emissions factor. 
* \TODO Need to add some way to flag the input objects (or perhaps give the input name to the GHG) to identify 
* which object is the driver to be used for emissions.
*/
    const double totalInput = FunctionUtils::getPhysicalDemandSum( aInputs, aPeriod );
    const double emissDriver = emissionsDriver(totalInput, primaryOutput);
    if ( ghgMac.get() && aPeriod > scenario->getModeltime()->getFinalCalibrationPeriod() ){
        macReduction = ghgMac->findReduction( aRegionName, aPeriod);
    }

    double fControl = 0;
    double mAdjustedGDPMidControl = adjustControlParameters( gdpCap, emissDriver, macReduction, aPeriod );
    // -1000 is used as a flag to indicate maxCntrl is not being used. This is problematic
    // because maxCntrl can and is currently being computed to take on a value smaller
    // than -1000. This could be a data error, but is happening now.
    // TODO: Separate the switch and the value from maxCntrl ( possibly add a boolean )
    if ( maxCntrl > -999 ) {
        fControl = controlFunction( maxCntrl, mGDPControlRange, mAdjustedGDPMidControl, gdpCap );
    }

    double adjEmissDriver = emissDriver * ( 1 - emAdjust ) * ( 1 - fControl );

	// apply mac reductions after updating coef which means we assume that if the user
	// used an "input-emissions" that was the emissions before any mac reductions
    mEmissionsCoef->updateCoef( adjEmissDriver );
	adjEmissDriver *= ( 1 - macReduction );
 
    // This will dynamically get the emissions value.
    mEmissions[ aPeriod ] = mEmissionsCoef->getEmissions( adjEmissDriver );
    mEmissionsByFuel[ aPeriod ] = mEmissions[ aPeriod ];

    addEmissionsToMarket( aRegionName, aPeriod );
}

bool AComplexEmissions::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ){
    if( nodeName == "GWP" ){
        gwp = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "input-emissions" ){
        mEmissionsCoef.reset( new InputEmissionsCoef( XMLHelper<double>::getValue( curr ) ) );
    }
    else if( nodeName == "emiss-adjust" ){
        emAdjust = XMLHelper<double>::getValue( curr );
    }
    else if( ( nodeName == "max-control" ) ){
        maxCntrl = XMLHelper<double>::getValue( curr );
     }
    else if( ( nodeName == "gdp-mid-control" ) ){
        mGDPMidControl = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "gdp-control-range" ){
        mGDPControlRange = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "final-emiss-coef" ){
        finalEmissCoef = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "tech-diffusion" ){
        techDiff = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "uncontrolled-emiss-coef" ){
        mEmissionsCoef.reset( new ReadEmissionsCoef( XMLHelper<double>::getValue( curr ) ) );
    }
    // Adjust max Control level, leaving current emissions constant
    else if( nodeName == "adjMaxCntrl" ){
        adjMaxCntrl = XMLHelper<double>::getValue( curr );
    }
    // Multiply maximum control level, changing current emissions
    else if( nodeName == "multMaxCntrl" ){
        multMaxCntrl = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == GhgMAC::getXMLNameStatic() ){
        parseSingleNode( curr, ghgMac, new GhgMAC );
    }
    else if( nodeName == "aggregate" ){
        mEmissionsCoef.reset( new AggrEmissionsCoef );
    }
    else{
        return false;
    }
    return true;
}

void AComplexEmissions::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    // Write out the EmissionsCoef
    mEmissionsCoef->toInputXML( out, tabs );
    XMLWriteElementCheckDefault( mEmissionsUnit, "emissions-unit", out, tabs, string("Tg") );
    XMLWriteElementCheckDefault( emAdjust, "emiss-adjust", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( maxCntrl, "max-control", out, tabs, -1000.0 );
    XMLWriteElementCheckDefault( mGDPMidControl, "gdp-mid-control", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mGDPControlRange, "gdp-control-range", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( finalEmissCoef, "final-emiss-coef", out, tabs, -1.0 );
    XMLWriteElementCheckDefault( adjMaxCntrl, "adjMaxCntrl", out, tabs, 1.0 );
    XMLWriteElementCheckDefault( multMaxCntrl, "multMaxCntrl", out, tabs, 1.0 );
    XMLWriteElementCheckDefault( techDiff, "tech-diffusion", out, tabs, 0.0 );

    // Write out the GHGMAC
    if( ghgMac.get() ){
        ghgMac->toInputXML( out, tabs );
    }
}

void AComplexEmissions::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    // Write out the EmissionsCoef
    mEmissionsCoef->toDebugXML( out, tabs );
    XMLWriteElement( emAdjust, "emiss-adjust", out, tabs );
    XMLWriteElement( maxCntrl, "max-control", out, tabs );
    XMLWriteElement( mGDPMidControl, "gdp-mid-control", out, tabs );
    XMLWriteElement( mGDPControlRange, "gdp-control-range", out, tabs );
    XMLWriteElement( finalEmissCoef, "final-emiss-coef", out, tabs );
    XMLWriteElement( adjMaxCntrl, "adjMaxCntrl", out, tabs );
    XMLWriteElement( techDiff, "tech-diffusion", out, tabs );

     // Write out the GHGMAC
    if( ghgMac.get() ){
        ghgMac->toDebugXML( period, out, tabs );
    }
}

/*! \brief finds an appropriate value for maxCntrl and adjusts mGDPMidControl as necessary
* The control function is needed in the calcEmission function and takes 4 parameters, maxCntrl, mGDPControlRange, mGDPMidControl, and gdpCap.
* mGDPControlRange and mGDPMidControl are necessary inputs to the control function, maxCntrl can either be inputed directly, or
* can be computed in this function using the input "finalEmissCoef."
* if either mGDPControlRange or mGDPMidControl are not input, or are 0, or if the emissions coefficient is 0, the function will set
* fControl to 0, and hence, there will be no emissions adjustment due to controls. In the case that both maxCntrl and
* finalEmissCoef are input (which does not make sense)  only finalEmissCoef will be used.
* The function additionally calls calcTechChange which reduces mGDPMidControl over time to account for technological change/diffusion.
* \author Nick Fernandez & Steve Smith
* \param gdpCap - The gdp per capita for the current period 
* \param emissDrive The amount of fuel that emissions are proportional to
* \param period the current period where calculations are occurring
* \return mAdjustedGDPMidControl
* \todo let initCalc know the period so that calcTechChange calculation can be moved there and will only have to be done once per period
*/
double AComplexEmissions::adjustControlParameters( const double gdpCap, const double emissDrive, const double macReduction, const int period ){
    double mAdjustedGDPMidControl = mGDPMidControl; //! mGDPMidControlRange used by the control function- adjusted for techDiffusion

    if (techDiff !=0){
        mAdjustedGDPMidControl = mGDPMidControl / calcTechChange(period);
    }
    if ( finalEmissCoef > 0 ){
        double B = 0;
        double multiplier = 0;
        // TODO: This is kind of a hack to avoid this computation.  There must be a better way.
        //       What about creating a struct to pass some of the needed variables?
        if( mEmissionsCoef->needsCalcForAdjustment() ){
            B = (1/controlFunction(100,mGDPControlRange,mAdjustedGDPMidControl,gdpCap));
            multiplier = emissDrive * ( 1 - emAdjust ) * ( 1 - macReduction );
        }
        maxCntrl = mEmissionsCoef->calcMaxCntrl( finalEmissCoef, B, multiplier );
        // Control values should never be larger than 100%.
        maxCntrl = min( maxCntrl, 100.0 );
    }
    return mAdjustedGDPMidControl;
}

/*! \brief Returns the control level for this gas in the current period
* \detailed The control function is a logistic exponential function.  It approaches 0 for values of gdpCap much less
* than mGDPMidControl, and approaches maxCntrl for values of gdpCap much greater than gdpcap0. CLOGIT is a constant equal
* to 2 times the natural log of 9, such that fControl is equal to 1/2 maxCntrl at mGDPMidControl.  
* the function returns the value of 0 in the case that either mGDPMidControl or mGDPControlRange are not input, or are equal to 0.
* \author Nick Fernandez
* \param maxCntrlIn maximum emissions reduction fraction due to controls
* \param aGDPControlRange the range over which the control percentage goes from 10% maxCntrl to 90% maxCntrl
* \param aGDPMidControl the midpoint gdpCap value of the control curve
* \param gdpCapIn the gdp per capita in PPP terms for the current period
*/
double AComplexEmissions::controlFunction( const double maxCntrlIn, const double aGDPControlRange, const double aGDPMidControl, const double gdpCapIn ){
    if( aGDPControlRange > util::getSmallNumber()  && aGDPMidControl > util::getSmallNumber() ){
        const double CLOGIT = 4.394; // See above for documentation.
        return (maxCntrlIn/100) / (1 + exp( -CLOGIT * (gdpCapIn - aGDPMidControl) / aGDPControlRange ));
    }
    return 0; 
}

/*! \brief adjusts maxCntrl (and then mGDPMidControl to recalibrate emissions) based on the read in multiplier adjMaxCntrl
* \ detailed adjMaxCntrl is a read in variable that represents a multiplier to maxCntrl.
			This is used to adjust the maximum emissions control level while leaving current emissions constant.
			the function multiplies maxCntrl by this value, and checks to make sure maxCntrl has not been
			given a multiplier that makes it greater that 100.  If this happens, it sets maxCntrl to 100
			and resets adjMaxCntrl to the value necessary to make maxCntrl 100.
			It then solves the control function for mGDPMidControl, keeping the base year emissions the same,
			so that changing maxCntrl does not mess with the base year calibrations.
			Note also that adjustMaxCntrl is run only once, in the base year when and if adjMaxCntrl != 1
* \author Nick Fernandez
* \param GDPcap the previous periods GDP per capita in PPP terms for this region
*/
void AComplexEmissions::adjustMaxCntrl( const double GDPcap ){
    if ( ( mGDPControlRange  > util::getSmallNumber() ) && 
         ( mGDPMidControl > util::getSmallNumber() ) && 
         ( adjMaxCntrl != 1 ) ) {
        // Note that maxCntrl is in percentage units
        maxCntrl *= adjMaxCntrl;
        if ( maxCntrl > 100 ) {
            adjMaxCntrl *= ( 100 / maxCntrl );
            maxCntrl = 100;
        }

        const double CLOGIT = 4.394;
        double factor1 =  1 + exp( -CLOGIT * ( GDPcap - mGDPMidControl ) / mGDPControlRange );
        
        if ( adjMaxCntrl != 1 ){
            mGDPMidControl = ( mGDPControlRange / CLOGIT ) * log( adjMaxCntrl * factor1 - 1 ) + GDPcap;
        }
        // After finished adjustments, adjMaxCntrl should be set to one so is not used anymore
        adjMaxCntrl = 1;
    }
}

/*! \brief Returns the value by which to adjust mGDPMidControl, based on technological diffusion
* The Variable TechCh represents the percent reduction in mGDPMidControl per year, due to technological change and diffusion.
* The overall reduction in mGDPMidControl is 1 + the techCh percentage raised to the power of the number of years after 
* the base year.  When applied to the control function, this will allow emissions controls to approach maxCntrl sooner.
* \ Author Nick Fernandez
* \param period the current period where calculations occur
* \returns amount to reduce the parameter mGDPMidControl by
*/
double AComplexEmissions::calcTechChange( const int period ){
    const Modeltime* modeltime = scenario->getModeltime();
    int year = modeltime->getper_to_yr( period ); 
    year -=  modeltime->getper_to_yr(1); // subtracts off base year to find number of years after base year
    return pow(1 + (techDiff / 100), year );
}

void AComplexEmissions::initCalc( const string& aRegionName,
                                  const IInfo* aLocalInfo,
                                  const int aPeriod )
{
    maxCntrl *= multMaxCntrl;
    // Make sure control percentage never goes above 100% so there are no negative emissions!
    maxCntrl = min( maxCntrl, 100.0 );

     // Perform any MAC initializations
    if( ghgMac.get() ){
        ghgMac->initCalc( getName() );
    }

    // Ensure the user set an emissions coefficient in the input
    if( !mEmissionsCoef.get() ){
        mEmissionsCoef.reset( new ReadEmissionsCoef( 0 ) );
    }

    mEmissionsCoef->initCalc( aLocalInfo, getName(), aPeriod );

    // If a finalEmissCoef or maxCntrl were read in, a mGDPControlRange and mGDPMidControl
    // must also have been read in.
    if( ( finalEmissCoef > 0 ) || ( maxCntrl > -999 ) ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        if ( mGDPControlRange == 0 ){
        mainLog << "Control function for " << getName() 
                << " requires an input of mGDPControlRange because either maxCntrl or "
                << " finalEmissCoef were read in." << endl;
        }
        if ( mGDPMidControl == 0 ){  
            mainLog << "Control function for " << getName()
                    << " requires an input of mGDPMidControl because either maxCntrl "
                    << " or finalEmissCoef were read in." << endl;
        }
    }

    mCachedMarket = scenario->getMarketplace()->locateMarket( getName(), aRegionName, aPeriod );
}

void AComplexEmissions::doInterpolations( const int aYear, const int aPreviousYear,
                                          const int aNextYear, const AGHG* aPreviousGHG,
                                          const AGHG* aNextGHG )
{
    const AComplexEmissions* prevComplexEmiss = static_cast<const AComplexEmissions*>( aPreviousGHG );
    const AComplexEmissions* nextComplexEmiss = static_cast<const AComplexEmissions*>( aNextGHG );
    
    /*!
     * \pre We are given a valid AComplexEmissions for the previous ghg.
     */
    assert( prevComplexEmiss );
    
    /*!
     * \pre We are given a valid AComplexEmissions for the next ghg.
     */
    assert( nextComplexEmiss );
    
    // tech diff is a tech change so we just copy from the next ghg to get the same
    // behavior
    techDiff = nextComplexEmiss->techDiff;
}
