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
* \file ghg.cpp
* \ingroup Objects
* \brief Ghg class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "emissions/include/aghg.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/gdp.h"
#include "functions/include/iinput.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/iinfo.h"
#include "util/logger/include/ilogger.h"
#include "technologies/include/ioutput.h"
#include "emissions/include/aemissions_driver.h"
#include "emissions/include/emissions_driver_factory.h"
#include "technologies/include/icapture_component.h"
#include "marketplace/include/cached_market.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;
// static initialize.

typedef vector<IInput*>::const_iterator CInputIterator;

//! Default constructor.
AGHG::AGHG():
// this is inefficient as it is greater than the lifetime
// but much simpler than converting period to lifetime period 
// TODO: Fix this so it has one spot per active period.
mEmissions( scenario->getModeltime()->getmaxper() ),
mEmissionsByFuel( scenario->getModeltime()->getmaxper() ),
mEmissionsSequestered( scenario->getModeltime()->getmaxper() )
{
}

//! Destructor
AGHG::~AGHG(){
}

//! Copy constructor.
AGHG::AGHG( const AGHG& other ){
    copy( other );
}

//! Assignment operator.
AGHG& AGHG::operator=( const AGHG& other ){
    if( this != &other ){
        // If there was a destructor it would need to be called here.
        copy( other );
    }
    return *this;
}

//! Copy helper function.
void AGHG::copy( const AGHG& other ){

    // Note results are never copied.
    mEmissions.resize( scenario->getModeltime()->getmaxper() );
    mEmissionsByFuel.resize( scenario->getModeltime()->getmaxper() );
    mEmissionsSequestered.resize( scenario->getModeltime()->getmaxper() );

    // Deep copy the auto_ptr
    if( other.mEmissionsDriver.get() ){
        mEmissionsDriver.reset( other.mEmissionsDriver->clone() );
    }
}

//! \brief initialize Ghg object with xml data
void AGHG::XMLParse(const DOMNode* node) {   
    /*! \pre Assume we are passed a valid node. */
    assert( node );

    DOMNodeList* nodeList = node->getChildNodes();

    // Parse the name attribute.
    parseName( XMLHelper<string>::getAttr( node, "name" ) );

    for( unsigned int i = 0; i < nodeList->getLength(); ++i ) {
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );      

        if( nodeName == "#text" ){
            continue;
        }
        else if( EmissionsDriverFactory::isEmissionsDriverNode( nodeName ) ){
            auto_ptr<AEmissionsDriver> newDriver = EmissionsDriverFactory::create( nodeName );
            setEmissionsDriver( newDriver );
        }
        else if( nodeName == "emissions-unit" ){
            mEmissionsUnit = XMLHelper<string>::getValue( curr );
        }
        else if( XMLDerivedClassParse( nodeName, curr ) ){
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing GHG." << endl;
        }
    }
}

//! Writes datamembers to datastream in XML format.
void AGHG::toInputXML( ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLName(), out, tabs, getName() );

    // write xml for data members
    toInputXMLDerived( out, tabs );
    // done writing xml for data members.

    XMLWriteClosingTag( getXMLName(), out, tabs );

}
//! Writes datamembers to debugging datastream in XML format.
void AGHG::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLName(), out, tabs, getName() );

    // write xml for data members
    XMLWriteElement( mEmissions[ period ], "emission", out, tabs );
    XMLWriteElement( mEmissionsByFuel[ period ], "emissFuel", out, tabs );
    XMLWriteElement( mEmissionsSequestered[ period ], "emissions-sequestered", out, tabs );

    toDebugXMLDerived( period, out, tabs );
    // done writing xml for data members.

    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*!
 * \brief Calculate the aggregate output emissions coefficient for the gas.
 * \details The output coefficient is the sum of all output coefficients of all
 *          the outputs.
 * \param aOutputs Vector of Technology outputs.
 * \param aPeriod Period.
 * \return Aggregate output coefficient.
 */
double AGHG::calcOutputCoef( const vector<IOutput*>& aOutputs, const int aPeriod ) const {
    // The output coefficient is the sum of the output coefficients of all outputs.
    double outputCoef = 0;
    for( unsigned int i = 0; i < aOutputs.size(); ++i ){
        outputCoef += aOutputs[ i ]->getEmissionsPerOutput( getName(), aPeriod );
    }
    return outputCoef;
}

/*!
 * \brief Sets the emissions as the demand side of the gas market.
 * \param aRegionName the region to set
 * \param aPeriod the period
 */
void AGHG::addEmissionsToMarket( const string& aRegionName, const int aPeriod ){
    // Emissions can be positive or negative.
    if( !util::isEqual( mEmissions[ aPeriod ], 0.0 ) ){
        // set emissions as demand side of gas market
        mCachedMarket->addToDemand( getName(), aRegionName,
                                               mEmissions[ aPeriod ],
                                               aPeriod, false );
    }
}

/*! Second Method: Convert GHG tax and any storage costs into energy units using
*   GHG coefficients and return the value or cost of the tax and storage for the
*   GHG. Apply taxes only if emissions occur. Emissions occur if there is a
*   difference in the emissions coefficients.
*  \param aInput Input for which to calculate the carbon tax.
*  \param aRegionName The name of the current region.
*  \param aGoodName The name of the output product.
*  \param aSequestrationDevice A capture component which will adjust the cost.
*  \param aPeriod The period in which this calculation is occurring. 
*  \return Generalized cost or value of the GHG
*  \todo Collapsing two methods.
*/
double AGHG::getGHGValue( const IInput* aInput, const string& aRegionName,
                          const string& aGoodName,
                          const ICaptureComponent* aSequestrationDevice,
                          const int aPeriod ) const
{
    // Determine if there is a tax.
    double ghgTax = mCachedMarket->getPrice( getName(), aRegionName, aPeriod, false );
    if( ghgTax == Marketplace::NO_MARKET_PRICE ){
        ghgTax = 0;
    }

    // Retrieve proportional tax rate.
    const IInfo* marketInfo = mCachedMarket->getMarketInfo( getName(), aRegionName, aPeriod, false );
    // Note: the key includes the region name.
    const double proportionalTaxRate = 
        ( marketInfo && marketInfo->hasValue( "proportional-tax-rate" + aRegionName ) ) 
        ? marketInfo->getDouble( "proportional-tax-rate" + aRegionName, true )
        : 1.0;
    // Adjust greenhouse gas tax with the proportional tax rate.
    ghgTax *= proportionalTaxRate;

    // Get the emissions coef for the input.
    double currInputGasCoef = aInput->getCO2EmissionsCoefficient( getName(), aPeriod );

    // Get the remove fraction
    double removeFract = aSequestrationDevice ? aSequestrationDevice->getRemoveFraction( getName() ) : 0;

    // Get the storage cost of sequestered emissions
    double storageCost = aSequestrationDevice ? aSequestrationDevice->getStorageCost( aRegionName, getName(), 
        aPeriod ) : 0;

    // Return the rate.
    return ( ( 1 - removeFract ) * ghgTax + removeFract * storageCost ) * currInputGasCoef;
}

/*! Second Method: Convert GHG tax and any storage costs into energy units using
*   GHG coefficients and return the value or cost of the tax and storage for the
*   GHG. Apply taxes only if emissions occur. Emissions occur if there is a
*   difference in the emissions coefficients.
*  \param aOutput Output for which to calculate the carbon tax.
*  \param aRegionName The name of the current region.
*  \param aGoodName The name of the output product.
*  \param aSequestrationDevice A capture component which will adjust the cost.
*  \param aPeriod The period in which this calculation is occurring. 
*  \return Generalized cost or value of the GHG
*  \todo Collapsing two methods.
*/
double AGHG::getGHGValue( const IOutput* aOutput, const string& aRegionName,
                          const string& aGoodName,
                          const ICaptureComponent* aSequestrationDevice,
                          const int aPeriod ) const
{
    // Determine if there is a tax.
    double ghgTax = mCachedMarket->getPrice( getName(), aRegionName, aPeriod, false );
    if( ghgTax == Marketplace::NO_MARKET_PRICE ){
        ghgTax = 0;
    }

    // Retrieve proportional tax rate.
    const IInfo* ghgMarketInfo = mCachedMarket->getMarketInfo( getName(), aRegionName, aPeriod, false );
    // Note: the key includes the region name.
    const double proportionalTaxRate = 
        ( ghgMarketInfo && ghgMarketInfo->hasValue( "proportional-tax-rate" + aRegionName ) ) 
        ? ghgMarketInfo->getDouble( "proportional-tax-rate" + aRegionName, true )
        : 1.0;
    // Adjust greenhouse gas tax with the proportional tax rate.
    ghgTax *= proportionalTaxRate;

    // Get the emissions coef for the output.
    double currOutputGasCoef = aOutput->getEmissionsPerOutput( getName(), aPeriod );

    // Get the remove fraction
    double removeFract = aSequestrationDevice ? aSequestrationDevice->getRemoveFraction( getName() ) : 0;

    // Get the storage cost of sequestered emissions
    double storageCost = aSequestrationDevice ? aSequestrationDevice->getStorageCost( aRegionName, getName(), 
        aPeriod ) : 0;

    // Return the rate.
    return ( ( 1 - removeFract ) * ghgTax + removeFract * storageCost ) * currOutputGasCoef;
}

/*! \brief Calculate the input CO2 emissions for a good.
* \details Calculates the sum of all emissions contained in the inputs to the production of a good. This is calculated
* by looping over all the inputs and for each input, determining its carbon by multiplying its coefficient and its
* physical demand. This amount of carbon is then added to the total, which is returned by the function. This carbon
* may not all be emitted, as a portion may remain in the output good. This function may or may not work for non-CO2
* gases, depending on when it is called as their emissions coefficients are not available 
* \param aInputs Vector of inputs to determine the amount of carbon in.
* \param aRegionName Name of the region in which the emission is occurring.
* \param aPeriod Period in which the emission is occurring. 
*/
double AGHG::calcInputCO2Emissions( const vector<IInput*>& aInputs, const string& aRegionName, const int aPeriod ) const {
    double totalEmissions = 0;

    // Loop over the inputs calculating the amount of carbon in each.
    for( CInputIterator input = aInputs.begin(); input != aInputs.end(); ++input ){
        // Add on the physical amount of the input multplied by the amount of
        // emissions per unit of physical output.
        totalEmissions += (*input)->getPhysicalDemand( aPeriod ) 
                        * (*input)->getCO2EmissionsCoefficient( getName(), aPeriod );
    }
    return totalEmissions;
}

/*!
 * \brief Calculate the sum of all emissions contained in all outputs.
 * \details Determines the emissions in each output by multiplying the output's
 *          coefficient by its physical output. These emissions are then summed.
 * \param aOutputs Vector of technology outputs.
 * \param aPeriod Period.
 * \return Sum of emissions in all outputs
 */
double AGHG::calcOutputEmissions( const vector<IOutput*>& aOutputs,
                                  const int aPeriod ) const
{
    double emissions = 0;
    for( unsigned int i = 0; i < aOutputs.size(); ++i ){
        emissions += aOutputs[ i ]->getEmissionsPerOutput( getName(), aPeriod )
                   * aOutputs[ i ]->getPhysicalOutput( aPeriod );
    }
    return emissions;
}

/*!
 * \brief Calculate the aggregate input emissions coefficient for the gas.
 * \details The input coefficient is the weighted sum of all input coefficients
 *          of all the inputs.
 * \param aOutputs Vector of Technology inputs.
 * \param aPeriod Period.
 * \return Aggregate input coefficient.
 */
double AGHG::calcInputCoef( const vector<IInput*>& aInputs, const int aPeriod ) const {
    // Calculate an aggregate coefficient.
    double coefFuel = 0;
    for( unsigned int i = 0; i < aInputs.size(); ++i ){
        // Input coefficients must be greater than zero if they contribute
        // to the aggregate emissions.
        if( aInputs[ i ]->getCoefficient( aPeriod ) > 0 ){
            coefFuel += aInputs[ i ]->getCO2EmissionsCoefficient( getName(), aPeriod )
                      * aInputs[ i ]->getCoefficient( aPeriod );
        }
    }
    return coefFuel;
}

//! Return Ghg emissions.
double AGHG::getEmission( const int aPeriod ) const {
    assert( aPeriod < static_cast<int>( mEmissions.size() ) );
    return mEmissions[ aPeriod ];
}

//! Return ghg emissions implicit in fuel.
double AGHG::getEmissFuel( const int aPeriod ) const {
    return mEmissionsByFuel[ aPeriod ];
}

//! Return sequestered amount of GHG emissions.
double AGHG::getEmissionsSequestered( const int aPeriod ) const {
    assert( aPeriod < static_cast<int>( mEmissionsSequestered.size() ) );
    return mEmissionsSequestered[ aPeriod ];
}

//! returns the emissions Driver value. emissions are proportional to input minus output.
double AGHG::emissionsDriver( const double inputIn, const double outputIn ) const {
    return mEmissionsDriver->calcEmissionsDriver( inputIn, outputIn );
}

/*!
 * \brief Sets the emissions driver.
 * \note Ownership of the auto_ptr is transfered.
 * \param aEmissionsDriver New emissions driver.
 */
void AGHG::setEmissionsDriver( auto_ptr<AEmissionsDriver>& aEmissionsDriver ){
    mEmissionsDriver = aEmissionsDriver;
}

/*! \brief Update a visitor with information from a GHG for a given period.
* \param aVisitor The visitor to update.
* \param aPeriod The period for which to update.
*/
void AGHG::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitGHG( this, aPeriod );
    aVisitor->endVisitGHG( this, aPeriod );
}

/*! \return the GHG Driver tag for each period from the input file
*/
std::string AGHG::getGHGDriverName() const
{
    return ( mEmissionsDriver->getXMLName() );
}

/*!
 * \brief Hook for a ghg to do interpolations to fill in any data that
 *        should be interpolated to a newly created ghg for the missing
 *        technology.
 * \param aYear the year to be filled in.
 * \param aPreviousYear The year of the last parsed ghg.
 * \param aNextYear The year of the next closest parsed ghg.
 * \param aPreviousGHG The previous parsed ghg.
 * \param aNextGHG The next parsed ghg.
 */
void AGHG::doInterpolations( const int aYear, const int aPreviousYear,
                             const int aNextYear, const AGHG* aPreviousGHG,
                             const AGHG* aNextGHG )
{
    // the default is to not interpolate anything
}
