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
 * \file non_energy_input.cpp
 * \ingroup Objects
 * \brief The NonEnergyInput class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "functions/include/non_energy_input.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string NonEnergyInput::XML_REPORTING_NAME = "input-non-energy";

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& NonEnergyInput::getXMLNameStatic() {
    const static string XML_NAME = "minicam-non-energy-input";
    return XML_NAME;
}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& NonEnergyInput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

//! Constructor
NonEnergyInput::NonEnergyInput()
: mAdjustedCosts( scenario->getModeltime()->getmaxper() ),
  mAdjustedCoefficients( scenario->getModeltime()->getmaxper() ){
}

/*! \brief Constructor that sets name attribute.
*
* This allows non-energy input objects to be created by technology
* objects.
* \author Steve Smith
*/
NonEnergyInput::NonEnergyInput( const std::string& aName ) 
: mAdjustedCosts( scenario->getModeltime()->getmaxper() ),
  mAdjustedCoefficients( scenario->getModeltime()->getmaxper() )
{
    mName = aName;
}

//! Clone the input.
NonEnergyInput* NonEnergyInput::clone() const {
    return new NonEnergyInput( *this );
}

bool NonEnergyInput::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

void NonEnergyInput::copyParam( const IInput* aInput,
                                const int aPeriod )
{
    aInput->copyParamsInto( *this, aPeriod );
}

void NonEnergyInput::copyParamsInto( NonEnergyInput& aInput,
                                     const int aPeriod ) const
{
    // Copy the coefficients forward. This is done to adjust for technical
    // change which already occurred.
    assert( aPeriod > 0 );
    aInput.mAdjustedCoefficients[ aPeriod ] = mAdjustedCoefficients[ aPeriod - 1 ];
}

void NonEnergyInput::XMLParse( const xercesc::DOMNode* node ) {
    // TODO: Replace this with the restructured XMLParse.
    // Make sure we were passed a valid node.
    assert( node );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( node, "name" );

    // get all child nodes.
    const DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() == DOMNode::TEXT_NODE ){
            continue;
        }

        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if ( nodeName == "input-cost" ) {
            mCost = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "tech-change" ){
            mTechChange = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                    << getXMLNameStatic() << "." << endl;
        }
    }
}

void NonEnergyInput::toInputXML( ostream& aOut,
                                 Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mCost, "input-cost", aOut, aTabs );
    XMLWriteElementCheckDefault( mTechChange, "tech-change", aOut, aTabs, Value( 0 ) );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void NonEnergyInput::toDebugXML( const int aPeriod,
                                 ostream& aOut,
                                 Tabs* aTabs ) const
{
    XMLWriteOpeningTag ( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mCost, "input-cost", aOut, aTabs );
    XMLWriteElement( mTechChange, "tech-change", aOut, aTabs );
    XMLWriteElement( mAdjustedCosts[ aPeriod ], "adjusted-cost", aOut, aTabs );
    XMLWriteElement( mAdjustedCoefficients[ aPeriod ], "adjusted-coef", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void NonEnergyInput::completeInit( const string& aRegionName,
                                   const string& aSectorName,
                                   const string& aSubsectorName,
                                   const string& aTechName,
                                   DependencyFinder* aDependencyFinder,
                                   const IInfo* aTechInfo )
{
    // Initialize the adjusted costs in all periods to the base read-in costs.
    // These costs may be adjusted by the Technology, for instance for capture
    // penalties.
    mAdjustedCosts.assign( mAdjustedCosts.size(), mCost );
}

void NonEnergyInput::initCalc( const string& aRegionName,
                               const string& aSectorName,
                               const bool aIsNewInvestmentPeriod,
                               const bool aIsTrade,
                               const int aPeriod )
{
    // Initialize the current coefficient to 1 if it has not 
    // been initialized through copyParam. It may be adjusted
    // later when coefficients are copied forward.

    mAdjustedCosts[ aPeriod ] = mCost;
    mAdjustedCoefficients[ aPeriod ] = 1;
}

double NonEnergyInput::getPrice( const string& aRegionName,
                                 const int aPeriod ) const
{
    assert( mAdjustedCosts[ aPeriod ].isInited() );
    return mAdjustedCosts[ aPeriod ];
}

void NonEnergyInput::setPrice( const string& aRegionName,
                               const double aPrice,
                               const int aPeriod ) 
{
    mAdjustedCosts[ aPeriod ] = aPrice;
}

double NonEnergyInput::getPhysicalDemand( const int aPeriod ) const {
    return 0;
}

void NonEnergyInput::setPhysicalDemand( double aPhysicalDemand,
                                        const string& aRegionName,
                                        const int aPeriod )
{
    // Does not add to the marketplace.
}

double NonEnergyInput::getCO2EmissionsCoefficient( const string& aGHGName,
                                                const int aPeriod ) const
{
    // Non-energy inputs cannot have emissions coefficients.
    return 0;
}

double NonEnergyInput::getCoefficient( const int aPeriod ) const {
    assert( mAdjustedCoefficients[ aPeriod ].isInited() );
    return mAdjustedCoefficients[ aPeriod ];
}

void NonEnergyInput::setCoefficient( const double aCoefficient,
                                     const int aPeriod )
{
    mAdjustedCoefficients[ aPeriod ] = aCoefficient;
}

double NonEnergyInput::getCalibrationQuantity( const int aPeriod ) const
{
    // Non-energy inputs are not calibrated.
    return -1;
}

bool NonEnergyInput::hasTypeFlag( const int aTypeFlag ) const {
    return (  ( aTypeFlag & ~( IInput::CAPITAL | IInput::MATERIAL ) ) == 0 );
}

double NonEnergyInput::getIncomeElasticity() const {
    return 0;
}

double NonEnergyInput::getPriceElasticity() const {
    return 0;
}

double NonEnergyInput::getTechChange( const int aPeriod ) const
{
    return mTechChange;
}

void NonEnergyInput::doInterpolations( const int aYear, const int aPreviousYear,
                                       const int aNextYear, const IInput* aPreviousInput,
                                       const IInput* aNextInput )
{
    const NonEnergyInput* prevNonEneInput = static_cast<const NonEnergyInput*>( aPreviousInput );
    const NonEnergyInput* nextNonEneInput = static_cast<const NonEnergyInput*>( aNextInput );
    
    /*!
     * \pre We are given a valid NonEnergyInput for the previous input.
     */
    assert( prevNonEneInput );
    
    /*!
     * \pre We are given a valid NonEnergyInput for the next input.
     */
    assert( nextNonEneInput );
    
    // interpolate the base cost
    mCost.set( util::linearInterpolateY( aYear, aPreviousYear, aNextYear,
                                         prevNonEneInput->mCost, nextNonEneInput->mCost ) );
}


