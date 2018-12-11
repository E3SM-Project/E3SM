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
 * \file input_tax.cpp
 * \ingroup Objects
 * \brief The InputTax class source file.
 * \author Kate Calvin
 */

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <cmath>

#include "functions/include/input_tax.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/xml_helper.h"
#include "technologies/include/icapture_component.h"
#include "functions/include/icoefficient.h"
#include "functions/include/efficiency.h"
#include "functions/include/intensity.h"
#include "containers/include/dependency_finder.h"
#include "containers/include/iinfo.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string InputTax::XML_REPORTING_NAME = "input-tax";

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Sonny Kim
* \return The constant XML_NAME as a static.
*/
const string& InputTax::getXMLNameStatic() {
    const static string XML_NAME = "input-tax";
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
const string& InputTax::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

//! Constructor
InputTax::InputTax()
: mPhysicalDemand( scenario->getModeltime()->getmaxper() ),
  mAdjustedCoefficients( scenario->getModeltime()->getmaxper(), 1.0 )
{
}

/*!
 * \brief Destructor.
 * \note An explicit constructor must be defined to avoid the compiler inlining
 *       it in the header file before the header file for the type contained in
 *       the auto_ptr is included.
 */
InputTax::~InputTax() {
}

/*!
 * \brief Copy constructor.
 * \note This class requires a copy constructor because it has dynamically
 *          allocated memory.
 * \param aOther tax input from which to copy.
 */
InputTax::InputTax( const InputTax& aOther ){
    // Do not clone the input coefficient as the calculated
    // coeffient will be filled out later.

    // Do not copy calibration values into the future
    // as they are only valid for one period.
    mName = aOther.mName;
    
    // Resize vectors to the correct size.
    mPhysicalDemand.resize( scenario->getModeltime()->getmaxper() );
    mAdjustedCoefficients.resize( scenario->getModeltime()->getmaxper() );
    
    // copy keywords
    mKeywordMap = aOther.mKeywordMap;
}

InputTax* InputTax::clone() const {
    return new InputTax( *this );
}

bool InputTax::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

void InputTax::XMLParse( const xercesc::DOMNode* node ) {
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

        if( nodeName == "keyword" ){
            DOMNamedNodeMap* keywordAttributes = curr->getAttributes();
            for( unsigned int attrNum = 0; attrNum < keywordAttributes->getLength(); ++attrNum ) {
                DOMNode* attrTemp = keywordAttributes->item( attrNum );
                mKeywordMap[ XMLHelper<string>::safeTranscode( attrTemp->getNodeName() ) ] = 
                    XMLHelper<string>::safeTranscode( attrTemp->getNodeValue() );
            }
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                    << getXMLNameStatic() << "." << endl;
        }
    }
}

void InputTax::toInputXML( ostream& aOut,
                               Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    if( !mKeywordMap.empty() ) {
        XMLWriteElementWithAttributes( "", "keyword", aOut, aTabs, mKeywordMap );
    }
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void InputTax::toDebugXML( const int aPeriod,
                               ostream& aOut,
                               Tabs* aTabs ) const
{
    XMLWriteOpeningTag ( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mAdjustedCoefficients[ aPeriod ], "current-coef", aOut, aTabs );
    XMLWriteElement( mPhysicalDemand[ aPeriod ], "physical-demand", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void InputTax::completeInit( const string& aRegionName,
                                 const string& aSectorName,
                                 const string& aSubsectorName,
                                 const string& aTechName,
                                 DependencyFinder* aDependencyFinder,
                                 const IInfo* aTechInfo )
{

    // Add the input dependency to the dependency finder.
    aDependencyFinder->addDependency( aSectorName, mName );
    mSectorName = aSectorName;
    
}

void InputTax::initCalc( const string& aRegionName,
                             const string& aSectorName,
                             const bool aIsNewInvestmentPeriod,
                             const bool aIsTrade,
                             const int aPeriod )
{
    // There must be a valid region name.
    assert( !aRegionName.empty() );
    mAdjustedCoefficients[ aPeriod ] = 1.0;
}

void InputTax::copyParam( const IInput* aInput,
                             const int aPeriod )
{
    aInput->copyParamsInto( *this, aPeriod );
}

void InputTax::copyParamsInto( InputTax& aInput,
                                  const int aPeriod ) const
{
    // do nothing 
}


double InputTax::getCO2EmissionsCoefficient( const string& aGHGName,
                                             const int aPeriod ) const
{
    return 0;
}

double InputTax::getPhysicalDemand( const int aPeriod ) const {
    assert( mPhysicalDemand[ aPeriod ].isInited() );
    return mPhysicalDemand[ aPeriod ];
}

double InputTax::getCarbonContent( const int aPeriod ) const {
    return 0;
}

void InputTax::setPhysicalDemand( double aPhysicalDemand,
                                     const string& aRegionName,
                                     const int aPeriod )
{

    Marketplace* marketplace = scenario->getMarketplace();
    IInfo* marketInfo = marketplace->getMarketInfo( mName, aRegionName, 0, true );

    // If tax is shared based, then divide by sector output.
    // Check if marketInfo exists and has the "isShareBased" boolean.
    if( marketInfo && marketInfo->hasValue( "isShareBased" ) ){
        if( marketInfo->getBoolean( "isShareBased", true ) ){
            // Each share is additive
            aPhysicalDemand/= marketplace->getDemand( mSectorName, aRegionName, aPeriod );
        }
    }
    // mPhysicalDemand can be a share if tax is share based.
    mPhysicalDemand[ aPeriod ].set( aPhysicalDemand );
    // Each technology share is additive.
    marketplace->addToDemand( mName, aRegionName, mPhysicalDemand[ aPeriod ],
                              aPeriod, true );
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
}

double InputTax::getCoefficient( const int aPeriod ) const {
    // Check that the coefficient has been initialized.
    assert( mAdjustedCoefficients[ aPeriod ].isInited() );

    return mAdjustedCoefficients[ aPeriod ];
}

void InputTax::setCoefficient( const double aCoefficient,
                                  const int aPeriod )
{
    // Do nothing.
}

double InputTax::getPrice( const string& aRegionName,
                              const int aPeriod ) const
{
    // A high tax decreases demand.
    return scenario->getMarketplace()->getPrice( mName, aRegionName, aPeriod, true );
}

void InputTax::setPrice( const string& aRegionName,
                            const double aPrice,
                            const int aPeriod )
{
    // Not hooking this up yet, it could work.
}

double InputTax::getCalibrationQuantity( const int aPeriod ) const
{
    return 0;
}

bool InputTax::hasTypeFlag( const int aTypeFlag ) const {
    return ( ( aTypeFlag & ~IInput::TAX ) == 0 );
}

double InputTax::getIncomeElasticity() const {
    return 0;
}

double InputTax::getPriceElasticity() const {
    return 0;
}

double InputTax::getTechChange( const int aPeriod ) const
{
    return 0;
}

