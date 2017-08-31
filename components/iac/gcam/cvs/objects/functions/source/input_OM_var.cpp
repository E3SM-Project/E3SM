/*
 * This software, which is provided in confidence, was prepared by employees of
 * Pacific Northwest National Labratory operated by Battelle Memorial Institute.
 * Battelle has certain unperfected rights in the software which should not be
 * copied or otherwise disseminated outside your organization without the
 * express written authorization from Battelle. All rights to the software are
 * reserved by Battelle. Battelle makes no warranty, express or implied, and
 * assumes no liability or responsibility for the use of this software.
 */

/*! 
 * \file input_OM_fixed.cpp
 * \ingroup Objects
 * \brief The InputOMVar class source file.
 * \author Sonny Kim
 */

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "functions/include/input_OM_var.h"
#include "functions/include/function_utils.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string InputOMVar::XML_REPORTING_NAME = "input-OM-var";

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, 
* \return The constant XML_NAME as a static.
*/
const string& InputOMVar::getXMLNameStatic() {
    const static string XML_NAME = "input-OM-var";
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
const string& InputOMVar::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

//! Constructor
InputOMVar::InputOMVar()
: mAdjustedCosts( scenario->getModeltime()->getmaxper() ),
  mAdjustedCoefficients( scenario->getModeltime()->getmaxper() ){
}

//! Clone the input.
InputOMVar* InputOMVar::clone() const {
    return new InputOMVar( *this );
}

bool InputOMVar::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

void InputOMVar::copyParam( const IInput* aInput,
                            const int aPeriod )
{
    aInput->copyParamsInto( *this, aPeriod );
}

void InputOMVar::copyParamsInto( InputOMVar& aInput,
                                 const int aPeriod ) const
{
    // Copy the coefficients forward. This is done to adjust for technical
    // change which already occurred.
    assert( aPeriod > 0 );
    aInput.mAdjustedCoefficients[ aPeriod ] = mAdjustedCoefficients[ aPeriod - 1 ];
}

void InputOMVar::XMLParse( const xercesc::DOMNode* node ) {
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
        if ( nodeName == "OM-var" ) {
            mOMVar = XMLHelper<double>::getValue( curr );
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

void InputOMVar::toInputXML( ostream& aOut,
                                 Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mOMVar, "OM-var", aOut, aTabs );
    XMLWriteElementCheckDefault( mTechChange, "tech-change", aOut, aTabs, Value( 0 ) );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void InputOMVar::toDebugXML( const int aPeriod,
                             ostream& aOut,
                             Tabs* aTabs ) const
{
    XMLWriteOpeningTag ( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( calcOMVarCost(), "levelized-OM-var", aOut, aTabs );
    XMLWriteElement( mOMVar, "OM-var", aOut, aTabs );
    XMLWriteElement( mTechChange, "tech-change", aOut, aTabs );
    XMLWriteElement( mAdjustedCosts[ aPeriod ], "adjusted-cost", aOut, aTabs );
    XMLWriteElement( mAdjustedCoefficients[ aPeriod ], "adjusted-coef", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void InputOMVar::completeInit( const string& aRegionName,
                               const string& aSectorName,
                               const string& aSubsectorName,
                               const string& aTechName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo )
{   
    // Initialize the adjusted costs in all periods to the base calculate
    // levelized OM-var cost.
    // These costs may be adjusted by the Technology, for instance for capture
    // penalties.
    mAdjustedCosts.assign( mAdjustedCosts.size(), calcOMVarCost() );
}

/** Calculate the levelizd OM_fixed cost.
 *
 * \param void 
 * \return Levelized OM_fixed costs.
 * \author Sonny Kim
 */
double InputOMVar::calcOMVarCost( void ) const
{
	//read in as $/MWh
    double OMVarCost = mOMVar /(1000 * FunctionUtils::GJ_PER_KWH());
	
    return OMVarCost; // 1975$/GJ
}

void InputOMVar::initCalc( const string& aRegionName,
                           const string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod )
{
    // Initialize the current coefficient to 1 if it has not 
    // been initialized through copyParam. It may be adjusted
    // later when coefficients are copied forward.
    mAdjustedCoefficients[ aPeriod ] = 1;
}

double InputOMVar::getPrice( const string& aRegionName,
                             const int aPeriod ) const
{
    assert( mAdjustedCosts[ aPeriod ].isInited() );
    return mAdjustedCosts[ aPeriod ];
}

void InputOMVar::setPrice( const string& aRegionName,
                           const double aPrice,
                           const int aPeriod ) 
{
    mAdjustedCosts[ aPeriod ] = aPrice;
}

double InputOMVar::getPhysicalDemand( const int aPeriod ) const {
    return 0;
}

void InputOMVar::setPhysicalDemand( double aPhysicalDemand,
                                    const string& aRegionName,
                                    const int aPeriod )
{
    // Does not add to the marketplace.
}

double InputOMVar::getCO2EmissionsCoefficient( const string& aGHGName,
                                            const int aPeriod ) const
{
    // Capital cost inputs cannot have emissions coefficients.
    return 0;
}

double InputOMVar::getCoefficient( const int aPeriod ) const {
    assert( mAdjustedCoefficients[ aPeriod ].isInited() );
    return mAdjustedCoefficients[ aPeriod ];
}

void InputOMVar::setCoefficient( const double aCoefficient,
                                 const int aPeriod )
{
    mAdjustedCoefficients[ aPeriod ] = aCoefficient;
}

void InputOMVar::tabulateFixedQuantity( const string& aRegionName,
                                        const double aFixedOutput,
                                        const bool aIsInvestmentPeriod,
                                        const int aPeriod )
{
}

void InputOMVar::scaleCalibrationQuantity( const double aScaleFactor ){
    // Capital cost inputs are not calibrated.
}

double InputOMVar::getCalibrationQuantity( const int aPeriod ) const
{
    // Capital cost inputs are not calibrated.
    return -1;
}

bool InputOMVar::hasTypeFlag( const int aTypeFlag ) const {
    return ( ( aTypeFlag & ~( IInput::OM_VAR ) ) == 0 );
}

double InputOMVar::getIncomeElasticity() const {
    return 0;
}

double InputOMVar::getPriceElasticity() const {
    return 0;
}

double InputOMVar::getTechChange( const int aPeriod ) const
{
    return mTechChange;
}

void InputOMVar::doInterpolations( const int aYear, const int aPreviousYear,
                                    const int aNextYear, const IInput* aPreviousInput,
                                    const IInput* aNextInput )
{
    const InputOMVar* prevOMInput = static_cast<const InputOMVar*>( aPreviousInput );
    const InputOMVar* nextOMInput = static_cast<const InputOMVar*>( aNextInput );
    
    /*!
     * \pre We are given a valid InputOMVar for the previous input.
     */
    assert( prevOMInput );
    
    /*!
     * \pre We are given a valid InputOMVar for the next input.
     */
    assert( nextOMInput );
    
    // tech change is just copied from the next input
    mTechChange = nextOMInput->mTechChange;
    
    // interpolate the costs
    mOMVar = util::linearInterpolateY( aYear, aPreviousYear, aNextYear,
                                       prevOMInput->mOMVar, nextOMInput->mOMVar );
}
