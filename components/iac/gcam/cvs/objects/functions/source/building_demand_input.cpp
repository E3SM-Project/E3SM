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
 * \file energy_input.cpp
 * \ingroup Objects
 * \brief The BuildingDemandInput class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "functions/include/building_demand_input.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/dependency_finder.h"
#include "containers/include/iinfo.h"
#include "functions/include/function_utils.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string BuildingDemandInput::XML_REPORTING_NAME = "input-building";

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
const string& BuildingDemandInput::getXMLNameStatic() {
    const static string XML_NAME = "building-demand-input";
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
const string& BuildingDemandInput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

//! Constructor
BuildingDemandInput::BuildingDemandInput():
mPhysicalDemand( scenario->getModeltime()->getmaxper() ),
mType( eEnd ),
mSaturation( 1 )
{
}

//! Destructor required so that destructor is not inlined.
BuildingDemandInput::~BuildingDemandInput() {
}

BuildingDemandInput* BuildingDemandInput::clone() const {
    return new BuildingDemandInput( *this );
}

bool BuildingDemandInput::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

void BuildingDemandInput::XMLParse( const DOMNode* aNode ) 
{
    // TODO: Replace this with the restructured XMLParse.
    // Make sure we were passed a valid node.
    assert( aNode );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( aNode, "name" );

    // get all child nodes.
    const DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() == DOMNode::TEXT_NODE ){
            continue;
        }

        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "saturation" ){
            mSaturation = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "type" ){
            BuildingInputType newType =
                convertStringToType( XMLHelper<string>::getValue( curr ) );
            if( newType == eEnd ){
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::ERROR );
                mainLog << "Invalid type was read in for building demand input "
                    << mName << "." << endl;
            }
            else {
                mType = newType;
            }
        }
        else if( nodeName == "coef-adjustment" ){
            mCoefficientAdjustment = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "price-elasticity" ){
            mPriceElasticity = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                << getXMLNameStatic() << "." << endl;
        }
    }
}

void BuildingDemandInput::toInputXML( ostream& aOut,
                                      Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElementCheckDefault( mSaturation.get(), "saturation", aOut,
                                 aTabs, 1.0 );
    XMLWriteElement( getTypeInfo( mType ).mTypeString, "type", aOut, aTabs );
    XMLWriteElementCheckDefault( mCoefficientAdjustment, "coef-adjustment",
                                 aOut, aTabs, Value( 1.0 ) );
    XMLWriteElementCheckDefault( mPriceElasticity, "price-elasticity",
                                 aOut, aTabs, Value( 0 ) );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void BuildingDemandInput::toDebugXML( const int aPeriod,
                                      ostream& aOut,
                                      Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( getTypeInfo( mType ).mTypeString, "type", aOut, aTabs );
    XMLWriteElement( mSaturation, "saturation", aOut, aTabs );
    XMLWriteElement( mPhysicalDemand[ aPeriod ], "physical-demand", aOut, aTabs );
    XMLWriteElement( mCoefficient, "coefficient", aOut, aTabs );
    XMLWriteElement( mCoefficientAdjustment, "coef-adjustment", aOut, aTabs );
    XMLWriteElement( mPriceElasticity, "price-elasticity", aOut, aTabs );
    XMLWriteElement( mPhysicalDemand[ aPeriod ], "physical-demand", aOut, aTabs );
    XMLWriteElement( mAverageInsulation, "average-insulation", aOut, aTabs );
    XMLWriteElement( mFloorToSurfaceArea, "floor-to-surface", aOut, aTabs );
    XMLWriteElement( mDegreeDays, "degree-days", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void BuildingDemandInput::completeInit( const string& aRegionName,
                                        const string& aSectorName,
                                        const string& aSubsectorName,
                                        const string& aTechName,
                                        DependencyFinder* aDependencyFinder,
                                        const IInfo* aTechInfo )
{
    // Check for an invalid type.
    if( mType == eEnd ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "A valid type was not read in for building demand input "
                << mName << "." << endl;
        // Reset to the generic type to avoid crashing.
        mType = eGeneric;
    }
    initializeParameters( aSectorName, aSubsectorName, aTechName, aTechInfo );

    aDependencyFinder->addDependency( aSectorName, mName );
    
    // Adjust the coefficient for various other input specific parameters. This 
    // coefficient will later be scaled for calibration values.
    mCoefficient = calcDemandAdjustmentFactor();
}

void BuildingDemandInput::initCalc( const string& aRegionName,
                                    const string& aSectorName,
                                    const bool aIsNewInvestmentPeriod,
                                    const bool aIsTrade,
                                    const int aPeriod )
{
    if( FunctionUtils::getCO2Coef( aRegionName, mName, aPeriod ) > util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Building service supply sector " << mName
                << " has a non-zero CO2 coefficient." << endl;
    }
    if( aIsNewInvestmentPeriod && !mCoefficientAdjustment.isInited() ){
        mCoefficientAdjustment = 1;
    }
}

void BuildingDemandInput::copyParam( const IInput* aInput,
                                     const int aPeriod )
{
    aInput->copyParamsInto( *this, aPeriod );
}

double BuildingDemandInput::getCO2EmissionsCoefficient( const string& aGHGName,
                                                     const int aPeriod ) const
{
    // Buildings must only consume service goods which do not have CO2
    // coefficients.
    return 0;
}

double BuildingDemandInput::getPhysicalDemand( const int aPeriod ) const {
    assert( mPhysicalDemand[ aPeriod ].isInited() );
    return mPhysicalDemand[ aPeriod ];
}

void BuildingDemandInput::setPhysicalDemand( double aPhysicalDemand,
                                             const string& aRegionName,
                                             const int aPeriod )
{
    mPhysicalDemand[ aPeriod ] = aPhysicalDemand;
    scenario->getMarketplace()->addToDemand( mName, aRegionName,
                                             mPhysicalDemand[ aPeriod ],
                                             aPeriod, true );
}

double BuildingDemandInput::getCoefficient( const int aPeriod ) const {
    assert( mCoefficient.isInited() && mCoefficientAdjustment.isInited() );

    // Return the base coefficient adjusted for calibration.
    return mCoefficient * mCoefficientAdjustment;
}

void BuildingDemandInput::setCoefficient( const double aCoefficient,
                                          const int aPeriod )
{
    // Coefficients must be positive.
    assert( aCoefficient >= 0 && mCoefficient >= 0 );

    // Store the adjusted coefficient locally. Does not modify the base
    // coefficient.
    if( mCoefficient > 0 ){
        mCoefficientAdjustment = aCoefficient / mCoefficient;
    }
}

double BuildingDemandInput::getPrice( const string& aRegionName,
                                      const int aPeriod ) const 
{
    // The market price for building services will be in units of $ per GJ_out.
    // The building demand function needs a service cost in units of $ per EJ_out.
    // However because the service is in billions of m^2 the EJ to GJ will cancel
    // itself out with the billions of m^2 so we do not want to adjust here.
    
    return scenario->getMarketplace()->getPrice( mName, aRegionName, aPeriod );
}

void BuildingDemandInput::setPrice( const string& aRegionName,
                                    const double aPrice,
                                    const int aPeriod )
{
    // Not hooking this up yet, it could work.
}

double BuildingDemandInput::getCalibrationQuantity( const int aPeriod ) const
{
    // Calibration for detailed buildings is a special case where we are
    // calibrating the coef by using upstream and downstream calibrated
    // values.  So to avoid the potential for inconsistencies we are not
    // allowing calibrated values here.
    return -1;
}

bool BuildingDemandInput::hasTypeFlag( const int aTypeFlag ) const {
    return ( ( aTypeFlag & ~IInput::ENERGY ) == 0 );
}

double BuildingDemandInput::getIncomeElasticity() const {
    return 0;
}

double BuildingDemandInput::getPriceElasticity() const {
    return mPriceElasticity;
}

double BuildingDemandInput::getTechChange( const int aPeriod ) const {
    // Building inputs do not currently have input specific technical change.
    return 0;
}

/*!
 * \brief Initialize parameters only known after the parent Technology is
 *        initialized.
 * \param aSectorName Sector name.
 * \param aSubsectorName Subsector name.
 * \param aTechName Technology name.
 * \param aTechInfo Technology info object.
 */
void BuildingDemandInput::initializeParameters( const string& aSectorName,
                                                const string& aSubsectorName,
                                                const string& aTechName,
                                                const IInfo* aTechInfo )
{
    // Store necessary technology values. Generic demand technologies only need
    // the defaults.
    if( mType == eGeneric ){
        mAverageInsulation = mFloorToSurfaceArea = mDegreeDays = 1;
    }
    else {
        mAverageInsulation = aTechInfo->getDouble( "average-insulation", true );
        mFloorToSurfaceArea = aTechInfo->getDouble( "floor-to-surface-area", true );

        // Store the heating or cooling degree days depending on the type.
        if( mType == eHeating ){
            mDegreeDays = aTechInfo->getDouble( "heatingDegreeDays", true );
        }
        else if( mType == eCooling ){
            mDegreeDays = aTechInfo->getDouble( "coolingDegreeDays", true );
        }
        else {
            assert( false );
            mDegreeDays = 1;
        }
    }
}

/*!
 * \brief Calculate an adjustment factor which modifies the input demand as
 *        given by the Technology to a service specific demand.
 * \details TODO
 * \return Input demand adjustment factor to convert from Technology input
 *         demand to service demand.
 */
double BuildingDemandInput::calcDemandAdjustmentFactor() const {
    // Adjust the demand for any potential heating or cooling service
    // parameters. For generic demands these parameters will all be one.
    double demandAdjustment = mSaturation * mAverageInsulation
                              * mFloorToSurfaceArea * mDegreeDays;

    // Don't allow the adjustment factor to be zero.
    return demandAdjustment > util::getSmallNumber() ? demandAdjustment : 1;
}

/*!
 * \brief Return a structure containing the type information.
 * \param aType Enum type.
 * \return The type information.
 */
const BuildingDemandInput::TypeInfo&
BuildingDemandInput::getTypeInfo( const BuildingInputType aType ){
    // There is no type info for the end struct.
    assert( aType != eEnd );

    // Setup the information for each type. The format for the structure is enum
    // type, type string.
    static const TypeInfo types[] = {
                                      { eGeneric, "generic" },
                                      { eHeating, "heating" },
                                      { eCooling, "cooling" }
                                    };
    return types[ aType ];
}

/*!
 * \brief Convert a string to a building input type enum.
 * \param aTypeString The string to convert.
 * \return The enum type, the end marker if the type is unknown.
 */
BuildingDemandInput::BuildingInputType
BuildingDemandInput::convertStringToType( const string& aTypeString ){
    // Iterate over the type list to find the string.
    unsigned int i = 0;
    for( ; i < eEnd; ++i ){
        if( getTypeInfo( BuildingInputType( i ) ).mTypeString == aTypeString ){
            break;
        }
    }
    return BuildingInputType( i );
}

void BuildingDemandInput::copyParamsInto( BuildingDemandInput& aInput,
                                          const int aPeriod ) const
{
    // Check that this is the right input.
    assert( aInput.mName == mName );

    // Only copy the coefficient adjustment factor, not the entire coefficient
    // as that includes components such as saturation that may change over time.
    // do not override a future adjustment if it was already read in
    // we do not need to check calibration conditions since the user will get
    // warned elsewhere if they set up something that will lead to trouble
    if( !aInput.mCoefficientAdjustment.isInited() ) {
        aInput.mCoefficientAdjustment = mCoefficientAdjustment;
    }
}

void BuildingDemandInput::doInterpolations( const int aYear, const int aPreviousYear,
                                            const int aNextYear, const IInput* aPreviousInput,
                                            const IInput* aNextInput )
{
    const BuildingDemandInput* prevBldInput = static_cast<const BuildingDemandInput*>( aPreviousInput );
    const BuildingDemandInput* nextBldInput = static_cast<const BuildingDemandInput*>( aNextInput );
    
    /*!
     * \pre We are given a valid BuildingDemandInput for the previous input.
     */
    assert( prevBldInput );
    
    /*!
     * \pre We are given a valid BuildingDemandInput for the next input.
     */
    assert( nextBldInput );
    
    // interpolate the saturation
    mSaturation = util::linearInterpolateY( aYear, aPreviousYear, aNextYear,
                                            prevBldInput->mSaturation,
                                            nextBldInput->mSaturation );
}
