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
* \file technology.cpp
* \ingroup Objects
* \brief Technology class source file.
* \author Sonny Kim
*/
// Standard Library headers
#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>

// User headers
#include "util/base/include/configuration.h"
#include "technologies/include/technology.h"
#include "containers/include/scenario.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/model_time.h"
#include "emissions/include/indirect_emiss_coef.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/gdp.h"
#include "util/logger/include/ilogger.h"
#include "technologies/include/icapture_component.h"
#include "technologies/include/capture_component_factory.h"
#include "technologies/include/ishutdown_decider.h"
#include "technologies/include/shutdown_decider_factory.h"
#include "functions/include/iinput.h"
#include "functions/include/non_energy_input.h"
#include "functions/include/input_capital.h"
#include "functions/include/ifunction.h"
#include "technologies/include/input_factory.h"
#include "functions/include/function_manager.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/iinfo.h"
#include "containers/include/info_factory.h"

#include "technologies/include/ioutput.h"
#include "technologies/include/output_factory.h"
#include "technologies/include/primary_output.h"
#include "technologies/include/output_factory.h"

#include "technologies/include/ical_data.h"
#include "technologies/include/iproduction_state.h"
#include "technologies/include/production_state_factory.h"

#include "technologies/include/marginal_profit_calculator.h"
#include "emissions/include/ghg_factory.h"
#include "emissions/include/co2_emissions.h"

// TODO: Factory for cal data objects.
#include "technologies/include/cal_data_output.h"
#include "technologies/include/itechnical_change_calc.h"
#include "technologies/include/standard_technical_change_calc.h"
#include "functions/include/function_utils.h"
#include "marketplace/include/marketplace.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

typedef vector<IOutput*>::iterator OutputIterator;
typedef vector<IOutput*>::const_iterator COutputIterator;
typedef vector<AGHG*>::iterator GHGIterator;
typedef vector<AGHG*>::const_iterator CGHGIterator;

typedef vector<IProductionState*>::iterator ProductionStateIterator;
typedef vector<IProductionState*>::const_iterator CProductionStateIterator;

typedef vector<IShutdownDecider*>::iterator ShutdownDeciderIterator;
typedef vector<IShutdownDecider*>::const_iterator CShutdownDeciderIterator;

typedef vector<IInput*>::const_iterator CInputIterator;

/*! 
 * \brief Constructor.
 * \param aName Technology name.
 * \param aYear Technology year.
 */
Technology::Technology( const string& aName, const int aYear ): mName( aName ), year ( aYear ){
    init();
}

//! Destructor
Technology::~Technology()
{
    clear();
}

//! Copy constructor
Technology::Technology( const Technology& aTech )
    : mName( aTech.mName ), year( aTech.year )
{
    init();
    copy( aTech );
}

//! Assignment operator.
Technology& Technology::operator =( const Technology& techIn )
{
    if( this != &techIn ) {
        // check for self assignment 
        clear();
        init();
        copy( techIn );
    }
    return *this;
}

//! Helper copy function to avoid replicating code.
void Technology::copy( const Technology& techIn ) {
    mName = techIn.mName;
    mLifetimeYears = techIn.mLifetimeYears;
    mShareWeight = techIn.mShareWeight;
    mParsedShareWeight = techIn.mParsedShareWeight;
    mPMultiplier = techIn.mPMultiplier;

    year = techIn.year;
    mCosts = techIn.mCosts;
    mFixedOutput = techIn.mFixedOutput;
    mAlphaZero = techIn.mAlphaZero;

    if( mCaptureComponent.get() ){
        mCaptureComponent.reset( techIn.mCaptureComponent->clone() );
    }

    // Copy the input vector.
    for( vector<IInput*>::const_iterator iter = techIn.mInputs.begin(); iter != techIn.mInputs.end(); ++iter ) {
        mInputs.push_back( ( *iter )->clone() );
    }

    if( techIn.mCaptureComponent.get() ) {
        mCaptureComponent.reset( techIn.mCaptureComponent->clone() );
    }

    if( techIn.mTechChangeCalc.get() ){
        mTechChangeCalc.reset( techIn.mTechChangeCalc->clone() );
    }
    
    for (CGHGIterator iter = techIn.ghg.begin(); iter != techIn.ghg.end(); ++iter) {
        ghg.push_back( (*iter)->clone() );
    }

    for ( CShutdownDeciderIterator iter = techIn.mShutdownDeciders.begin();
        iter != techIn.mShutdownDeciders.end(); ++iter)
    {
        mShutdownDeciders.push_back( (*iter)->clone() );
    }
    
    for( COutputIterator iter = techIn.mOutputs.begin(); iter != techIn.mOutputs.end(); ++iter ) {
        mOutputs.push_back( ( *iter )->clone() );
    }
    
    // copy keywords for reporting as well
    mKeywordMap = techIn.mKeywordMap;
}

//! Clear member variables.
void Technology::clear()
{
    // Delete the Inputs, GHGs and shutdown deciders.
    for( vector<IInput*>::const_iterator iter = mInputs.begin(); iter != mInputs.end(); ++iter ) {
        delete *iter;
    }
    for( GHGIterator iter = ghg.begin(); iter != ghg.end(); ++iter ) {
        delete *iter;
    }
    for( ShutdownDeciderIterator iter = mShutdownDeciders.begin(); iter != mShutdownDeciders.end(); ++iter ) {
        delete *iter;
    }
    for( ProductionStateIterator iter = mProductionState.begin(); iter != mProductionState.end(); ++iter ) {
        delete *iter;
    }
    for( OutputIterator iter = mOutputs.begin(); iter != mOutputs.end(); ++iter ) {
        delete *iter;
    }
}

//! Initialize elemental data members.
void Technology::init()
{
    // This will be reinitialized in completeInit once the technologies start
    // year is known.
    mLifetimeYears = -1;

    const Modeltime* modeltime = scenario->getModeltime();
    mCosts.resize( modeltime->getmaxper(), -1.0 );
    mProductionState.resize( modeltime->getmaxper(), 0 );
    mProductionFunction = 0;
    mPMultiplier = 1;
    mFixedOutput = -1;
    mAlphaZero = 1;
}

bool Technology::isSameType( const string& aType ) const {
    return aType == getXMLName();
}

/*! \brief Default value for mFixedOutput;
* \author Steve Smith
*/
double Technology::getFixedOutputDefault()
{
    return -1.0;
}

bool Technology::XMLParse( const DOMNode* node )
{
    /*! \pre Assume we are passed a valid node. */
    assert( node );

    const DOMNodeList* nodeList = node->getChildNodes();
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ) {
        const DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() != DOMNode::ELEMENT_NODE ) {
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );        
        if( nodeName == "lifetime" ) {
            mLifetimeYears = XMLHelper<int>::getValue( curr );
        }
        else if( nodeName == "share-weight" ) {
            mParsedShareWeight.set( XMLHelper<double>::getValue( curr ) );
        }
        else if( nodeName == "pMultiplier" ){
           mPMultiplier = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "fixedOutput" ) {
            mFixedOutput = XMLHelper<double>::getValue( curr );
        }
        else if( InputFactory::isOfType( nodeName ) ) {
            parseContainerNode( curr, mInputs, InputFactory::create( nodeName ).release() );
        }
        else if( CaptureComponentFactory::isOfType( nodeName ) ) {
            // Check if a new capture component needs to be created because
            // there is not currently one or the current type does not match the
            // new type.
            if( !mCaptureComponent.get() || !mCaptureComponent->isSameType( nodeName ) ) {
                mCaptureComponent = CaptureComponentFactory::create( nodeName );
            }
            mCaptureComponent->XMLParse( curr );
        }
        else if( nodeName == StandardTechnicalChangeCalc::getXMLNameStatic() ){
            if( !mTechChangeCalc.get() || !mTechChangeCalc->isSameType( nodeName ) ) {
                mTechChangeCalc.reset( new StandardTechnicalChangeCalc );
            }
            mTechChangeCalc->XMLParse( curr );
        }
        else if( ShutdownDeciderFactory::isOfType( nodeName ) ) {
            parseContainerNode( curr, mShutdownDeciders, ShutdownDeciderFactory::create( nodeName ).release() );
        }
        else if( GHGFactory::isGHGNode( nodeName ) ) {
            parseContainerNode( curr, ghg, GHGFactory::create( nodeName ).release() );
        }
        else if( nodeName == CalDataOutput::getXMLNameStatic() ) {
            parseSingleNode( curr, mCalValue, new CalDataOutput );
        }
        else if( OutputFactory::isOfType( nodeName ) ) {
            parseContainerNode( curr, mOutputs, OutputFactory::create( nodeName ).release() );
        }
        else if( nodeName == "keyword" ){
            DOMNamedNodeMap* keywordAttributes = curr->getAttributes();
            for( unsigned int attrNum = 0; attrNum < keywordAttributes->getLength(); ++attrNum ) {
                DOMNode* attrTemp = keywordAttributes->item( attrNum );
                mKeywordMap[ XMLHelper<string>::safeTranscode( attrTemp->getNodeName() ) ] = 
                    XMLHelper<string>::safeTranscode( attrTemp->getNodeValue() );
            }
        }
        // parse derived classes
        else if( !XMLDerivedClassParse( nodeName, curr ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName
                    << " found while parsing " << getXMLName() << "." << endl;
        }
    }

    // TODO: Improve error handling.
    return true;
}

/*!
* \brief Complete the initialization of the technology.
* \note This routine is only called once per model run
* \param aSectorName Sector name, also the name of the product.
* \param aDepDefinder Regional dependency finder.
* \param aSubsectorInfo Subsector information object.
* \param aLandAllocator Regional land allocator.
* \author Josh Lurz
* \warning Markets are not necessarily set when completeInit is called
*/
void Technology::completeInit( const string& aRegionName,
                               const string& aSectorName,
                               const string& aSubsectorName,
                               DependencyFinder* aDepFinder,
                               const IInfo* aSubsectorInfo,
                               ILandAllocator* aLandAllocator )
{
    // Inititalize the technology info object
    mTechnologyInfo.reset( InfoFactory::constructInfo( aSubsectorInfo, mName ) );

    /*! \pre There must be at least one input. */
   // assert( !mInputs.empty() ); //sjs remove this for now since ag techs don't have any inputs at present
    // Check for an unset or invalid year.
    if( year == 0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Technology " << mName << " in sector " << aSectorName
            << " has an invalid year attribute." << endl;
    }

    // If the technology did not read a lifetime calculate a default.
    if( mLifetimeYears == -1 ) {
        mLifetimeYears = calcDefaultLifetime();
    }
    
    // Check if both the original MiniCAM non-energy-input and the new input-capital
    // are in the vector.  If so, eliminate the non-energy-input and use input-capital 
    // only so that non-energy costs are not double accounted.
    // Does not check for fixed and variable O&M, however.
    vector<IInput*>::iterator iterNonEnergy = mInputs.end();
    vector<IInput*>::iterator iterCapital = mInputs.end();
	
    // First look for input-capital since most technologies will have non-energy-input.
    for( vector<IInput*>::iterator iter = mInputs.begin(); iter != mInputs.end(); ++iter ) {
        // Cannot use hasTypeFlag() as both have same type.
        if( ( *iter )->isSameType( InputCapital::getXMLNameStatic() ) ){
            iterCapital = iter;
        }
    }
    // Only look for non-energy-input iterator if input-capital iterator is found.
    if( iterCapital != mInputs.end() ){
        for( vector<IInput*>::iterator iter = mInputs.begin(); iter != mInputs.end(); ++iter ) {
            // Cannot use hasTypeFlag() as both have same type.
            if( ( *iter )->isSameType( NonEnergyInput::getXMLNameStatic() ) ){
                iterNonEnergy = iter;
            }
        }
    }
    // If both are found, then eliminate the orginal non-energy-input since the 
    // more detailed levelized capital calculation is intended to be used.
    //if( iterNonEnergy != mInputs.end() && iterCapital != mInputs.end() ){
    // mInputs.erase( iterNonEnergy );
	//}
    
    // Complete the initialization of the inputs. Pass the inputs and outputs
    // the most local info object available.
    const IInfo* localInfo = getTechInfo() != 0 ? getTechInfo() : aSubsectorInfo;
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        mInputs[ i ]->completeInit( aRegionName, aSectorName, aSubsectorName, mName, aDepFinder, localInfo );
    }

    // Price multiplier must be positive.
    if( mPMultiplier < 0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Price Multiplier >pMultiplier< is " << mPMultiplier << ". Must be positive. "
                << "Value reset to 1. " << endl;
        mPMultiplier = 1;
    }

    // Create the primary output for this technology. All technologies will have
    // a primary output. Always insert the primary output at position 0.
    mOutputs.insert( mOutputs.begin(), new PrimaryOutput( aSectorName ) );

    // Check for attempts to calibrate fixed output.
    if( mFixedOutput != getFixedOutputDefault() && mCalValue.get() ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Cannot calibrate a fixed output Technology. Turning off calibration." << endl;
        mCalValue.reset( 0 );
    }
    for( unsigned int i = 0; i < mOutputs.size(); ++i ) {
        // TODO: Change this when dependencies are determined by period.
        mOutputs[ i ]->completeInit( aSectorName, aDepFinder, localInfo, !hasNoInputOrOutput( 0 ) );
    }

    // Initialize the production function. Uses a virtual method so that
    // subclasses can override the production function used.
    mProductionFunction = getProductionFunction();

    /*! \post The production function reference is initialized. */
    assert( mProductionFunction );

    // initialize the working share-weight to the parsed value even
    // if it was not set
    mShareWeight = mParsedShareWeight;

    // TODO: Calibrating to zero does not work correctly so reset the shareweights
    // to zero and remove the calibration input. This could be improved.
    if( mCalValue.get() && mCalValue->getCalOutput() <= 0.0 ){
        mShareWeight = 0;
    }

    // Clear shareweights for fixed output technologies.
    if( mFixedOutput != getFixedOutputDefault() ) {
        mShareWeight = mParsedShareWeight = 0;
    }

    // Accidentally missing CO2 is very easy to do, and would cause big
    // problems. Add it automatically if it does not exist. Warn the user so
    // they remember to add it.
    const string CO2 = "CO2";
    if( util::searchForValue( ghg, CO2 ) == ghg.end() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::DEBUG );
        mainLog << "Adding CO2 to Technology " << mName << " in region " << aRegionName << " in sector " << aSectorName << "." << endl;
        AGHG* CO2Ghg = new CO2Emissions;
        ghg.push_back( CO2Ghg );
    }

    // Initialize the capture component.
    if( mCaptureComponent.get() ) {
        mCaptureComponent->completeInit( aRegionName, aSectorName, aDepFinder );
    }

    // Initialize the technical change calculator.
    if( mTechChangeCalc.get() ){
        mTechChangeCalc->completeInit();
    }

    // Initialize the cal data object.
    if( mCalValue.get() ) {
        mCalValue->completeInit();
    }

    if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
        const Modeltime* modeltime = scenario->getModeltime();
        bool hasCalInput = false;
        const int periodForYear = modeltime->getyr_to_per( year );
        for( CInputIterator it = mInputs.begin(); it != mInputs.end() && !hasCalInput; ++it ) {
            hasCalInput = (*it)->getCalibrationQuantity( periodForYear ) >= 0;
        }
        if( mCalValue.get() && hasCalInput ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            ILogger& calLog = ILogger::getLogger( "calibration_log" );
            calLog.setLevel( ILogger::WARNING );
            mainLog << "Technology " << getName() << " year " << year << "has both calibrated inputs and outputs.";
            mainLog << "Note that we are assuming the coefficient is consistent and the calibrated output takes precedence.";
            calLog << "Technology " << getName() << " year " << year << "has both calibrated inputs and outputs.";
            calLog << "Note that we are assuming the coefficient is consistent and the calibrated output takes precedence.";
        }
    }
}

//! write object to xml output stream
void Technology::toInputXML( ostream& out,
                             Tabs* tabs ) const
{
    XMLWriteOpeningTag( getXMLVintageNameStatic(), out, tabs, "", year );
    
    // write the xml for the class members.
    if( mParsedShareWeight.isInited() ) {
        XMLWriteElement( mParsedShareWeight, "share-weight", out, tabs );
    }
    int defualtLifetimeYears = calcDefaultLifetime();
    XMLWriteElementCheckDefault( mLifetimeYears, "lifetime", out, tabs, defualtLifetimeYears );
    XMLWriteElementCheckDefault( mFixedOutput, "fixedOutput", out, tabs, getFixedOutputDefault() );
    XMLWriteElementCheckDefault( mPMultiplier, "pMultiplier", out, tabs, 1.0 );
    if( !mKeywordMap.empty() ) {
        XMLWriteElementWithAttributes( "", "keyword", out, tabs, mKeywordMap );
    }

    if( mCalValue.get() ) {
        mCalValue->toInputXML( out, tabs );
    }

    if( mCaptureComponent.get() ) {
        mCaptureComponent->toInputXML( out, tabs );
    }

    if( mTechChangeCalc.get() ) {
        mTechChangeCalc->toInputXML( out, tabs );
    }

    for( vector<IShutdownDecider*>::const_iterator i = mShutdownDeciders.begin(); i != mShutdownDeciders.end(); ++i ) {
        ( *i )->toInputXML( out, tabs );
    }

    for( COutputIterator iter = mOutputs.begin(); iter != mOutputs.end(); ++iter ) {
        ( *iter )->toInputXML( out, tabs );
    }
    for( CGHGIterator iter = ghg.begin(); iter != ghg.end(); ++iter ) {
        ( *iter )->toInputXML( out, tabs );
    }
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        mInputs[ i ]->toInputXML( out, tabs );
    }

    // finished writing xml for the class members.
    toInputXMLDerived( out, tabs );
    XMLWriteClosingTag( getXMLVintageNameStatic(), out, tabs );
}

//! write object to xml debugging output stream
void Technology::toDebugXML( const int period,
                             ostream& out,
                             Tabs* tabs ) const
{
    // Only output technologies that are operating.
    if( !mProductionState[ period ]->isOperating() ){
        return;
    }

    XMLWriteOpeningTag( getXMLName(), out, tabs, mName, year );
    // write the xml for the class members.

    XMLWriteElement( mShareWeight, "share-weight", out, tabs );
    XMLWriteElement( mFixedOutput, "fixedOutput", out, tabs );
    XMLWriteElement( mLifetimeYears, "lifetime", out, tabs );
    XMLWriteElement( mAlphaZero, "alpha-zero", out, tabs );
    XMLWriteElement( mCosts[ period ], "cost", out, tabs );
    XMLWriteElement( mPMultiplier, "pMultiplier", out, tabs );
    if( mCalValue.get() ) {
        mCalValue->toDebugXML( out, tabs );
    }

    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        mInputs[ i ]->toDebugXML( period, out, tabs );
    }


    if( mCaptureComponent.get() ) {
        mCaptureComponent->toDebugXML( period, out, tabs );
    }

    if( mTechChangeCalc.get() ) {
        mTechChangeCalc->toDebugXML( period, out, tabs );
    }

    for( CShutdownDeciderIterator i = mShutdownDeciders.begin(); i != mShutdownDeciders.end(); ++i ){
        ( *i )->toDebugXML( period, out, tabs );
    }
    mProductionState[ period ]->toDebugXML( period, out, tabs );

    for( COutputIterator iter = mOutputs.begin(); iter != mOutputs.end(); ++iter ) {
        ( *iter )->toDebugXML( period, out, tabs );
    }

    // write our ghg object, vector is of number of gases
    for( CGHGIterator i = ghg.begin(); i != ghg.end(); i++ ) {
        ( *i )->toDebugXML( period, out, tabs );
    }

    // finished writing xml for the class members.
    toDebugXMLDerived( period, out, tabs );
    XMLWriteClosingTag( getXMLName(), out, tabs );
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
const string& Technology::getXMLVintageNameStatic()
{
    const static string XML_VINTAGE_NAME = "period";
    return XML_VINTAGE_NAME;
}

/*! \brief Perform initializations that only need to be done once per period.
* \param aRegionName Region name.
* \param aSectorName Sector name, also the name of the product.
* \param aSubsectorInfo Parent information container.
* \param aDemographics Regional demographics container.
* \param aPrevPeriodInfo A structure containing information about the same
*        Technology in a previous period. When the Technology has completed
*        accessing the information stored in this structure it should update it
*        with information about the current period so that the Technology in
*        next period may access it.
* \param aPeriod Model period.
*/
void Technology::initCalc( const string& aRegionName,
                           const string& aSectorName,
                           const IInfo* aSubsectorInfo,
                           const Demographic* aDemographics,
                           PreviousPeriodInfo& aPrevPeriodInfo,
                           const int aPeriod )
{
    if( mCalValue.get() ) {
        mCalValue->initCalc( aDemographics, aPeriod );
    }

    if( mCalValue.get() && ( mCalValue->getCalOutput() < 0 ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::DEBUG );
        mainLog << "Negative calibration value for technology " << mName << ". Calibration removed." << endl;
        mCalValue.reset( 0 );
    }

    // Setup the technology production state which represents how the technology
    // decides to produce output.
    setProductionState( aPeriod );

    for( unsigned int i = 0; i < ghg.size(); i++ ) {
        ghg[ i ]->initCalc( aRegionName, getTechInfo(), aPeriod );
    }

    if( !aPrevPeriodInfo.mIsFirstTech && !aPrevPeriodInfo.mInputs ){
        // The first period technology, which is not necessarily in the base year should
        // not have any previous technology information so do not print the warning.
        if( year != scenario->getModeltime()->getper_to_yr( 0 ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog << "Previous period technology from technology " << mName << " in year " << year
                    << " did not pass forward the required information." << endl;
        }
    }

    // Only copy inputs forward in the starting year of the technology.
    else if( !aPrevPeriodInfo.mIsFirstTech && mProductionState[ aPeriod ]->isOperating() &&
             mProductionState[ aPeriod ]->isNewInvestment() ){
        // Copy information from the previous inputs forward.
        FunctionUtils::copyInputParamsForward( *aPrevPeriodInfo.mInputs, mInputs, aPeriod );
    }
    
    // Setup the structure for copying forward with information about the technology in this period.
    aPrevPeriodInfo.mInputs = &mInputs;

    // Initialize the inputs.
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        mInputs[ i ]->initCalc( aRegionName, aSectorName, mProductionState[ aPeriod ]->isNewInvestment(), false, aPeriod );
    }

    if( mCaptureComponent.get() ) {
        mCaptureComponent->adjustInputs( aRegionName, mInputs, aPeriod );
    }

    for( unsigned int i = 0; i < mOutputs.size(); ++i ) {
        mOutputs[ i ]->initCalc( aRegionName, aSectorName, aPeriod );
    }

    // Determine cumulative technical change. Alpha zero defaults to 1.
    if( mTechChangeCalc.get() ){
        mAlphaZero = mTechChangeCalc->calcAndAdjustForTechChange( mInputs,
                     aPrevPeriodInfo, mProductionFunction, aRegionName,
                     aSectorName, aPeriod );
    }

    // If Calibration is Active, reinitialize share weights for calibration.
    if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
        // For new technology vintages up to and including final calibration period.
        const int FinalCalibrationPeriod = scenario->getModeltime()->getFinalCalibrationPeriod();

        if( ( aPeriod <= FinalCalibrationPeriod ) &&
            mProductionState[ aPeriod ]->isNewInvestment() ){
            // Reinitialize share weight to 1 for competing technology with non-zero read-in share weight
            // for calibration periods only.
            if( mShareWeight > 0 && mShareWeight != 1.0 ){
                // Reinitialize to 1 to remove bias, calculate new share weights and
                // normalize in postCalc to anchor to dominant technology.
                mShareWeight = 1.0;
            }

            // If there is a calibration value re-set 0 shareweight to 1 so that calibration 
            // can occur.
            if( getCalibrationOutput( false, "", aPeriod ) > 0.0 
                && mShareWeight == 0 && mFixedOutput == getFixedOutputDefault() ) {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::NOTICE );
                mainLog << "Resetting zero shareweight for Technology " << mName
                    << " in sector " << aSectorName << " in region " << aRegionName
                    << " since calibration value was present." << endl;

                mShareWeight = 1.0;
            }
        }
    }
}

/*!
 * \brief Initializes the production state for the period.
 * \details Sets up the production state in the given period. Responsibility for
 *          selecting the production state is delegated to the
 *          ProductionStateFactory.
 * \param aPeriod Model period.
 */
void Technology::setProductionState( const int aPeriod ){
    // Check that the state for this period has not already been initialized.
    // Note that this is the case when the same scenario is run multiple times
    // for instance when doing the policy cost calculation.  In which case
    // we must delete the memory to avoid a memory leak.
    if( mProductionState[ aPeriod ] ) {
        delete mProductionState[ aPeriod ];
    }
    
    double initialOutput = 0;
    const Modeltime* modeltime = scenario->getModeltime();
    initialOutput = mOutputs[ 0 ]->getPhysicalOutput( modeltime->getyr_to_per( year ) );
    
    mProductionState[ aPeriod ] =
        ProductionStateFactory::create( year, mLifetimeYears, mFixedOutput,
                                        initialOutput, aPeriod ).release();
}

/*!
 * \brief Calculates all technology benefits and costs not accounted for by the
 *        primary output.
 * \details Technologies may contain greenhouse gases and secondary output,
 *          which incur both costs and benefits to the technology. Costs can be
 *          incurred if the emissions are taxed, or if the secondary output has
 *          a cost. Benefits may accrue if a the emissions are negative or if
 *          the secondary output is has a positive value.
 * \author Sonny Kim, Josh Lurz
 * \param aRegionName The region containing this technology.
 * \param aPeriod The period to calculate this value for.
 * \return Total secondary value.
 */
double Technology::calcSecondaryValue( const string& aRegionName,
                                       const int aPeriod ) const
{
    double totalValue = 0;
    // Add all costs from the GHGs.
    for( unsigned int i = 0; i < ghg.size(); ++i ) {
        totalValue -= ghg[ i ]->getGHGValue( aRegionName, mInputs, mOutputs, mCaptureComponent.get(), aPeriod );
    }

    // Add all values from the outputs. The primary output is included in this
    // loop but will have a value of zero.
    for( unsigned int i = 0; i < mOutputs.size(); ++i ) {
        totalValue += mOutputs[ i ]->getValue( aRegionName, mCaptureComponent.get(), aPeriod );
    }
    return totalValue;
}

/*!
* \brief Perform calculations that need to be done after the solution is found
*        for the period.
* \param aRegionName Region name.
* \param aPeriod Model period that has solved.
*/
void Technology::postCalc( const string& aRegionName,
                           const int aPeriod )
{
    if( mProductionState[ aPeriod ]->isOperating() ) {
        for( unsigned int i = 0; i < mOutputs.size(); ++i ) {
            mOutputs[ i ]->postCalc( aRegionName, aPeriod );
        }
    }

    // make sure calibrated share weights get written out
    if( aPeriod <= scenario->getModeltime()->getFinalCalibrationPeriod() &&
        mProductionState[ aPeriod ]->isNewInvestment() ){
            mParsedShareWeight = mShareWeight;
    }
}

/*! \brief This function calculates the sum of the Carbon Values for all GHG's
*          in this Technology.
* \details The function first checks if a carbon tax exists for the Technology,
*          and if it does loops through all GHGs to calculate a sum carbon
*          value. The GHG function which it calls, getGHGValue() calculates the
*          carbon equivalent of all GHG's contained in this Technology. The
*          totalGHGCost attribute of the Technology is then set to this new
*          value.
* \author Sonny Kim, Josh Lurz
* \param aRegionName The region containing this Technology.
* \param aSectorName The sector containing this Technology.
* \param aPeriod The period to calculate this value for.
* \return The total emissions and storage cost of all ghgs.
*/
double Technology::getTotalGHGCost( const string& aRegionName,
                                    const string& aSectorName,
                                    const int aPeriod ) const
{
    double totalGHGCost = 0;
    // totalGHGCost and carbontax must be in same unit as fuel price
    for( unsigned int i = 0; i < ghg.size(); i++ ) {
        totalGHGCost += ghg[ i ]->getGHGValue( aRegionName, mInputs, mOutputs, mCaptureComponent.get(), aPeriod );
    }
    return totalGHGCost;
}

/*!
* \brief Calculate unnormalized technology unnormalized shares.
* \author Sonny Kim, Steve Smith
* \param aRegionName Region name.
* \param aSectorName Sector name, also the name of the product.
* \param aGDP Regional GDP container.
* \param aLogitExp The logit exponent for this technology nest.
* \param aPeriod Model period.
* \todo Check to see if power function for trivial values really wastes time
* \return The Technology share.
*/
double Technology::calcShare( const string& aRegionName,
                              const string& aSectorName,
                              const GDP* aGDP,
                              const double aLogitExp,
                              const int aPeriod ) const
{
    // A Technology which is not operating does not have a share.
    if( !mProductionState[ aPeriod ] || !mProductionState[ aPeriod ]->isOperating() ){
        return 0;
    }

    // Vintages and fixed output technologies should never have a share.
    if( !mProductionState[ aPeriod ]->isNewInvestment() ||
        mFixedOutput != IProductionState::fixedOutputDefault() ){
        return 0;
    }

    // Calculate the new vintage share.
    assert( getCost( aPeriod ) > 0 );
    double share = mShareWeight* pow( getCost( aPeriod ), aLogitExp );

    double fuelPrefElasticity = calcFuelPrefElasticity( aPeriod );
    // This is rarely used, so probably worth it to not to waste cycles on the
    // power function.
    if( fuelPrefElasticity != 0 ) {
        double scaledGdpPerCapita = aGDP->getBestScaledGDPperCap( aPeriod );
        share *= pow( scaledGdpPerCapita, fuelPrefElasticity );
    }

    assert( share >= 0 && util::isValidNumber( share ) );
    return share;
}

/*! \brief Return true if technology is fixed for no output or input
* 
* returns true if this technology is set to never produce output or input
* At present, this can only be guaranteed by assigning a fixedOutput value of zero.
*
* \author Steve Smith
* \return Returns whether this technology will always have no output or input
*/
bool Technology::hasNoInputOrOutput( const int aPeriod ) const
{
    // Technology has zero fixed output if fixed output was read-in as zero.
    // TODO: There are several other ways this could happen.
    return( util::isEqual( mFixedOutput, 0.0 ) );
}

/*!
 * \brief Get the production function to use for this Technology.
 * \details This is a virtual function which gets a pointer to the production
 *          function to use. This function should be used instead of directly
 *          getting the production function from the FunctionManager so that
 *          derived classes can override the function type. The Technology holds
 *          a weak reference to the function, so memory management must be
 *          handled by the FunctionManager or the class that implements this
 *          function.
 * \return A pointer to a production function.
 */
const IFunction* Technology::getProductionFunction() const
{
    return FunctionManager::getFunction( "minicam-leontief" );
}

/*! \brief Return fixed Technology output
* \details Returns the current value of fixed output.
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aHasRequiredInput Whether the technology should check what the required
*        input is.
* \param aRequiredInput The input the technology is required to have if it
*        returns a fixed output value.
* \param aPeriod model period
* \return Value of fixed output for this Technology
* \author Steve Smith
*/
double Technology::getFixedOutput( const string& aRegionName,
                                   const string& aSectorName,
                                   const bool aHasRequiredInput,
                                   const string& aRequiredInput,
                                   const int aPeriod ) const
{
    /*! \pre If the caller requests only output for a specific fixed input, the
       *        specific input name must be passed. 
       */
    assert( !aHasRequiredInput || !aRequiredInput.empty() );
    // Check that a state has been created for the period.
    assert( mProductionState[ aPeriod ] );

    // Construct a marginal profit calculator. This allows the calculation of 
    // marginal profits to be lazy.
    MarginalProfitCalculator marginalProfitCalc( this );
    return mProductionState[ aPeriod ]->calcProduction( aRegionName,
                                                        aSectorName,
                                                        0, // No variable output.
                                                        &marginalProfitCalc,
                                                        1, // Not shutting down any fixed output using
                                                           // the scale factor.
                                                        mShutdownDeciders,
                                                        aPeriod );
}

/*! \brief Calculates fuel input and Technology output.
* \details Adds demands for fuels and ghg emissions to markets in the
*          marketplace.
* \param aRegionName name of the region
* \param aSectorName name of the product for this sector
* \param aVariableDemand Total variable demand for this subsector.
* \param aFixedOutputScaleFactor Scale factor by which to reduce production when
*        fixed output is greater than demand.
* \param aGDP Regional GDP container.
* \param aPeriod Model period.
*/
void Technology::production( const string& aRegionName,
                             const string& aSectorName,
                             double aVariableDemand,
                             double aFixedOutputScaleFactor,
                             const GDP* aGDP,
                             const int aPeriod )
{
    // Can't have a scale factor and positive demand.
    assert( aFixedOutputScaleFactor == 1 || aVariableDemand == 0 );

    // Can't have negative variable demand.
    assert( aVariableDemand >= 0 && util::isValidNumber( aVariableDemand ) );

    // Check for positive variable demand and positive fixed output.
    assert( mFixedOutput == IProductionState::fixedOutputDefault() || util::isEqual( aVariableDemand, 0.0 ) );

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

    // Use the production state to determine output.
    double primaryOutput =
        mProductionState[ aPeriod ]->calcProduction( aRegionName,
                                                     aSectorName,
                                                     aVariableDemand,
                                                     &marginalProfitCalc,
                                                     aFixedOutputScaleFactor,
                                                     mShutdownDeciders,
                                                     aPeriod );

    // Calculate input demand.
    mProductionFunction->calcDemand( mInputs, primaryOutput, aRegionName, aSectorName,
                                     1, aPeriod, 0, mAlphaZero );

    calcEmissionsAndOutputs( aRegionName, primaryOutput, aGDP, aPeriod );
}

/*!
 * \brief Calculate the emissions, primary and secondary outputs for the
 *        Technology.
 * \details Determines the output levels and emissions for the Technology once
 *          the primary output and input quantities are known. Emissions and
 *          outputs are added to the marketplace by the Output and GHG objects.
 * \param aRegionName Region name.
 * \param aPrimaryOutput Primary output quantity.
 * \param aGDP Regional GDP container.
 * \param aPeriod Period.
 */
void Technology::calcEmissionsAndOutputs( const string& aRegionName,
                                          const double aPrimaryOutput,
                                          const GDP* aGDP,
                                          const int aPeriod )
{
    for( unsigned int i = 0; i < mOutputs.size(); ++i ) {
        mOutputs[ i ]->setPhysicalOutput( aPrimaryOutput, aRegionName, mCaptureComponent.get(), aPeriod );
    }

    // calculate emissions for each gas after setting input and output amounts
    for( unsigned int i = 0; i < ghg.size(); ++i ) {
        ghg[ i ]->calcEmission( aRegionName, mInputs, mOutputs, aGDP, mCaptureComponent.get(), aPeriod );
    }
}

//! calculate GHG emissions from Technology use
/* \brief Get a map containing emissions by gas from the Technology.
* \param aGoodName Name of the sector.
* \param aPeriod Period for which to get emissions.
* \return A map of gas name to emissions. There are more gas names than ghgs.
*/
const map<string, double> Technology::getEmissions( const string& aGoodName,
                                                    const int aPeriod ) const
{
    map<string, double> emissions;
    for( unsigned int i = 0; i < ghg.size(); ++i ) {
        // emissions by gas name only
        emissions[ ghg[ i ]->getName() ] = ghg[ i ]->getEmission( aPeriod );
        // emissions by gas and fuel names combined
        // used to calculate emissions by fuel
        for( unsigned int j = 0; j < mInputs.size(); ++j ) {
            // Calculate the per fuel emissions. TODO: This is very difficult to
            // do as different fuels may have different input coefficients and
            // different emissions coefficients. emissions[ghg[i]->getName() +
            // mInputs[ j ]->getName() ] = ghg[i]->getEmission( aPeriod );
        }
        // add sequestered amount to emissions map used to calculate emissions
        // by fuel if there are sequestered emissions.
        if( mCaptureComponent.get() ) {
            emissions[ ghg[ i ]->getName() + "sequestGeologic" ] = mCaptureComponent->
            getSequesteredAmount( ghg[ i ]->getName(), true, aPeriod );
            emissions[ ghg[ i ]->getName() + "sequestNonEngy" ] = mCaptureComponent->
            getSequesteredAmount( ghg[ i ]->getName(), false, aPeriod );
        }
    }

    return emissions;
}

/* \brief Get a map containing emissions by fuel from the Technology.
* \param aGoodName Name of the sector.
* \param aPeriod Period for which to get emissions.
* \return A map of fuel name to emissions.
*/
const map<string, double> Technology::getEmissionsByFuel( const string& aGoodName,
                                                          const int aPeriod ) const
{
    map<string, double> emissionsByFuel;
    for( unsigned int i = 0; i < ghg.size(); ++i ) {
        emissionsByFuel[ ghg[ i ]->getName() ] = ghg[ i ]->getEmissFuel( aPeriod );
    }
    return emissionsByFuel;
}

/*! \brief Returns Technology name
*
* \author Sonny Kim
* \return sector name as a string
*/
const string& Technology::getName() const
{
    return mName;
}

/*! \brief returns share weight for this Technology
*
* \author Steve Smith
* \return share weight
*/
Value Technology::getShareWeight() const
{
    /*! \post Share weight is a valid number and greater than or equal to zero. */
    assert( util::isValidNumber( mShareWeight ) && mShareWeight >= 0 );

    return mShareWeight;
}

/*!
 * \brief Gets the share weight value which was parsed by the user.
 * \return Parsed share weight value.
 */
Value Technology::getParsedShareWeight() const {
    // The parsed share weight does not have to be valid

    return mParsedShareWeight;
}

/*! \brief scale share weight for this Technology
*
* \author Steve Smith
* \param shareWeightValue new value for share weight
*/
void Technology::setShareWeight( double shareWeightValue )
{
    /*! \pre Share weight value is greater than or equal to zero. */
    assert( shareWeightValue >= 0 );

    mShareWeight = shareWeightValue;
}

/*!
 * \brief Returns true if all output is either fixed or calibrated
 * \param aHasRequiredInput Whether the technology should check what the required
 *        input is.
 * \param aRequiredInput The input the technology is required to have if it
 *        returns that is has a fixed output value.
 */
bool Technology::isOutputFixed( const bool aHasRequiredInput,
                                const string& aRequiredInput,
                                const int aPeriod ) const
{
    /*! \pre If the caller requests only output for a specific fixed input, the
       *        specific input name must be passed. 
       */
    assert( !aHasRequiredInput || !aRequiredInput.empty() );

    /*! \pre Ensure that allInputs has not reached technology as an input name. */
    assert( aRequiredInput != "allInputs" );

    // If the technology does not used the required input, it must be fixed for it.
    if( aHasRequiredInput && !hasInput( aRequiredInput ) ) {
        return true;
    }

    // The technology has fixed output if it is an existing vintage, has
    // exogenously specified output, or a zero share weight.
    if( !mProductionState[ aPeriod ] || !mProductionState[ aPeriod ]->isNewInvestment() || mFixedOutput != -1 || mShareWeight == 0 ) {
        return true;
    }

    // If the technology has a calibrated output value than the output is fixed.
    if( mCalValue.get() ) {
        return true;
    }

    // Search the inputs for a fixed input using the required input. If one
    // input is fixed than the technology is fixed.
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        if( !aHasRequiredInput ||
             ( mInputs[ i ]->getName() == aRequiredInput &&
             mInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ) )
        {
            if( mInputs[ i ]->getCalibrationQuantity( aPeriod ) != -1 ) {
                return true;
            }
        }
    }
    // Otherwise output is not fixed.
    return false;
}

/*!\brief Returns a boolean for whether a technology is a fixed output 
 * technology for new investments.
 * \param aPeriod Model period.
 * \author Sonny Kim
 * \return Boolean for whether the new investment technology has a fixed output.
 */
bool Technology::isFixedOutputTechnology( const int aPeriod ) const
{
    // Technology is a fixed output technology if there is an 
    // exogenously specified output for new investments.
    if( mProductionState[ aPeriod ]->isNewInvestment() && 
        mFixedOutput != getFixedOutputDefault() ) {
        return true;
    }
    // Otherwise new investment is not a fixed output.
    return false;
}

/*!
 * \brief Returns true if this Technology is available for production and not
 *          fixed.
 * \details A true value means that this Technology is available to respond to a
 *          demand and vary its output
 * \author Steve Smith
 * \param aPeriod Model Period
 * \return Whether the Technology is available
 * \todo Check that this is correct.
 */
bool Technology::isAvailable( const int aPeriod ) const
{
    // If it is an existing vintage, has exogenously specified ouput, or a
    // sharweight of zero, it is not available for calibration adjustments.
    if( !mProductionState[ aPeriod ]->isNewInvestment() || mFixedOutput != -1 || mShareWeight == 0 ) {
        return false;
    }

    if( mCalValue.get() && mCalValue->getCalOutput() > 0 ) {
        return true;
    }

    // Otherwise search for a calibrated input.
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        if( mInputs[ i ]->getCalibrationQuantity( aPeriod ) > 0 ) {
            return true;
        }
    }

    // There isn't a calibrated input, so this technology cannot be adjusted.
    return false;
}

/*! \brief Returns whether a technology uses a specific input.
* \details Loops through the input set and checks if the input set exists.
* \param aInputName The name of the input.
* \return Whether the technology uses the given input.
*/
bool Technology::hasInput( const string& aInput ) const
{
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        if( mInputs[ i ]->getName() == aInput ) {
            return true;
        }
    }
    // Input was not found.
    return false;
}

/*! \brief Return the Technology's output for a given period.
* \details Return the previously calculated total technology output for the
*          period.
* \pre production has been called for the iteration.
* \sa production
* \param aPeriod The period for which to get output.
* \return The output for the period.
*/
double Technology::getOutput( const int aPeriod ) const
{
    // Primary output is at position zero.
    return mOutputs[ 0 ]->getPhysicalOutput( aPeriod );
}

/*! \brief Return Technology input cost.
* \param aRegionName The region containing the Technology.
* \param aSectorName The sector containing the Technology.
* \param aPeriod Period in which to calculate the input cost.
* \return A calculated input cost for the Technology.
*/
double Technology::getTotalInputCost( const string& aRegionName,
                                      const string& aSectorName,
                                      const int aPeriod ) const
{
    /*! \pre The technology must have a production function. */
    assert( mProductionFunction );
    double cost = mProductionFunction->calcCosts( mInputs, aRegionName,
                                                  mAlphaZero, aPeriod );
    assert( cost >= 0 );
    return cost;
}

/*! \brief Return the total input cost minus any non-energy costs.
* \todo This assumes a leontief production function.
* \param aRegionName The region containing the Technology.
* \param aSectorName The sector containing the Technology.
* \param aPeriod Period in which to calculate the energy cost.
* \return A calculated energy cost for the Technology.
*/
double Technology::getEnergyCost( const string& aRegionName,
                                  const string& aSectorName,
                                  const int aPeriod ) const
{
    // Calculates the energy cost by first calculating the total cost including
    // all inputs and then removing the non-energy costs.
    double cost = getTotalInputCost( aRegionName, aSectorName, aPeriod );

    // Deduct non-energy costs.
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        if( !mInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ) {
            // TODO: Leontief assumption.
            cost -= mInputs[ i ]->getPrice( aRegionName, aPeriod )
                    * mInputs[ i ]->getCoefficient( aPeriod )
                    / mAlphaZero;
        }
    }
    assert( cost >= -util::getSmallNumber() );
    return cost;
}

/*! \brief Sum the total energy input into the Technology.
* \param aPeriod Period.
* \return A calculated energy cost for the Technology.
*/
double Technology::getEnergyInput( const int aPeriod ) const
{
    double totalEnergy = 0;
    // If technology is not operating return with zero total energy.
    if( mProductionState[ aPeriod ]->isOperating() ) {
        for( unsigned int i = 0; i < mInputs.size(); ++i ) {
            if( mInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ) {
                totalEnergy += mInputs[ i ]->getPhysicalDemand( aPeriod );
            }
        }
        /*! \post totalEnergy must still be positive. */
        assert( totalEnergy >= 0 );
    }
    return totalEnergy;
}

/*!
 * \brief Get the calibration output value.
 * \param aHasRequiredInput Whether the technology should check what the
 *        required input is.
 * \param aRequiredInput The input the technology is required to have if it
 *        returns a fixed output value.
 * \param aPeriod Model period.
 * \return Calibrated output.
 */
double Technology::getCalibrationOutput( const bool aHasRequiredInput,
                                         const string& aRequiredInput,
                                         const int aPeriod ) const
{
    /*! \pre If the caller requests only output for a specific fixed input, the
       *        specific input name must be passed. 
       */
    assert( !aHasRequiredInput || ( !aRequiredInput.empty() && aRequiredInput != "allInputs" ) );

    // Check if this is an existing vintage which cannot have a calibration value.
    if( !mProductionState[ aPeriod ]->isNewInvestment() ) {
        return -1;
    }

    // If an input is required and the technology does not have it return early.
    if( aHasRequiredInput && !hasInput( aRequiredInput ) ) {
        return -1;
    }

    // Check if the technology has a calibrated output value.
    if( mCalValue.get() ) {
        return mCalValue->getCalOutput();
    }

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
                                 * mAlphaZero;
                break;
            }
        }
    }
    return totalCalOutput;
}

/*! \brief Check output and all inputs for new investments for any calibrated values.
* \param aPeriod Model period.
* \return Returns true if calibration value is found.
*/
bool Technology::hasCalibratedValue( const int aPeriod ) const {
    bool hasRequiredInputs = false;
    string requiredInputName = "";
    if ( getCalibrationOutput( hasRequiredInputs, requiredInputName, aPeriod ) == -1 ){
        return false;
    }
    return true;
}

/*! \brief Calculate Technology fuel cost and total cost.
* \details This calculates the cost (per unit output) of this specific
*          Technology. The cost includes fuel cost, carbon value, and non-fuel
*          costs. Conversion efficiency, and optional fuel cost and total price
*          multipliers are used if specified.
* \author Sonny Kim, Steve Smith
* \param aRegionName Region name.
* \param aSectorName SectorName
*/
void Technology::calcCost( const string& aRegionName,
                           const string& aSectorName,
                           const int aPeriod )
{
    // If technology is not operating set cost zero.
    if( !mProductionState[ aPeriod ]->isOperating() )
    {
        mCosts[ aPeriod ] = 0;
    }
    else {
        const double techCost = getTotalInputCost( aRegionName, aSectorName, aPeriod )
                                * mPMultiplier -
                                calcSecondaryValue( aRegionName, aPeriod );

        // We want to allow regional and delivered biomass prices to be negative
        if ( ( aSectorName == "regional biomass" || aSectorName == "delivered biomass" ) ){
            mCosts[ aPeriod ] = techCost;
        }
        else {
            // techcost can drift below zero in disequilibrium.
            mCosts[ aPeriod ] = max( techCost, util::getSmallNumber() );
        }
    }

    assert( util::isValidNumber( mCosts[ aPeriod ] ) );
}

/*!
* \brief Get the total cost of the technology for a period.
* \details Returns the previously calculated cost for a period.
* \pre calcCost has been called for the iteration.
* \param aPeriod Model period.
* \return The total Technology cost.
*/
double Technology::getCost( const int aPeriod ) const
{
    // Check that the cost has been calculated for the period. This could still
    // be a stale cost however if the cost has not been calculated for the
    // iteration.
//    assert( mCosts[ aPeriod ] != -1 ); //sjsTEMP
    return mCosts[ aPeriod ];
}

/*!
 * \brief Calculates the average fuel preference elasticity of all inputs.
 * \details Calculates the average fuel preference elasticity of all energy
 *          inputs. This is equal to the sum of the energy coefficients
 *          multiplied by their elasticities divided by the sum of the energy
 *          coefficients.
 * \param aPeriod Model period.
 * \return The average fuel preference elasticity of all inputs.
 */
double Technology::calcFuelPrefElasticity( const int aPeriod ) const
{
    double totalElas = 0;
    double totalEnergyCoefficients = 0;
    for( CInputIterator i = mInputs.begin(); i != mInputs.end(); ++i ){
        if( (*i)->hasTypeFlag( IInput::ENERGY ) ){
            totalElas += (*i)->getCoefficient( aPeriod ) *
                         (*i)->getIncomeElasticity();
            totalEnergyCoefficients += (*i)->getCoefficient( aPeriod );
        }
    }

    // Normalize the elasticity sum.
    if( totalEnergyCoefficients > util::getSmallNumber() ){
        totalElas /= totalEnergyCoefficients;
    }
    return totalElas;
}

/*! \brief Return a vector listing the names of all the GHGs within the Technology.
* \details This function returns all GHG names the Technology contains. It does 
* this by searching the underlying ghgNameMap.
* \author Josh Lurz
* \return A vector of GHG names contained in the Technology.
*/
const vector<string> Technology::getGHGNames() const
{
    vector<string> names;
    for( CGHGIterator i = ghg.begin(); i != ghg.end(); ++i ){
        names.push_back( ( *i )->getName() );
    }
    return names;
}

/*! \brief returns the number of ghg objects.
*
* Calcuation is done using length of GHG string to be consistant with use of ghg names to access GHG information.
*
*
* \author Steve Smith
*/
int Technology::getNumbGHGs()  const {
    vector<string> ghgNames = getGHGNames();
    return static_cast<int>( ghgNames.size() ); 
}


/*! \brief Copies parameters across periods for a specific GHG 
* \param prevGHG Pointer to the previous GHG object that needs to be passed to
*        the corresponding object this period.
* \warning Assumes there is only one GHG object with any given name
*/
void Technology::copyGHGParameters( const AGHG* prevGHG ) {
    bool found = false;
    for( GHGIterator i = ghg.begin(); i != ghg.end() && !found; ++i ){
        if( (*i)->getName() == prevGHG->getName() ){
            ( *i )->copyGHGParameters( prevGHG );
            found = true;
        }
    }
    if( !found ) {
        ghg.push_back( prevGHG->clone() );
    }
}

/*! \brief Returns the pointer to a specific GHG 
* \param aGHGName Name of GHG 
*/
const AGHG* Technology::getGHGPointer( const string& aGHGName ) const {
    for( CGHGIterator i = ghg.begin(); i != ghg.end(); ++i ){
        if( (*i)->getName() == aGHGName ){
            return *i;
        }
    }
    return 0;
}

//! return value for ghg
double Technology::getEmissionsByGas( const string& aGasName,
                                      const int aPeriod ) const
{
    for( unsigned int i = 0; i < ghg.size(); ++i ) {
        if( ghg[ i ]->getName() == aGasName ) {
            return ghg[ i ]->getEmission( aPeriod );
        }
    }
    return 0;
}

/*!
 * \brief Update the subsectors input fuel map with information from the technology.
 * \param aPeriod Model period.
 * \return Mapping of input name to fuel consumption.
 */
const map<string, double> Technology::getFuelMap( const int aPeriod ) const
{
    // Loop over the inputs and record their fuel usage.
    map<string, double> inputMap;
    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        inputMap[ mInputs[ i ]->getName() ] = mInputs[ i ]->getPhysicalDemand( aPeriod );
    }
    return inputMap;
}

/*! \brief Test to see if calibration worked for this Technology.
* \author Josh Lurz
* \param aPeriod The model period.
* \param aCalAccuracy Accuracy (fraction) to check if calibrations are within.
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aSubsectorName Subsector name.
* \param aPrintWarnings Whether to print a warning.
* \return Whether calibration was successful.
*/
bool Technology::isAllCalibrated( const int aPeriod,
                                  double aCalAccuracy,
                                  const string& aRegionName,
                                  const string& aSectorName,
                                  const string& aSubsectorName,
                                  const bool aPrintWarnings ) const
{
    // Check that the period is the new vintage period and does not have any calibrated values.
    if( !mProductionState[ aPeriod ]->isNewInvestment() && !hasCalibratedValue( aPeriod ) ){
        return true;
    }

    // Check if the technology has any calibrated output.
    double calOutput = getCalibrationOutput( false, "", aPeriod );
    if( calOutput < 0 ){
        return true;
    }

    // Compare calibrated to actual output.
    double output = getOutput( aPeriod );
    double relativeDiff;
    double sectorOutput = scenario->getMarketplace()->getSupply( aSectorName, aRegionName, aPeriod );

    // Do not write warning to main log if the calibration value or the
    // relative difference is smaller than the criteria for calibration accuracy.
    if( calOutput > aCalAccuracy ) {
        relativeDiff = fabs( output - calOutput ) / calOutput;
    }
    else {
        // Use absolute accuracy since the calibrated output level is zero.
        relativeDiff = fabs( output - calOutput );
    }
    // Return false (not calibrated) and print warning only if relativeDiff is
    // greater than the calibration accuracy.
    /*!
     * \warning We must also allow calibration values to be off when it is inconsequential
     *          to the sector output.  This is due to the solver not being able to solve
     *          to a tight enough accuracy.  For example global Rice has an order of magnitude
     *          of 10^9 however Rice in Australia_NZ is in the 10^5 range.  So unless we used
     *          a very tight tolerence we probably won't be able to calibrate exactly when scales
     *          are so different.
     */
    if( relativeDiff > aCalAccuracy && ( fabs( output - calOutput ) / sectorOutput ) > aCalAccuracy ) {
        // Print warning then return false.
        if( aPrintWarnings ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog.setf(ios_base::left,ios_base::adjustfield); // left alignment
            mainLog << "Calibration failed by ";
            mainLog.precision(2); // for floating-point
            mainLog.width(4); mainLog << relativeDiff * 100; mainLog << " %";
            mainLog << " Technology: "; mainLog.width(18); mainLog << mName;
            mainLog << " Region: "; mainLog.width(14); mainLog << aRegionName;
            mainLog << " Sector: "; mainLog.width(12); mainLog << aSectorName;
            mainLog << " Subsector: "; mainLog.width(12); mainLog << aSubsectorName;
            mainLog.precision(4); // for floating-point
            mainLog << " Output: "; mainLog.width(8); mainLog << output;
            mainLog << " Calibration: "; mainLog.width(8); mainLog << calOutput;
            mainLog << " relativeDiff: "; mainLog.width(8); mainLog << relativeDiff;
            mainLog << " SectorOutput: "; mainLog.width(8); mainLog << sectorOutput;
            mainLog << " SectorShare: "; mainLog.width(8); mainLog << calOutput/sectorOutput;
            mainLog << endl;
            mainLog.setf(ios_base::fmtflags(0),ios_base::floatfield); //reset to default
        }
        return false;
    }

    // Calibration at the Technology level was successful.
    return true;
}

//! Set the technology year.
void Technology::setYear( const int aYear )
{
    // This is called through parsing, so report an error to the user.
    if( aYear <= 0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid year passed to set year for technology " << mName << "." << endl;
    }
    else {
        year = aYear;
    }
}

/*! \brief Return the marginal revenue for this Technology's output.
* \details The marginal revenue for the Technology is defined as the market
*          price for the good in the given period divided by the price
*          multiplier.
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aPeriod Model period.
* \return The marginal revenue.
*/
double Technology::getMarginalRevenue( const string& aRegionName,
                                       const string& aSectorName,
                                       const int aPeriod ) const
{
    // TODO: Change to true when below is fixed.
    double marginalRevenue = scenario->getMarketplace()->getPrice( aSectorName, aRegionName,
                                                                   aPeriod, false );

    // Demand sectors won't have markets so the price could be wrong here. This
    // will be fixed by splitting demand and supply sectors.
    // TODO: This will prevent vintaging from working for end use sectors.
    // assert( price != Marketplace::NO_MARKET_PRICE );
    if( marginalRevenue == Marketplace::NO_MARKET_PRICE ){
        return 0;
    }
    
    // Adjust for the price multiplier.
    marginalRevenue /= mPMultiplier;

    // Add any value or costs of secondary good.
    marginalRevenue += calcSecondaryValue( aRegionName, aPeriod );

    return marginalRevenue;
}

/*!
 * \brief Return the Technology info object if one exists.
 * \details By default, Technologies do not have an info object. This function
 *          allows derived classes to supply an IInfo object. The Technology
 *          must always access the Info object through this function, and check
 *          for null.
 * \return The IInfo object if one exists.
 */
const IInfo* Technology::getTechInfo() const
{
    return 0;
}

void Technology::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitTechnology( this, aPeriod );
    acceptDerived( aVisitor, aPeriod );

    for( unsigned int i = 0; i < mOutputs.size(); ++i ) {
        mOutputs[ i ]->accept( aVisitor, aPeriod );
    }

    for( unsigned int i = 0; i < ghg.size(); ++i ) {
        ghg[ i ]->accept( aVisitor, aPeriod );
    }

    for( unsigned int i = 0; i < mInputs.size(); ++i ) {
        mInputs[ i ]->accept( aVisitor, aPeriod );
    }

    aVisitor->endVisitTechnology( this, aPeriod );
}

/**
 * \brief Method for visiting derived technologies.
 * \param aVisitor 
 * \param aPeriod 
 */
void Technology::acceptDerived( IVisitor* aVisitor, const int aPeriod ) const {
    // do nothing
}

/*!
 * \brief Call doInterpolations for each object in aInterpolated with the object that
 *        has the same name as in aPrev and aNext.
 * \details Any value in aPrev or aNext which does not have a corresponding object
 *          in the other will be skipped.  Note that since aInterpolated was cloned
 *          from aPrev in the case of a skip it will be left the same as aPrev.
 * \param aInterpolated A vector of pointers to objects which need to be interpolated.
 * \param aPrev A vector of pointers to objects that can be interpolated from.
 * \param aNext A vector of pointers to objects that can be interpolated to.
 * \param aYear The year which needs to be interpolated.
 * \param aPrevYear The year to interpolate from.
 * \param aNextYear The year to interpolate to.
 */
template<class T>
void interpolateChildVector( std::vector<T*> aInterpolated, std::vector<T*> aPrev,
                             std::vector<T*> aNext, const int aYear, const int aPrevYear,
                             const int aNextYear )
{
    // Sort each vector by name to help identify mismatches.
    // Note that the vectors are passed by value on purpose so that these sorts
    // do not mess with ordering.
    util::NameComparator<T> comp;
    sort( aInterpolated.begin(), aInterpolated.end(), comp );
    sort( aPrev.begin(), aPrev.end(), comp );
    sort( aNext.begin(), aNext.end(), comp );
    
    /*!
     * \pre We are assuming that the prev and the current contain the same
     *      elements since the current should have been cloned from the prev.
     */
    assert( aInterpolated.size() == aPrev.size() );
    unsigned int prevIndex = 0;
    unsigned int nextIndex = 0;
    while( prevIndex < aPrev.size() && nextIndex < aNext.size() ) {
        // If the previous name is less than the next then that means the previous
        // had an additional element which must be skipped.
        if( aPrev[ prevIndex ]->getName() < aNext[ nextIndex ]->getName() ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Skipping interpolation for " << aPrev[ prevIndex ]->getName()
                    << " in year " << aYear << " due to mismatch." << endl;
            ++prevIndex;
            continue;
        }
        
        // If the next name is less than the previous then that means the next
        // had an additional element which must be skipped.        
        if( aNext[ nextIndex ]->getName() < aPrev[ prevIndex ]->getName() ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Skipping interpolation for " << aNext[ nextIndex ]->getName()
                    << " in year " << aYear << " due to mismatch." << endl;
            ++nextIndex;
            continue;
        }
        
        // The previous and next match up so this value can be interpolated.
        assert( aInterpolated[ prevIndex ]->getName() == aPrev[ prevIndex ]->getName() );
        aInterpolated[ prevIndex ]->doInterpolations( aYear, aPrevYear, aNextYear,
                                                      aPrev[ prevIndex ], aNext[ nextIndex ] );
        ++prevIndex;
        ++nextIndex;
    }
    
    // additional error checks
    for( ; prevIndex < aPrev.size(); ++prevIndex ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Skipping interpolation for " << aPrev[ prevIndex ]->getName()
                << " in year " << aYear << " due to mismatch." << endl;
    }
    for( ; nextIndex < aNext.size(); ++nextIndex ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Skipping interpolation for " << aNext[ nextIndex ]->getName()
                << " in year " << aYear << " due to mismatch." << endl;
    }
}

/*!
 * \brief Perform any interpolations for any data members that should be interpolated.
 * \details This technology was created by cloning aPrevTech however some of the data
 *          members should be interpolated between the previous value and the next.  All
 *          of those values should be interploated here.
 * \param aPrevTech The previous technology this one was cloned from.
 * \param aNextTech The next parsed technology to which values should be interpolated to.
 */
void Technology::doInterpolations( const Technology* aPrevTech, const Technology* aNextTech ) {
    // interpolate share weights
    // Zero share weights for prevTech or NextTech are valid end points for interpolation.
    // Do not interpolate from a share-weight in a calibration period as that is
    // likely to be replaced during calibration.
    const Modeltime* modeltime  = scenario->getModeltime();
    if( modeltime->getyr_to_per( aPrevTech->year ) > modeltime->getFinalCalibrationPeriod() 
        && aPrevTech->mParsedShareWeight.isInited() && aNextTech->mParsedShareWeight.isInited() )
    {
        mParsedShareWeight = util::linearInterpolateY( year, aPrevTech->year, aNextTech->year,
            aPrevTech->mParsedShareWeight, aNextTech->mParsedShareWeight );
    }
    else if( modeltime->getyr_to_per( aPrevTech->year ) == modeltime->getFinalCalibrationPeriod() )
    {
        // Make sure the share weight gets interpolated from the calibrated value.
        mParsedShareWeight = Value();
    }
    
    // have inputs do any interpolations
    interpolateChildVector( mInputs, aPrevTech->mInputs, aNextTech->mInputs,
                            year, aPrevTech->year, aNextTech->year );
    
    // have outputs do any interpolations
    interpolateChildVector( mOutputs, aPrevTech->mOutputs, aNextTech->mOutputs,
                            year, aPrevTech->year, aNextTech->year );
    
    // have ghgs do any interpolations
    interpolateChildVector( ghg, aPrevTech->ghg, aNextTech->ghg,
                            year, aPrevTech->year, aNextTech->year );

    // only copy fixed output if the technology did not explicitly set a lifetime
    // which indicates that the user intended to exogenously set a output path as
    // apposed to just a single chunk
    if( mFixedOutput != -1 ) {
        mFixedOutput = mLifetimeYears == -1 ?
            util::linearInterpolateY( year, aPrevTech->year, aNextTech->year,
                                      aPrevTech->mFixedOutput, aNextTech->mFixedOutput ) : 0;
    }
    if( aPrevTech->mLifetimeYears != aNextTech->mLifetimeYears ) {
        mLifetimeYears = aNextTech->mLifetimeYears;
    }
}

/*!
 * \brief Calculates what the default lifetime should be if one was not read-in.
 * \details The default will be the number of years it takes to get from this
 *          technology year to the next model year.
 * \return The appropriate default lifetime for this technology.
 */
int Technology::calcDefaultLifetime() const {
    const Modeltime* modeltime = scenario->getModeltime();
    // TODO: worry about non-aligned technologies?
    const int nextTechPeriod = modeltime->getyr_to_per( year ) + 1;
    return nextTechPeriod < modeltime->getmaxper()
        ? modeltime->getper_to_yr( nextTechPeriod ) - year
        : modeltime->gettimestep( modeltime->getmaxper() - 1 );
}
