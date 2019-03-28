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
* \file sgm_output.cpp
* \ingroup Objects
* \brief SGMOutput class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "technologies/include/sgm_output.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/iinfo.h"
#include "util/base/include/xml_helper.h"
#include "functions/include/function_utils.h"
#include "marketplace/include/cached_market.h"

extern Scenario* scenario;

using namespace std;
using namespace xercesc;

SGMOutput::SGMOutput( const string& aSectorName )
    : mName( aSectorName ), mPhysicalOutputs( scenario->getModeltime()->getmaxper() ),
      mCachedCO2Coef( 0 )
{
}

SGMOutput::~SGMOutput() {
}

SGMOutput* SGMOutput::clone() const
{
    // we do not actually want to copy any of the stored
    // output quantaties so just create a new one with the
    // same name
    return new SGMOutput( mName );
}

bool SGMOutput::isSameType( const string& aType ) const
{
    return aType == "sgm-output";
}

const string& SGMOutput::getName() const
{
    // Make sure the name is initialized.
    assert( !mName.empty() );

    return mName;
}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& SGMOutput::getXMLReportingName() const{
    return getXMLReportingNameStatic();
}

const string& SGMOutput::getXMLReportingNameStatic() {
    const static string XML_REPORTING_NAME = "output-SGM";
    return XML_REPORTING_NAME;
}

bool SGMOutput::XMLParse( const DOMNode* aNode )
{
    /*! \pre make sure we were passed a valid node. */
    assert( aNode );

    const Modeltime* modeltime = scenario->getModeltime();

    mName = XMLHelper<string>::getAttr( aNode, "name" );
    // get all child nodes.
    DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if ( nodeName == "physical-output" ) {
            XMLHelper<Value>::insertValueIntoVector( curr, mPhysicalOutputs, modeltime );
        }
        // TODO: do I need a derived class parse?
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " 
                << getXMLReportingName() << "." << endl;
        }
    }
    // TODO: why is there a return value here?
    return true;
}

void SGMOutput::toInputXML( ostream& aOut,
                            Tabs* aTabs ) const
{
    const Modeltime* modeltime = scenario->getModeltime();
    XMLWriteOpeningTag( getXMLReportingName(), aOut, aTabs, mName );
    XMLWriteVector( mPhysicalOutputs, "physical-output", aOut, aTabs, modeltime );
    XMLWriteClosingTag( getXMLReportingName(), aOut, aTabs );
}

void SGMOutput::toDebugXML( const int aPeriod,
                            ostream& aOut,
                            Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLReportingName(), aOut, aTabs, mName );
    XMLWriteElement( mPhysicalOutputs[ aPeriod ], "physical-output", aOut, aTabs );
    XMLWriteElement( mCachedCO2Coef, "cached-co2-coef", aOut, aTabs );
    XMLWriteClosingTag( getXMLReportingName(), aOut, aTabs );
}

void SGMOutput::completeInit( const string& aSectorName,
                              DependencyFinder* aDependencyFinder,
                              const IInfo* aTechInfo,
                              const bool aIsTechOperating )
{
    // SGM outputs do not have any additional dependencies.
}

void SGMOutput::initCalc( const string& aRegionName,
                          const string& aSectorName,
                          const int aPeriod )
{
    // Make sure the sgm output has a name.
    assert( !mName.empty() );

    mCachedMarket = scenario->getMarketplace()->locateMarket( mName, aRegionName, aPeriod );

    // Initialize the cached CO2 coefficient.  If this output IsPrimaryEnergyGood
    // then we want to use a co2 coef of zero since emissions of fuels are accounted
    // for by use instead of production.
    const IInfo* marketInfo = mCachedMarket->getMarketInfo( aSectorName,
        aRegionName, aPeriod, false );
    if( marketInfo && marketInfo->getBoolean( "IsPrimaryEnergyGood", false ) ){
        mCachedCO2Coef.set( 0 );
    }
    else {
        mCachedCO2Coef.set( FunctionUtils::getCO2Coef( aRegionName, mName, aPeriod, false ) );
    }

    mRegionName = aRegionName;
}

void SGMOutput::postCalc( const string& aRegionName,
                          const int aPeriod )
{
}

void SGMOutput::scaleCoefficient( const double aScaler ){
    // SGMOutput outputs do not support scaling.
    assert( false );
}

IOutput::OutputList SGMOutput::calcPhysicalOutput( const double aSGMOutput,
                                                   const string& aRegionName,
                                                   const ICaptureComponent* aCaptureComponent,
                                                   const int aPeriod ) const
{
    assert( false );
    // TODO?
    return OutputList();
}

void SGMOutput::setPhysicalOutput( const double aSGMOutput,
                                   const string& aRegionName,
                                   ICaptureComponent* aCaptureComponent,
                                   const int aPeriod )
{
    assert( aSGMOutput >= 0 );
    
    mPhysicalOutputs[ aPeriod ] = aSGMOutput;
    
    // add the output to the marketplace.
    // note that this does not make sense for consumers
    // and so we say that the market does not need to exist
    if( aSGMOutput > util::getSmallNumber() ) {
        mCachedMarket->addToSupply( mName, aRegionName, aSGMOutput, aPeriod, false );
    }
}

double SGMOutput::getPhysicalOutput( const int aPeriod ) const
{
    return mPhysicalOutputs[ aPeriod ];
}
void SGMOutput::setCurrencyOutput( const std::string& aRegionName,
                                    const double aOutput,
                                    const int aPeriod )
{
    assert( false );
    // could get this to work
}

double SGMOutput::getCurrencyOutput( const int aPeriod ) const {
    /*!
     * \warning Can not use FunctionUtil to get the price since it would not exist
     *          for consumers.
     */
    // we can use the price received which is typically in $ / GJ to convert from
    // physical to currency
    // note that this will get the full currency demand including any production taxes
    // which is why we just the price rather than price recieved
    double priceRecieved = mCachedMarket->getPrice( getName(), mRegionName, aPeriod, false );
    if( priceRecieved == Marketplace::NO_MARKET_PRICE ) {
        priceRecieved = 1;
    }
    
    return getPhysicalOutput( aPeriod ) * priceRecieved;
}

double SGMOutput::getValue( const string& aRegionName,
                            const ICaptureComponent* aCaptureComponent,
                            const int aPeriod ) const
{
    assert( false );
    // The value of the sgm output is fully accounted by the technology.
    return 0;
}

double SGMOutput::getEmissionsPerOutput( const string& aGHGName,
                                         const int aPeriod ) const
{
    /*!
     * \warning Currently assumes that the ghg name will
     *          be CO2 and no other gasses will be supported.
     */
    assert( aGHGName == "CO2" );
    
    assert( mCachedCO2Coef.isInited() );
    return mCachedCO2Coef;
}

void SGMOutput::accept( IVisitor* aVisitor,
                        const int aPeriod ) const
{
    aVisitor->startVisitOutput( this, aPeriod );
    aVisitor->endVisitOutput( this, aPeriod );
}
