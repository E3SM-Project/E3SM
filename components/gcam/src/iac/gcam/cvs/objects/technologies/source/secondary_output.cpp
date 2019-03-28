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
* \file secondary_output.cpp
* \ingroup Objects
* \brief SecondaryOutput class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "technologies/include/secondary_output.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "containers/include/iinfo.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/dependency_finder.h"
#include "functions/include/function_utils.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string SecondaryOutput::XML_REPORTING_NAME = "output-secondary";

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& SecondaryOutput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

const string& SecondaryOutput::getXMLNameStatic()
{
    const static string XML_NAME = "secondary-output";
    return XML_NAME;
}

SecondaryOutput::SecondaryOutput()
    : mPhysicalOutputs( scenario->getModeltime()->getmaxper() ),
      mPriceMult( 1.0 )
{
}

SecondaryOutput* SecondaryOutput::clone() const
{
    return new SecondaryOutput( *this );
}

bool SecondaryOutput::isSameType( const string& aType ) const
{
    return aType == getXMLNameStatic();
}

const string& SecondaryOutput::getName() const
{
    // Make sure the name is initialized.
    assert( !mName.empty() );

    return mName;
}

void SecondaryOutput::setName( const string& aName )
{
    // Make sure the name is initialized.
    assert( !aName.empty() );

    mName = aName;
}

bool SecondaryOutput::XMLParse( const DOMNode* aNode )
{
    // assume we are passed a valid node.
    assert( aNode );

    // get all the children.
    DOMNodeList* nodeList = aNode->getChildNodes();

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( aNode, "name" );

    for( unsigned int i = 0; i < nodeList->getLength(); ++i ) {
        const DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "output-ratio" ) {
            mOutputRatio.set( XMLHelper<double>::getValue( curr ) );
        }
        else if( nodeName == "pMultiplier" ) {
            mPriceMult = XMLHelper<double>::getValue( curr );;
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLNameStatic() << "." << endl;
        }
    }

    // TODO: Improve error handling.
    return true;
}

void SecondaryOutput::toInputXML( ostream& aOut,
                                  Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mOutputRatio, "output-ratio", aOut, aTabs );
    XMLWriteElement( mPriceMult, "pMultiplier", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void SecondaryOutput::toDebugXML( const int aPeriod,
                                  ostream& aOut,
                                  Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mOutputRatio, "output-ratio", aOut, aTabs );
    XMLWriteElement( mPriceMult, "pMultiplier", aOut, aTabs );
    XMLWriteElement( mPhysicalOutputs[ aPeriod ], "output", aOut, aTabs );
    XMLWriteElement( mCachedCO2Coef, "cached-co2-coef", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void SecondaryOutput::completeInit( const string& aSectorName,
                                    DependencyFinder* aDependencyFinder,
                                    const IInfo* aTechInfo,
                                    const bool aIsTechOperating )
{
    // Secondary output is removed from demand, so add a dependency.
    if( aIsTechOperating ) {
        aDependencyFinder->addDependency( aSectorName, mName );
    }
}

void SecondaryOutput::initCalc( const string& aRegionName,
                                const string& aSectorName,
                                const int aPeriod )
{
    // Initialize the cached CO2 coefficient. Output ratio is determined by the
    // CO2 coefficient and the ratio of output to the primary good.
    const double CO2Coef = FunctionUtils::getCO2Coef( aRegionName, mName, aPeriod );
    mCachedCO2Coef.set( CO2Coef * mOutputRatio );
}


void SecondaryOutput::postCalc( const string& aRegionName,
                          const int aPeriod )
{
}

void SecondaryOutput::scaleCoefficient( const double aScaler ){
    // Secondary outputs do not support scaling.
    // TODO: Should they? This interface is not great.
}

IOutput::OutputList SecondaryOutput::calcPhysicalOutput( const double aPrimaryOutput,
                                                         const string& aRegionName,
                                                         const ICaptureComponent* aCaptureComponent,
                                                         const int aPeriod ) const
{
    OutputList outputList;
    outputList.push_back( make_pair( mName, calcPhysicalOutputInternal( aPrimaryOutput ) ) );
    return outputList;
}

void SecondaryOutput::setPhysicalOutput( const double aPrimaryOutput,
                                         const string& aRegionName,
                                         ICaptureComponent* aCaptureComponent,
                                         const int aPeriod )
{
    // Secondary output is the primary output multiplied by the output ratio.
    mPhysicalOutputs[ aPeriod ] = calcPhysicalOutputInternal( aPrimaryOutput );

    // Remove the secondary output from demand instead of adding to supply
    // because the sector which has this output as a primary will attempt to
    // fill all of demand. If this technology also added to supply, supply would
    // not equal demand.
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->addToDemand( mName, aRegionName, -1 * mPhysicalOutputs[ aPeriod ], aPeriod, true );
}

double SecondaryOutput::getPhysicalOutput( const int aPeriod ) const
{
    assert( mPhysicalOutputs[ aPeriod ].isInited() );
    return mPhysicalOutputs[ aPeriod ];
}

double SecondaryOutput::getValue( const string& aRegionName,
                                  const ICaptureComponent* aCaptureComponent,
                                  const int aPeriod ) const
{
    double price = scenario->getMarketplace()->getPrice( mName, aRegionName, aPeriod, true );

    // Market price should exist or there is not a sector with this good as the
    // primary output. This can be caused by incorrect input files.
    assert( price != Marketplace::NO_MARKET_PRICE );
    if( price == Marketplace::NO_MARKET_PRICE ) {
        return 0;
    }

    // The value of the secondary output is the market price multiplied by the
    // output ratio.
    return price * mOutputRatio * mPriceMult;
}

double SecondaryOutput::getEmissionsPerOutput( const string& aGHGName,
                                               const int aPeriod ) const
{
    // Currently other GHGs do not use output emissions coefficients.
    assert( aGHGName == "CO2" );
    assert( mCachedCO2Coef.isInited() );
    return mCachedCO2Coef;
}

void SecondaryOutput::accept( IVisitor* aVisitor, const int aPeriod ) const
{
    aVisitor->startVisitOutput( this, aPeriod );
    aVisitor->endVisitOutput( this, aPeriod );
}

/*! 
 * \brief Calculate physical output.
 * \details Physical output of the secondary good is equal to the primary output
 *          multiplied by the coefficient.
 * \return Physical output.
 */
double SecondaryOutput::calcPhysicalOutputInternal( const double aPrimaryOutput ) const {
    return aPrimaryOutput * mOutputRatio;
}

void SecondaryOutput::doInterpolations( const int aYear, const int aPreviousYear,
                                        const int aNextYear, const IOutput* aPreviousOutput,
                                        const IOutput* aNextOutput )
{
    // TODO: do we really want to do this?
    const SecondaryOutput* prevSecOutput = static_cast<const SecondaryOutput*>( aPreviousOutput );
    const SecondaryOutput* nextSecOutput = static_cast<const SecondaryOutput*>( aNextOutput );
    
    /*!
     * \pre We are given a valid SecondaryOutput for the previous output.
     */
    assert( prevSecOutput );
    
    /*!
     * \pre We are given a valid SecondaryOutput for the next output.
     */
    assert( nextSecOutput );
    
    // interpolate the output ratio
    mOutputRatio.set( util::linearInterpolateY( aYear, aPreviousYear, aNextYear,
                                                prevSecOutput->mOutputRatio,
                                                nextSecOutput->mOutputRatio ) );
}
