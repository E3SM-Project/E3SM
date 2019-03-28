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
* \file internal_gains.cpp
* \ingroup Objects
* \brief InternalGains class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "technologies/include/internal_gains.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "containers/include/iinfo.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/ivisitor.h"
#include "containers/include/dependency_finder.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string InternalGains::XML_REPORTING_NAME = "output-internal-gains";

const string& InternalGains::getXMLNameStatic()
{
    const static string XML_NAME = "internal-gains";
    return XML_NAME;
}

InternalGains::InternalGains()
    : mPhysicalOutputs( scenario->getModeltime()->getmaxper() )
{
}

InternalGains* InternalGains::clone() const
{
    return new InternalGains( *this );
}

bool InternalGains::isSameType( const string& aType ) const
{
    return aType == getXMLNameStatic();
}

const string& InternalGains::getName() const
{
    // There can only be one internal gains object per technology so it does not
    // require a name.
    return getXMLNameStatic();
}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& InternalGains::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

bool InternalGains::XMLParse( const DOMNode* aNode )
{
    // assume we are passed a valid node.
    assert( aNode );

    // get all the children.
    DOMNodeList* nodeList = aNode->getChildNodes();

    for( unsigned int i = 0; i < nodeList->getLength(); ++i ) {
        const DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "output-ratio" ) {
            mOutputRatio = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "heating-market" ) {
            mHeatingMarket = XMLHelper<string>::getValue( curr );
        }
        else if( nodeName == "cooling-market" ) {
            mCoolingMarket = XMLHelper<string>::getValue( curr );
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

void InternalGains::toInputXML( ostream& aOut,
                                Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mOutputRatio, "output-ratio", aOut, aTabs );
    XMLWriteElement( mHeatingMarket, "heating-market", aOut, aTabs );
    XMLWriteElement( mCoolingMarket, "cooling-market", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void InternalGains::toDebugXML( const int aPeriod,
                                ostream& aOut,
                                Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mOutputRatio, "output-ratio", aOut, aTabs );
    XMLWriteElement( mHeatingMarket, "heating-market", aOut, aTabs );
    XMLWriteElement( mCoolingMarket, "cooling-market", aOut, aTabs );
    XMLWriteElement( mPhysicalOutputs[ aPeriod ], "output", aOut, aTabs );
    XMLWriteElement( mCoolingFractionOfYearActive, "cooling-fraction-of-year-active", aOut, aTabs );
    XMLWriteElement( mHeatingFractionOfYearActive, "heating-fraction-of-year-active", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void InternalGains::completeInit( const string& aSectorName,
                                  DependencyFinder* aDependencyFinder,
                                  const IInfo* aTechInfo,
                                  const bool aIsTechOperating )
{
    // Internal gains are removed or added to demand, so add a dependency.
    if( aIsTechOperating ) {
        aDependencyFinder->addDependency( aSectorName, mHeatingMarket );
        aDependencyFinder->addDependency( aSectorName, mCoolingMarket );
    }

    // Store the fraction of year active for heating and cooling.
    mHeatingFractionOfYearActive = aTechInfo->getDouble( "heating-fraction-of-year-active", true );
    mCoolingFractionOfYearActive = aTechInfo->getDouble( "cooling-fraction-of-year-active", true );
}

void InternalGains::initCalc( const string& aRegionName,
                              const string& aSectorName,
                              const int aPeriod )
{
}

void InternalGains::postCalc( const string& aRegionName,
                              const int aPeriod )
{
}

void InternalGains::scaleCoefficient( const double aScaler ){
    // InternalGains gains do not support scaling.
    // TODO: Should they?
}

IOutput::OutputList InternalGains::calcPhysicalOutput( const double aPrimaryOutput,
                                                       const string& aRegionName,
                                                       const ICaptureComponent* aCaptureComponent,
                                                       const int aPeriod ) const
{
    // Calculate total internal gains.
    double totalGains = mOutputRatio * aPrimaryOutput;

    OutputList outputList;
    
    // Add the heating internal gain.
    outputList.push_back( make_pair( mHeatingMarket, totalGains * mHeatingFractionOfYearActive ) );

    // Subtract the cooling interal gain.
    outputList.push_back( make_pair( mCoolingMarket, -1 * totalGains * mCoolingFractionOfYearActive ) );

    return outputList;
}

void InternalGains::setPhysicalOutput( const double aPrimaryOutput,
                                       const string& aRegionName,
                                       ICaptureComponent* aCaptureComponent,
                                       const int aPeriod )
{
    // Calculate total internal gains.
    double totalGains = mOutputRatio * aPrimaryOutput;

    mPhysicalOutputs[ aPeriod ] = totalGains;

    Marketplace* marketplace = scenario->getMarketplace();

    // Subtract from the heating market demand. The internal gains is actually
    // supplying heating services, but adding to demand is required due to
    // ordering constraints.
    marketplace->addToDemand( mHeatingMarket, aRegionName,
                              -1 * totalGains * mHeatingFractionOfYearActive,
                              aPeriod, true );

    // Add to the demand for cooling services. Internal gains cause additional
    // cooling services to be required.
    marketplace->addToDemand( mCoolingMarket, aRegionName,
                              totalGains * mCoolingFractionOfYearActive,
                              aPeriod, true );
}

double InternalGains::getPhysicalOutput( const int aPeriod ) const
{
    assert( mPhysicalOutputs[ aPeriod ].isInited() );
    return mPhysicalOutputs[ aPeriod ];
}

double InternalGains::getValue( const string& aRegionName,
                                const ICaptureComponent* aCaptureComponent,
                                const int aPeriod ) const
{
    const Marketplace* marketplace = scenario->getMarketplace();
    double heatingPrice = marketplace->getPrice( mHeatingMarket, aRegionName, aPeriod, true );
    double coolingPrice = marketplace->getPrice( mCoolingMarket, aRegionName, aPeriod, true );

    // The value of internal gains is equal to the ratio of output to the
    // primary good multiplied by the value of each unit of internal gains. The
    // value of each unit of internal gains is equal to the heating service
    // market price multiplied by the amount of internal gains during the
    // heating season(the positive value) minus the cooling price multiplied by
    // the amount of internal gains during the cooling season(the negative
    // value).
    double averagePrice = mOutputRatio *
                          ( heatingPrice * mHeatingFractionOfYearActive -
                            coolingPrice * mCoolingFractionOfYearActive );
    return averagePrice;
}

double InternalGains::getEmissionsPerOutput( const string& aGHGName,
                                             const int aPeriod ) const
{
    // Internal gains do not contain any emissions
    return 0;
}

void InternalGains::accept( IVisitor* aVisitor,
                            const int aPeriod ) const
{
    aVisitor->startVisitOutput( this, aPeriod );
    aVisitor->endVisitOutput( this, aPeriod );
}

void InternalGains::doInterpolations( const int aYear, const int aPreviousYear,
                                      const int aNextYear, const IOutput* aPreviousOutput,
                                      const IOutput* aNextOutput )
{
    // TODO: do we really want to do this?
    const InternalGains* prevInternalGains = static_cast<const InternalGains*>( aPreviousOutput );
    const InternalGains* nextInternalGains = static_cast<const InternalGains*>( aNextOutput );
    
    /*!
     * \pre We are given a valid InternalGains for the previous output.
     */
    assert( prevInternalGains );
    
    /*!
     * \pre We are given a valid InternalGains for the next output.
     */
    assert( nextInternalGains );
    
    // interpolate the output ratio
    mOutputRatio.set( util::linearInterpolateY( aYear, aPreviousYear, aNextYear,
                                                prevInternalGains->mOutputRatio,
                                                nextInternalGains->mOutputRatio ) );
}
