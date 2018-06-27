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
 * \file byProduct.cpp
 * \ingroup Objects
 * \brief ByProduct class source file.
 * \author Sonny Kim
 */

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "emissions/include/by_product.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "technologies/include/icapture_component.h"
#include "util/base/include/model_time.h"
#include "containers/include/iinfo.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string ByProduct::XML_REPORTING_NAME = "output-by-product";

ByProduct::ByProduct()
    : mQuantities( scenario->getModeltime()->getmaxper() )
{
}

ByProduct* ByProduct::clone() const
{
    return new ByProduct( *this );
}

bool ByProduct::isSameType( const string& aType ) const
{
    return aType == getXMLNameStatic();
}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& ByProduct::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}


bool ByProduct::XMLParse( const DOMNode* aNode )
{
    /*! \pre Assume we are passed a valid node. */
    assert( aNode );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( aNode, "name" );

    DOMNodeList* nodeList = aNode->getChildNodes();
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ) {
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );      

        if( nodeName == "#text" ) {
            continue;
        } else if( nodeName == "byProductCoef" ) {
            mCoef = XMLHelper<double>::getValue( curr );
        } else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLNameStatic() << "." << endl;
        }
    }
    // TODO: Improve error checking.
    return true;
}

void ByProduct::toInputXML( ostream& aOut,
                            Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElementCheckDefault( mCoef.get(), "byProductCoef", aOut, aTabs, 0.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void ByProduct::toDebugXML( const int aPeriod,
                            ostream& aOut,
                            Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mQuantities[ aPeriod ], "quantity", aOut, aTabs );
    XMLWriteElement( mCoef, "byProductCoef", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz 
* \return The constant XML_NAME as a static.
*/
const string& ByProduct::getXMLNameStatic()
{
    const static string XML_NAME = "byProduct";
    return XML_NAME;
}

void ByProduct::completeInit( const std::string& aSectorName,
                              DependencyFinder* aDependencyFinder,
                              const IInfo* aTechInfo,
                              const bool aIsTechOperating )
{
    // Byproducts do not have a dependency.
}

void ByProduct::initCalc( const string& aRegionName,
                          const  string& aSectorName,
                          const int aPeriod )
{
}


void ByProduct::postCalc( const string& aRegionName,
                          const int aPeriod )
{
    // Add byproduct quantity to the market info. 
    Marketplace* marketplace = scenario->getMarketplace();
    IInfo* byProductInfo = marketplace->getMarketInfo( mName, aRegionName, aPeriod, true );

    // Any output is in addition to the current level within the market info.
    double accumRscValue = byProductInfo->getDouble( "AccumulatedRsc", false ); 
    accumRscValue += mQuantities[ aPeriod ];
    byProductInfo->setDouble( "AccumulatedRsc", accumRscValue );
}

void ByProduct::scaleCoefficient( const double aScaler ){
    mAdjustedCoef = mCoef * aScaler;
}

double ByProduct::getValue( const string& aRegionName,
                            const ICaptureComponent* aSequestrationDevice,
                            const int aPeriod ) const
{
    // Convert byProduct tax and any storage costs into energy units using
    // by-product coefficients and return the value or cost of the tax and
    // storage for the by-product. Determine the applicable tax for the
    // by-product.
    // TODO: This is not the right name for the biproduct tax market.
    // const Marketplace* marketplace = scenario->getMarketplace();
    // double byProductTax = marketplace->getPrice( mName, aRegionName, aPeriod, false );
    // if( byProductTax == Marketplace::NO_MARKET_PRICE ) {
    //     byProductTax = 0;
    // }
    double byProductTax = 0;
    // Get the storage cost from the sequestrion device if there is one.
    double storageCost = aSequestrationDevice ? aSequestrationDevice->getStorageCost( mName, aRegionName, aPeriod ) : 0;

    // Get the remove fraction from the sequestration device. The remove
    // fraction is zero if there is no sequestration device.
    double removeFraction = aSequestrationDevice ? aSequestrationDevice->getRemoveFraction( mName ) : 0;

    // units for generalized cost is in 75$/gj
    const double CVRT90 = 2.212; // 1975 $ to 1990 $
    double generalizedCost = ( ( 1.0 - removeFraction ) * byProductTax + removeFraction * storageCost ) * mCoef / CVRT90;

    // Byproducts have a negative value due to the tax and sequestration costs.
    return -1 * generalizedCost;
}

IOutput::OutputList ByProduct::calcPhysicalOutput( const double aPrimaryOutput,
                                                   const string& aRegionName,
                                                   const ICaptureComponent* aCaptureComponent,
                                                   const int aPeriod ) const
{
    OutputList outputList;
    outputList.push_back( make_pair( mName, calcPhysicalOutputInternal( aPrimaryOutput, 
                                                                        aCaptureComponent ) ) );
    return outputList;
}

void ByProduct::setPhysicalOutput( const double aPrimaryOutput,
                                   const string& aRegionName,
                                   ICaptureComponent* aCaptureComponent,
                                   const int aPeriod )
{
    // Primary output is given by the technology.
    mQuantities[ aPeriod ] = calcPhysicalOutputInternal( aPrimaryOutput, aCaptureComponent );

    // Set the captured amount.
    if( aCaptureComponent ) {
        aCaptureComponent->calcSequesteredAmount( aRegionName, mName, aPrimaryOutput, aPeriod );
    }
}

const string& ByProduct::getName() const
{
    return mName;
}

double ByProduct::getPhysicalOutput( const int aPeriod ) const
{
    return mQuantities[ aPeriod ];
}

double ByProduct::getEmissionsPerOutput( const string& aGHGName,
                                         const int aPeriod ) const
{
    // Byproducts do not currently have emissions.
    return 0;
}

void ByProduct::accept( IVisitor* aVisitor,
                        const int aPeriod ) const
{
    aVisitor->startVisitOutput( this, aPeriod );
    aVisitor->endVisitOutput( this, aPeriod );
}

/*! 
 * \brief Calculate the physical output of the byproduct.
 * \details The physical output of the byproduct is equal to the primary output
 *          multiple by the coefficient minus any sequestered output as
 *          determined by the removal fraction of the capture component.
 * \todo When is a capture component useful here?
 * \param aPrimaryOutput Primary output level.
 * \param aCaptureComponent The capture component.
 * \return Output of the byproduct.
 */
double ByProduct::calcPhysicalOutputInternal( const double aPrimaryOutput,
                                              const ICaptureComponent* aCaptureComponent ) const
{
    // The byproduct from a technology is calculated by the byproduct
    // coefficient and a quantity passed in from the technology. This quantity
    // is an output energy amount of the technology.

    // The coefficient is a read in value. Determine the remove fraction.
    double removeFraction = aCaptureComponent ? aCaptureComponent->getRemoveFraction( mName ) : 0;

    // TODO: Fix sequestration device to have seperate calculate and set methods
    // for const correctness. This calculation may be wrong.
    return ( 1.0 - removeFraction ) * ( aPrimaryOutput * mCoef );
}

void ByProduct::doInterpolations( const int aYear, const int aPreviousYear,
                                  const int aNextYear, const IOutput* aPreviousOutput,
                                  const IOutput* aNextOutput )
{
    // TODO: do we really want to do this?
    const ByProduct* prevByProduct = static_cast<const ByProduct*>( aPreviousOutput );
    const ByProduct* nextByProduct = static_cast<const ByProduct*>( aNextOutput );
    
    /*!
     * \pre We are given a valid ByProduct for the previous output.
     */
    assert( prevByProduct );
    
    /*!
     * \pre We are given a valid ByProduct for the next output.
     */
    assert( nextByProduct );
    
    // interpolate the coefficient
    mCoef.set( util::linearInterpolateY( aYear, aPreviousYear, aNextYear,
                                         prevByProduct->mCoef, nextByProduct->mCoef ) );
}

