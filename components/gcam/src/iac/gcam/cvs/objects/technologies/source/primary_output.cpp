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
* \file primary_output.cpp
* \ingroup Objects
* \brief PrimaryOutput class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "technologies/include/primary_output.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "containers/include/iinfo.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/ivisitor.h"
#include "functions/include/function_utils.h"
#include "marketplace/include/cached_market.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string PrimaryOutput::XML_REPORTING_NAME = "output-primary";

PrimaryOutput::PrimaryOutput( const string& aSectorName )
    : mPhysicalOutputs( scenario->getModeltime()->getmaxper() ), // Name of the primary output is the sector name.
mName( aSectorName )
{
}

/*!
 * \brief Copy constructor.
 * \note This class requires a copy constructor because it has dynamically
 *          allocated memory.
 * \param aPrimaryOutput Primary output from which to copy.
 */
PrimaryOutput::PrimaryOutput( const PrimaryOutput& aPrimaryOutput )
    : mPhysicalOutputs( scenario->getModeltime()->getmaxper() ),
mName( aPrimaryOutput.mName )
{
    // not copying outputs or coefs is this fine?
}

PrimaryOutput::~PrimaryOutput() {
}

PrimaryOutput* PrimaryOutput::clone() const
{
    return new PrimaryOutput( *this );
}

bool PrimaryOutput::isSameType( const string& aType ) const
{
    return aType == "primary-output";
}

const string& PrimaryOutput::getName() const
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
const string& PrimaryOutput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

bool PrimaryOutput::XMLParse( const DOMNode* aNode )
{
    // Primary outputs are created and cannot be parsed.
    assert( false );
    return false;
}

void PrimaryOutput::toInputXML( ostream& aOut,
                                Tabs* aTabs ) const
{
    // Primary outputs are not parsed and so do not write themselves to XML.
}

void PrimaryOutput::toDebugXML( const int aPeriod,
                                ostream& aOut,
                                Tabs* aTabs ) const
{
    XMLWriteOpeningTag( "primary-output", aOut, aTabs, mName );
    XMLWriteElement( mPhysicalOutputs[ aPeriod ], "output", aOut, aTabs );
    XMLWriteElement( mCachedCO2Coef, "cached-co2-coef", aOut, aTabs );
    XMLWriteClosingTag( "primary-output", aOut, aTabs );
}

void PrimaryOutput::completeInit( const string& aSectorName,
                                  DependencyFinder* aDependencyFinder,
                                  const IInfo* aTechInfo,
                                  const bool aIsTechOperating )
{
    // Primary outputs do not have any additional dependencies.
}

void PrimaryOutput::initCalc( const string& aRegionName,
                              const string& aSectorName,
                              const int aPeriod )
{
    // Make sure the primary output has a name.
    assert( !mName.empty() );

    // Initialize the cached CO2 coefficient.
    mCachedCO2Coef.set( FunctionUtils::getCO2Coef( aRegionName, aSectorName, aPeriod ) );
    
    mCachedMarket = scenario->getMarketplace()->locateMarket( mName, aRegionName, aPeriod );
}

void PrimaryOutput::postCalc( const string& aRegionName,
                              const int aPeriod )
{
}

void PrimaryOutput::scaleCoefficient( const double aScaler ){
    // Primary outputs do not support scaling.
}

IOutput::OutputList PrimaryOutput::calcPhysicalOutput( const double aPrimaryOutput,
                                                       const string& aRegionName,
                                                       const ICaptureComponent* aCaptureComponent,
                                                       const int aPeriod ) const
{
    // Output of the primary output is given.
    OutputList outputList;
    outputList.push_back( make_pair( mName, aPrimaryOutput ) );
    return outputList;
}

void PrimaryOutput::setPhysicalOutput( const double aPrimaryOutput,
                                       const string& aRegionName,
                                       ICaptureComponent* aCaptureComponent,
                                       const int aPeriod )
{
    // Primary output cannot be negative.
    assert( aPrimaryOutput >= 0 );

    // Primary output is given by the technology.
    mPhysicalOutputs[ aPeriod ] = aPrimaryOutput;

    // Add the primary output to the marketplace.
    // TODO: Once demand sectors are separated this should be changed.
    if( aPrimaryOutput > 0.0 ){
        mCachedMarket->addToSupply( mName, aRegionName, aPrimaryOutput, aPeriod, false );
    }
}

double PrimaryOutput::getPhysicalOutput( const int aPeriod ) const {
    // TODO: Assert here that the output has been initialized for the current
    //       iteration. Currently this cannot check even that the period
    //       is initialized because the value is fetched to setup the production
    //       state, which may be before it is calculated.
    return mPhysicalOutputs[ aPeriod ];
}

double PrimaryOutput::getValue( const string& aRegionName,
                                const ICaptureComponent* aCaptureComponent,
                                const int aPeriod ) const
{
    // The value of the primary output is fully accounted by the technology.
    return 0;
}

double PrimaryOutput::getEmissionsPerOutput( const string& aGHGName,
                                             const int aPeriod ) const
{
    assert( mCachedCO2Coef.isInited() );
    return mCachedCO2Coef;
}

void PrimaryOutput::accept( IVisitor* aVisitor, const int aPeriod ) const
{
    aVisitor->startVisitOutput( this, aPeriod );
    aVisitor->endVisitOutput( this, aPeriod );
}


