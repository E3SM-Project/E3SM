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
* \file generic_output.cpp
* \ingroup Objects
* \brief GenericOutput class source file.
* \author Kate Calvin
*/
#include "util/base/include/definitions.h"
#include <string>

#include "util/base/include/xml_helper.h"
#include "technologies/include/generic_output.h"

using namespace std;

// static initialize.
const string GenericOutput::XML_REPORTING_NAME = "output-generic";

GenericOutput::GenericOutput( const string& aSectorName )
    : PrimaryOutput( aSectorName )
{
}

GenericOutput::~GenericOutput() {
}

GenericOutput* GenericOutput::clone() const
{
    return new GenericOutput( *this );
}

bool GenericOutput::isSameType( const string& aType ) const
{
    return aType == "generic-output";
}

/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& GenericOutput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

void GenericOutput::toDebugXML( const int aPeriod,
                                ostream& aOut,
                                Tabs* aTabs ) const
{
    XMLWriteOpeningTag( "generic-output", aOut, aTabs, mName );
    XMLWriteElement( mPhysicalOutputs[ aPeriod ], "output", aOut, aTabs );
    XMLWriteClosingTag( "generic-output", aOut, aTabs );
}

void GenericOutput::initCalc( const string& aRegionName,
                              const string& aSectorName,
                              const int aPeriod )
{
    // Do nothing
}

void GenericOutput::setPhysicalOutput( const double aPrimaryOutput,
                                       const string& aRegionName,
                                       ICaptureComponent* aCaptureComponent,
                                       const int aPeriod )
{
    // Primary output cannot be negative.
    assert( aPrimaryOutput >= 0 );

    // Primary output is given by the technology.
    mPhysicalOutputs[ aPeriod ] = aPrimaryOutput;
}
