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
 * \file intensity.cpp
 * \ingroup Objects
 * \brief The Intensity class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include "functions/include/intensity.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"

using namespace std;

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& Intensity::getXMLNameStatic() {
    const static string XML_NAME = "coefficient";
    return XML_NAME;
}

/*!
 * \brief Constructor
 * \param aIntensity An intensity.
 */
Intensity::Intensity( const double aIntensity )
: mReadInIntensity( aIntensity ) {
}

Intensity* Intensity::clone() const {
    return new Intensity( *this );
}

bool Intensity::isSameType( const string& aType ) const {
    return aType == getXMLNameStatic();
}

const string& Intensity::getName() const {
    // Coefficients do not have unique names.
    return getXMLNameStatic();
}

void Intensity::toInputXML( ostream& aOut,
                            Tabs* aTabs ) const
{
    XMLWriteElement( mReadInIntensity, getXMLNameStatic(), aOut, aTabs );
}

void Intensity::toDebugXML( const int aPeriod,
                            ostream& aOut,
                            Tabs* aTabs ) const
{
    XMLWriteElement( mReadInIntensity, getXMLNameStatic(), aOut, aTabs );
}

void Intensity::completeInit() {
    // Error check the intensity.
    if( mReadInIntensity < 0 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Resetting invalid intensity to 0." << endl;
        mReadInIntensity = 0;
    }
}

double Intensity::getCoefficient() const {
    return mReadInIntensity;
}
