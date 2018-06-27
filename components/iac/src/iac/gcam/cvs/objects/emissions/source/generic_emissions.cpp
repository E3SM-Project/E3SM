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
 * \file generic_emissions.h
 * \ingroup Objects
 * \brief GenericEmissions class header file.
 * \author Jim Naslund
 */

#include "util/base/include/definitions.h"

#include "emissions/include/generic_emissions.h"
#include "emissions/include/aemissions_driver.h"
#include "util/base/include/xml_helper.h"
#include "emissions/include/input_driver.h"

using namespace std;
using namespace xercesc;

//! Default Destructor.
GenericEmissions::~GenericEmissions(){
}

//! Clone operator.
GenericEmissions* GenericEmissions::clone() const {
    return new GenericEmissions( *this );
}

const string& GenericEmissions::getXMLNameStatic(){
    static const string XML_NAME = "generic-ghg";
    return XML_NAME;
}

const string& GenericEmissions::getName() const {
    return mName;
}

void GenericEmissions::initCalc( const string& aRegionName,
                                 const IInfo* aLocalInfo,
                                 const int aPeriod )
{
    AComplexEmissions::initCalc( aRegionName, aLocalInfo, aPeriod );
    if( !mEmissionsDriver.get() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "No driver set for generic-ghg: " << getName()
                << "defaulting to an input driver." << endl;
        mEmissionsDriver.reset( new InputDriver );
    }
}

const string& GenericEmissions::getXMLName() const {
    return getXMLNameStatic();
}

// Assign mName to passed in name.
void GenericEmissions::parseName( const string& aNameAttr ){
    mName = aNameAttr;
}

void GenericEmissions::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    AComplexEmissions::toInputXMLDerived( out, tabs );
    // There is no value because the name of the driver is what is signifigant
    XMLWriteElement( "", mEmissionsDriver->getXMLName(),
                     out, tabs );
}

void GenericEmissions::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    AComplexEmissions::toDebugXMLDerived( period, out, tabs );
    // There is no value because the name of the driver is what is signifigant
    XMLWriteElement( "", mEmissionsDriver->getXMLName(),
                     out, tabs );
}
