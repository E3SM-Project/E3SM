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
 * \file variable_production_state.cpp
 * \ingroup Objects
 * \brief VariableProductionState class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "technologies/include/variable_production_state.h"
#include "util/base/include/xml_helper.h"

using namespace std;

VariableProductionState::VariableProductionState(){
}

VariableProductionState* VariableProductionState::clone() const {
    return new VariableProductionState( *this );
}

bool VariableProductionState::isSameType( const std::string& aType ) const {
    return aType == getXMLNameStatic();
}

const string& VariableProductionState::getName() const {
    // Use the XML name as the name.
    return getXMLNameStatic();
}

void VariableProductionState::toDebugXML( const int aPeriod,
                                          ostream& aOut,
                                          Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

const string& VariableProductionState::getXMLNameStatic() {
    const static string XML_NAME = "variable-production-state";
    return XML_NAME;
}

double VariableProductionState::calcProduction( const string& aRegionName,
                                                const string& aSectorName,
                                                const double aVariableOutput,
                                                const MarginalProfitCalculator* aMarginalProfitCalc,
                                                const double aFixedOutputScaleFactor,
                                                const vector<IShutdownDecider*>& aShutdownDeciders,
                                                const int aPeriod ) const
{
    assert( aVariableOutput >= 0 );
    // Variable production states produce exactly what the subsector requested.
    return aVariableOutput;
}

void VariableProductionState::setBaseOutput( const double aBaseOutput,
                                             const int aPeriod )
{
    // This is ignored, variable production has no base level output.
}

bool VariableProductionState::isNewInvestment() const {
    return true;
}

bool VariableProductionState::isOperating() const {
    return true;
}
