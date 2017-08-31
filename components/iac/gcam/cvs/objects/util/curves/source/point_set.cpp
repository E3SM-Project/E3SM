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
* \file point_set.cpp
* \ingroup Util
* \brief PointSet class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include "util/curves/include/point_set.h"
#include "util/curves/include/explicit_point_set.h"

using namespace std;

const string PointSet::XML_NAME = "PointSet";

//! Constructor
PointSet::PointSet() {
}

//! Destructor
PointSet::~PointSet(){
}

//! Static function to return the name of the XML element associated with this object.
const string& PointSet::getXMLNameStatic(){
    return XML_NAME;
}

//! Function to return the name of the XML element associated with this object.
const string& PointSet::getXMLName() const {
    return XML_NAME;
}

//! Factory method which returns the correct type of PointSet based on the type string.
PointSet* PointSet::getPointSet( const string& type ) {
    PointSet* pointSet;

    if( type == ExplicitPointSet::getXMLNameStatic() ){
        pointSet = new ExplicitPointSet();
    } 
    else {
        cout << "Invalid type of " << getXMLNameStatic() << " requested: " << type << "." << endl;
        pointSet = 0;
    }
    return pointSet;
}
