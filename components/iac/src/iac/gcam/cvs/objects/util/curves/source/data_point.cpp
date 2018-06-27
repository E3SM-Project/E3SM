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
* \file data_point.cpp
* \ingroup Util
* \brief DataPoint class source file.
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <iostream>
#include "util/curves/include/data_point.h"
#include "util/curves/include/xy_data_point.h"

using namespace std;

const string DataPoint::XML_NAME = "DataPoint";

//! Constructor
DataPoint::DataPoint() {
}

//! Destructor
DataPoint::~DataPoint(){
}

//! Returns if two data points are equal.
bool DataPoint::operator==( const DataPoint& rhs ) const {
    return( ( getX() == rhs.getX() ) && ( ( getY() == rhs.getY() ) ) );
}

//! Returns if two data points are not equal.
bool DataPoint::operator!=( const DataPoint& rhs ) const {
    return !( *this == rhs );
}

//! Returns whether one DataPoint is less than another DataPoint, defined by a lesser x value.
bool DataPoint::operator<( const DataPoint& rhs ) const {
    return ( ( this->getX() < rhs.getX() ) || ( ( this->getX() == rhs.getX() ) && ( this->getY() < rhs.getY() ) ) );
}

//! Returns whether one DataPoint is greater than another DataPoint, defined by a greater x value.
bool DataPoint::operator>( const DataPoint& rhs ) const {
    return ( ( this->getX() > rhs.getX() ) || ( ( this->getX() == rhs.getX() ) && ( this->getY() > rhs.getY() ) ) );
}

//! Returns whether one DataPoint is less than or equal to another DataPoint, defined by a lesser or equal x value.
bool DataPoint::operator<=( const DataPoint& rhs ) const {
    return ( ( this->getX() < rhs.getX() ) || ( ( this->getX() == rhs.getX() ) && ( this->getY() <= rhs.getY() ) ) );
}

//! Returns whether one DataPoint is greater than or equal to another DataPoint, defined by a greater or equal x value.
bool DataPoint::operator>=( const DataPoint& rhs ) const {
    return ( ( this->getX() > rhs.getX() ) || ( ( this->getX() == rhs.getX() ) && ( this->getY() >= rhs.getY() ) ) );
}

//! Static function to return the name of the XML element associated with this object.
const string& DataPoint::getXMLNameStatic(){
    return XML_NAME;
}

//! Function to return the name of the XML element associated with this object.
const string& DataPoint::getXMLName() const {
    return XML_NAME;
}

//! Factory method which returns the correct type of DataPoint based on the type string.
DataPoint* DataPoint::getDataPoint( const string& type ) {
    DataPoint* dataPoint;

    if( type == XYDataPoint::getXMLNameStatic() ){
        dataPoint = new XYDataPoint();
    } 
    else {
        cout << "Invalid type of " << getXMLNameStatic() << " requested: " << type << "." << endl;
        dataPoint = 0;
    }
    return dataPoint;
}
