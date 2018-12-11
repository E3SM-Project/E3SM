#ifndef _XY_DATA_POINT_H_
#define _XY_DATA_POINT_H_
#if defined(_MSC_VER)
#pragma once
#endif

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
* \file xy_data_point.h
* \ingroup Util
* \brief The XYDataPoint class header file.
* \author Josh Lurz
*/

#include <iosfwd>
#include <xercesc/dom/DOMNode.hpp>
#include "util/curves/include/data_point.h"

/*!
* \ingroup Util
* \brief A DataPoint subclass which defines a simple x, y point.
* \author Josh Lurz
*/

class XYDataPoint: public DataPoint {
        friend std::ostream& operator<<( std::ostream& os, const XYDataPoint& dataPoint ){
            dataPoint.print( os );
            return os;
        }
    public:
        XYDataPoint( const double xIn = 0, const double yIn = 0 );
        ~XYDataPoint();
        const std::string& getXMLName() const;
        static const std::string& getXMLNameStatic();
        bool operator==( const XYDataPoint& rhs ) const;
        bool operator!=( const XYDataPoint& rhs ) const;
	    DataPoint* clone() const;
        double getX() const;
        double getY() const;
        void setX( const double xValue );
        void setY( const double yValue );
        void toInputXML( std::ostream& out, Tabs* tabs ) const;
        void XMLParse( const xercesc::DOMNode* node );
        void invertAxises();
    protected:
        static const std::string XML_NAME;  //!< The name of the XML tag associated with this object.
        double x;
        double y;
        void print( std::ostream& out ) const;
    };
#endif // _XY_POINT_H_
