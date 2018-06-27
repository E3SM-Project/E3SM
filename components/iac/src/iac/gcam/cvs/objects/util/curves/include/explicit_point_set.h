#ifndef _EXPLICIT_POINT_SET_H_
#define _EXPLICIT_POINT_SET_H_
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
* \file explicit_point_set.h
* \ingroup Util
* \brief The ExplicitPointSet class header file.
* \author Josh Lurz
* \todo Add an iterator??
*/

#include <vector>
#include <iosfwd>
#include <string>
#include <cfloat>
#include <xercesc/dom/DOMNode.hpp>
#include "util/curves/include/point_set.h"

class Tabs;
class DataPoint;

/*!
* \ingroup Util
* \brief A PointSet subclass which can be described as a set of points. 
* \author Josh Lurz
*/

class ExplicitPointSet: public PointSet {
    friend std::ostream& operator<<( std::ostream& os, const ExplicitPointSet& pointSet ) {
        pointSet.print( os );
        return os;
    }
public:
    ExplicitPointSet();
    ExplicitPointSet( const ExplicitPointSet& rhs );
    ~ExplicitPointSet();
    ExplicitPointSet& operator=( const ExplicitPointSet& rhs );
    bool operator==( const ExplicitPointSet& rhs ) const;
    bool operator!=( const ExplicitPointSet& rhs ) const;
    ExplicitPointSet* clone() const;
    static const std::string& getXMLNameStatic();
    bool addPoint( DataPoint* dataPoint );
    double getY( const double xValue ) const;
    double getX( const double yValue ) const;
    bool setY( const double xValue, const double yValue );
    bool setX( const double yValue, const double xValue );
    bool removePointFindX( const double xValue );
    bool removePointFindY( const double yValue );
    double getMaxX() const;
    double getMaxY() const;
    double getMinX() const;
    double getMinY() const;
    std::vector<std::pair<double,double> > getSortedPairs( const double lowDomain = -DBL_MAX, const double highDomain = DBL_MAX, const int minPoints = 0 ) const;
    bool containsX( const double x ) const;
    bool containsY( const double y ) const;
    double getNearestXBelow( const double x ) const;
    double getNearestXAbove( const double x ) const;
    double getNearestYBelow( const double x ) const;
    double getNearestYAbove( const double x ) const;
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void XMLParse( const xercesc::DOMNode* node );
    void invertAxises();
protected:   
    static const std::string XML_NAME; //!< The name of the XML tag associated with this object.
    std::vector<DataPoint*> points;
    typedef std::vector<DataPoint*>::iterator DataPointIterator;
    typedef std::vector<DataPoint*>::const_iterator DataPointConstIterator;
	const std::string& getXMLName() const;
	void copy( const ExplicitPointSet& rhs );
    void clear();
	const DataPoint* findX( const double xValue ) const;
    DataPoint* findX( const double xValue );
    const DataPoint* findY( const double yValue ) const;
    DataPoint* findY( const double yValue );
    void print( std::ostream& out, const double lowDomain = -DBL_MAX, const double highDomain = DBL_MAX,
        const double lowRange = -DBL_MAX, const double highRange = DBL_MAX, const int minPoints = 0 ) const;
};
#endif // _EXPLICIT_POINT_SET_H_
