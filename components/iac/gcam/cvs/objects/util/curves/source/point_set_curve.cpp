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
* \file point_set_curve.cpp
* \ingroup Util
* \brief PointSetCurve class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include "util/curves/include/point_set_curve.h"
#include "util/curves/include/point_set.h"
#include "util/curves/include/explicit_point_set.h" // I dont like this. 
#include "util/curves/include/data_point.h"
#include "util/base/include/xml_helper.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <vector>
#include <cassert>
#include <cfloat>

using namespace std;

const string PointSetCurve::XML_NAME = "PointSetCurve";

/*! \brief Constructor
* \warning This curve is takes responsibility for this PointSet one it is constructed. 
* \param pointSetIn The PointSet which defines this curve's data.
*/
PointSetCurve::PointSetCurve( PointSet* pointSetIn ) {
    pointSet = pointSetIn;
}

/*! \brief Constructor which uses a set of y coordinates, a starting X value, and an interval between X values.
* \details This constructor creates a PointSetCurve by creating a point for each input Y value
* with a companion X value starting at the input starting X, and incrementing by the X interval.
* This function uses the default PointSetType.
* \param pointSetType The type of PointSet to use.
* \param dataPointType The type of DataPoint to use. 
* \param yValues The Y Values to use to create the points.
* \param xStart The first X value to use.
* \param xInterval The amount to increment each X value by. 
* \todo This currently is forced to use a ExplicitPointSet.
*/
PointSetCurve::PointSetCurve( const string pointSetType, const string dataPointType, const vector<double> yValues, const double xStart, const double xInterval ) {
    
    // Create the PointSet
    ExplicitPointSet* exPointSet = new ExplicitPointSet();
    
    // Iterate through the yValues and add a point for each.
    typedef vector<double>::const_iterator YValueIter;
    
    double xVal = xStart;
    for( YValueIter yVal = yValues.begin(); yVal != yValues.end(); yVal++ ){
        // Get a new DataPoint
        DataPoint* currPoint = DataPoint::getDataPoint( dataPointType );
        
        // Set the X and Y values
        currPoint->setX( xVal );
        currPoint->setY( *yVal ); // Note: yVal is an iterator so the value must be dereferenced.
        
        // Add the point.
        exPointSet->addPoint( currPoint );

        // Increment the xValue
        xVal += xInterval;
    }
    // Set the pointSet member variable to the newly created ExplicitPointSet.
    pointSet = exPointSet;
}

//! Destructor
PointSetCurve::~PointSetCurve(){
    delete pointSet;
}

//! Copy Constructor
PointSetCurve::PointSetCurve( const PointSetCurve& otherCurve ){
    pointSet = otherCurve.pointSet->clone();
}

//! Return a copy of the Curve.
PointSetCurve* PointSetCurve::clone() const {
    return new PointSetCurve( *this );
}

//! Static function to return the name of the XML tag this objects uses to write itself out.
const string& PointSetCurve::getXMLNameStatic() {
    return XML_NAME;
}

//! Return the name of the XML tag this objects uses to write itself out.
const string& PointSetCurve::getXMLName() const {
    return XML_NAME;
}

//! Get the Y value corresponding to a given X value.
double PointSetCurve::getY( const double xValue ) const {
    double retValue;

    // First check if the point exists.
    if( pointSet->containsX( xValue ) ){
        retValue = pointSet->getY( xValue );
    } else {
        // Need to interpolate.
        double x1 = pointSet->getNearestXBelow( xValue );
        double x2 = pointSet->getNearestXAbove( xValue );

        // First check if both points are invalid. This means are curve has no points.
        if( x1 == -DBL_MAX && x2 == DBL_MAX ){
            retValue = -DBL_MAX;
        }
        else {
            // Check if the point below is invalid.
            if( x1 == -DBL_MAX ){
                // Get the next value above x2
                x1 = pointSet->getNearestXAbove( x2 );

                // Check if that is valid
                if( x1 == -DBL_MAX ){
                    // There is only one valid point.
                    x1 = x2;
                }
            }

            // Check if the point above is invalid
            if( x2 == DBL_MAX ){
                // Get the next value below x1
                x2 = pointSet->getNearestXBelow( x1 );

                // Check if that is valid.
                if( x2 == -DBL_MAX ){
                    // There is only one valid point.
                    x2 = x1;
                }
            }
            retValue = linearInterpolateY( xValue, x1, pointSet->getY( x1 ), x2, pointSet->getY( x2 ) );
        }
    }
    return retValue;
}

//! Get the X value corresponding to a given Y value.
double PointSetCurve::getX( const double yValue ) const {
      double retValue;

    // First check if the point exists.
    if( pointSet->containsY( yValue ) ){
        retValue = pointSet->getX( yValue );
    } else {
        // Need to interpolate.
        double y1 = pointSet->getNearestYBelow( yValue );
        double y2 = pointSet->getNearestYAbove( yValue );

        // First check if both points are invalid. This means are curve has no points.
        if( y1 == -DBL_MAX && y2 == -DBL_MAX ){
            retValue = -DBL_MAX;
        }
        else {
            // Check if the point below is invalid.
            if( y1 == -DBL_MAX ){
                // Get the next value above y2
                y1 = pointSet->getNearestYAbove( y2 );

                // Check if that is valid
                if( y1 == -DBL_MAX ){
                    // There is only one valid point.
                    y1 = y2;
                }
            }

            // Check if the point above is invalid
            if( y2 == -DBL_MAX ){
                // Get the next value below x1
                y2 = pointSet->getNearestXBelow( y1 );

                // Check if that is valid.
                if( y2 == -DBL_MAX ){
                    // There is only one valid point.
                    y2 = y1;
                }
            }
            retValue = linearInterpolateX( yValue, pointSet->getX( y1 ), y1, pointSet->getX( y2 ), y2 );
        }
    }
    return retValue;
}

//! Set the Y value for a point associated with an X value.
bool PointSetCurve::setY( const double xValue, const double yValue ){
    // Need to do more here I think. Add point?
    return pointSet->setY( xValue, yValue );
}

//! Set an X value for a point associated with a Y value.
bool PointSetCurve::setX( const double yValue, const double xValue ){
    // Need to do more here I think. Add point?
    return pointSet->setX( yValue, xValue );
}


//! Calculate the slope between two x coordinates.
double PointSetCurve::getSlope( const double x1, const double x2 ) const {
    return ( ( getY( x2 ) - getY( x1 ) ) / ( x2 - x1 ) );
}

//! Return the maximum X value contained in the underlying PointSet
double PointSetCurve::getMaxX() const {
    return pointSet->getMaxX();
}

//! Return the maximum Y value contained in the underlying PointSet
double PointSetCurve::getMaxY() const {
    return pointSet->getMaxY();
}

//! Return the minimum X value contained in the underlying PointSet
double PointSetCurve::getMinX() const {
    return pointSet->getMinX();
}

//! Return the minimum Y value contained in the underlying PointSet
double PointSetCurve::getMinY() const {
    return pointSet->getMinY();
}

//! Return a vector of pairs of x y coordinates sorted in increasing x order.
PointSet::SortedPairVector PointSetCurve::getSortedPairs( const double lowDomain, const double highDomain, const int minPoints ) const {
    return pointSet->getSortedPairs( lowDomain, highDomain, minPoints );
}

//! Integrate the curve. Currently uses a trapezoidal integration.
double PointSetCurve::getIntegral( const double lowDomain, const double highDomain ) const {
    typedef pair<double,double> PointPair;

    // Get the underlying points. 
    vector<PointPair> sortedPoints = pointSet->getSortedPairs( lowDomain, highDomain );
    double sum = 0;

    if( !sortedPoints.empty() ){
        // If the lowDomain is defined and is not the first point, add a point to represent represent that.
        if( lowDomain != -DBL_MAX && !util::isEqual( sortedPoints[ 0 ].first, lowDomain ) ){
            PointPair pointA( lowDomain, getY( lowDomain ) );
            sortedPoints.insert( sortedPoints.begin(), pointA );
        }

        // If the highDomain is defined and is not the last point, add a point to represent that.
        if( highDomain != DBL_MAX && !util::isEqual( sortedPoints[ sortedPoints.size() - 1 ].first, highDomain ) ){
            PointPair pointB( highDomain, getY( highDomain ) );
            sortedPoints.push_back( pointB );
        }

        // Perform the integration.
        // Note: If there was only one point, and no low domain or high domain, this loop will be skipped.
        // This is correct, the integral of a single point is 0.
        for( unsigned int i = 1; i < sortedPoints.size(); i++ ){
            double height = sortedPoints.at( i ).first - sortedPoints.at( i - 1 ).first;
            sum += 0.5 * height * ( sortedPoints.at( i ).second + sortedPoints.at( i - 1 ).second );
        }
    }
    return sum;

}

/*!
 * \brief Return a discounted sum of the values of the curve back to the year
 *        aLowDomain assuming the X values are years.
 * \details Computes a sum of the discounted y values of the curve from
 *          aLowDomain to aHighDomain inclusive using the given interest rate.
 *          The present value in any year is calculated as:
 *          pv(x, yearsInFuture, discountRate) = x / ( 1 + discountRate )^n
 * \param aLowDomain Year to begin discounting costs from inclusively. This is the base year.
 * \param aHighDomain Last year to add the discounted value for, inclusive.
 * \param aDiscountRate Discount rate to use, should be in the range 0-1.
 */
double PointSetCurve::getDiscountedValue( const double aLowerBound,
                                          const double aUpperBound,
                                          const double aDiscountRate ) const
{
    assert( aLowerBound <= aUpperBound );
    assert( aDiscountRate >= 0 && aDiscountRate <= 1 );

    double sum = 0;
    for( double year = aLowerBound; year <= aUpperBound; ++year ){
        sum += ( getY( year ) / pow( 1 + aDiscountRate, year - aLowerBound ) );
    }

    assert( sum >= 0 && util::isValidNumber( sum ) );
    return sum;
}

//! Print out the curve to an XML File
void PointSetCurve::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    pointSet->toInputXML( out, tabs );
}

//! Parse a curve from a DOM tree.
bool PointSetCurve::XMLParseDerived( const xercesc::DOMNode* node ) {

    bool nodeParsed = false;

    // assume node is valid.
    assert( node );
    
    // Get the name of the node.
    string nodeName = XMLHelper<void>::safeTranscode( node->getNodeName() );

    if ( nodeName == PointSet::getXMLNameStatic() ){
        nodeParsed = true;
        // First clear the existing pointset to prevent a memory leak.
        delete pointSet;
        pointSet = PointSet::getPointSet( XMLHelper<string>::getAttr( node, "type" ) );
        pointSet->XMLParse( node );
    }
    return nodeParsed;
}

//! Swap the X and the Y axises.
void PointSetCurve::invertAxises() {
    swap( xAxisLabel, yAxisLabel );
    swap( xAxisUnits, yAxisUnits );
    pointSet->invertAxises();
}

//! Perform a linear interpolation determining a y value.
double PointSetCurve::linearInterpolateY( const double xVal, const double x1, const double y1, const double x2, const double y2 ) {
   return ( xVal - x1 ) * ( y2 - y1 )/( x2 - x1 ) + y1;
}

//! Perform a linear interpolation determining an x value.
double PointSetCurve::linearInterpolateX( const double yVal, const double x1, const double y1, const double x2, const double y2 ) {
    return ( yVal - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) + x1;
}

//! Determine the slope of a line.
double PointSetCurve::getSlope( const double x1, const double y1, const double x2, const double y2 ) {
    return ( y2 - y1 ) / ( x2 - x1 );
}

//! Determine the y intercept of a line based on a slope and a point.
double PointSetCurve::getYIntercept( const double slope, const double x1, const double y1 ){
    return ( y1 - slope * x1 );
}

//! Determine the x intercept of a line based on a slope and a point.
double PointSetCurve::getXIntercept( const double slope, const double x1, const double y1 ){
    const double yIntercept = getYIntercept( slope, x1, y1 );
    return ( -1 * yIntercept / slope );
}

void PointSetCurve::print( std::ostream& out, const double lowDomain, const double highDomain, const double lowRange, const double highRange, const int minPoints ) const {

    // Get the sorted points in the specified range.
    vector<pair<double,double> > sortedPairs = pointSet->getSortedPairs( lowDomain, highDomain, minPoints );
    typedef vector<pair<double,double> >::iterator pairIter;
    typedef vector<pair<double,double> >::const_iterator constPairIter;

    // Remove any points not in the specified range.
    for( pairIter iter = sortedPairs.begin(); iter != sortedPairs.end(); iter++ ){
        // If it is out of range remove it.
        if( ( iter->second < lowRange ) || ( iter->second > highRange ) ){
            sortedPairs.erase( iter );
        }
    }

    // Now print the labels.
    out << title << " " << name << " " << numericalLabel << endl;
    out << xAxisLabel << "," << yAxisLabel << endl;

    // Print the points
    for( constPairIter iter = sortedPairs.begin(); iter != sortedPairs.end(); iter++ ){
        out << iter->first << "," << iter->second << endl;
    }
    out << endl;
}
