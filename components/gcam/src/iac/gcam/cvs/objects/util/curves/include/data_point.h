#ifndef _DATA_POINT_H_
#define _DATA_POINT_H_
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
* \file data_point.h
* \ingroup Util
* \brief The DataPoint class header file.
* \author Josh Lurz
*/

#include <iosfwd>
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <functional>

class Tabs;

/*!
* \ingroup Util
* \brief An abstract base class which defines a single point. 
* \author Josh Lurz
*/

class DataPoint {
        friend std::ostream& operator<<( std::ostream& os, const DataPoint& dataPoint ){
            dataPoint.print( os );
            return os;
        }
    public:
        DataPoint();
        virtual ~DataPoint();
        static const std::string& getXMLNameStatic();
        virtual const std::string& getXMLName() const;
        static DataPoint* getDataPoint( const std::string& type );
        virtual bool operator==( const DataPoint& rhs ) const;
        virtual bool operator!=( const DataPoint& rhs ) const;
        virtual bool operator<( const DataPoint& rhs ) const;
        virtual bool operator>( const DataPoint& rhs ) const;
        virtual bool operator<=( const DataPoint& rhs ) const;
        virtual bool operator>=( const DataPoint& rhs ) const;
        virtual DataPoint* clone() const = 0;
        virtual double getX() const = 0;
        virtual double getY() const = 0;
        virtual void setX( const double xValue ) = 0;
        virtual void setY( const double yValue ) = 0;
        virtual void toInputXML( std::ostream& out, Tabs* tabs ) const = 0;
        virtual void XMLParse( const xercesc::DOMNode* node ) = 0;
        virtual void invertAxises() = 0;
        
        /*!
        * \brief Binary comparison operator used for DataPoint pointers to order by increasing values. 
        * \author Josh Lurz
        */  
        struct Lesser : public std::binary_function<DataPoint*,DataPoint*,bool>
        {
            //! Operator which performs comparison. 
            bool operator()( const DataPoint* lhs, const DataPoint* rhs ) const
            {   
                return ( ( lhs->getX() < rhs->getX() ) || ( ( lhs->getX() == rhs->getX() ) && ( lhs->getY() < rhs->getY() ) ) );
            }
        };
        /*!
        * \brief Binary comparison operator which compares two DataPoints by least X value.
        * \author Josh Lurz
        */  
        struct LesserX : public std::binary_function<DataPoint*,DataPoint*,bool>
        {
            //! Operator which performs comparison. 
            bool operator()( const DataPoint* lhs, const DataPoint* rhs ) const
            {   
                return ( lhs->getX() < rhs->getX() );
            }
        };
        /*!
        * \brief Binary comparison operator which compares two DataPoints by least Y value.
        * \author Josh Lurz
        */  
        struct LesserY : public std::binary_function<DataPoint*,DataPoint*,bool>
        {
            //! Operator which performs comparison. 
            bool operator()( const DataPoint* lhs, const DataPoint* rhs ) const
            {   
                return ( lhs->getY() < rhs->getY() );
            }
        };

    protected:
        static const std::string XML_NAME; //!< The name of the XML tag associated with this object.
        virtual void print( std::ostream& out ) const = 0;
    };
#endif // _DATA_POINT_H_
