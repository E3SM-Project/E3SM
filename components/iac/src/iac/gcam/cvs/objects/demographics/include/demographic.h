#ifndef _DEMOGRAPHIC_H_
#define _DEMOGRAPHIC_H_
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
* \file demographic.h
* \ingroup Objects
* \brief The Demographic class header file.
* \author Sonny Kim
* \author Katherine Chung
*/

#include <vector>
#include <map>
#include "util/base/include/ivisitable.h"
#include "util/base/include/iround_trippable.h"
#include "demographics/include/population.h"

/*! 
* \ingroup Objects
* \brief Demographics model that calculates population by gender and age cohort.
*/

class Demographic: public IVisitable, public IRoundTrippable {
    friend class XMLDBOutputter; // For getXMLName()
public:
    Demographic();
    ~Demographic();

    void XMLParse( const xercesc::DOMNode* node );
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    void completeInit();
    void initCalc();

    static const std::string& getXMLNameStatic();

    void calcPop();
   
    double getWorkingAgePopulation( const int period ) const;
    double getWorkingAgePopulationMales( const int per ) const;
    double getWorkingAgePopulationFemales( const int per ) const; 
    double getTotal( const int per ) const;
    const std::vector<double> getTotalPopVec() const;

    void csvOutputFile( const std::string& regionName ) const; 
    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    void dbOutput( const std::string& regionName ) const; 
    void accept( IVisitor* aVisitor, const int aPeriod ) const;

private:
    void clear();
    int convertPeriodToPopulationIndex( int aPeriod ) const;
    const std::string& getXMLName() const;
    
    //! Vector of Population objects by period.
    std::vector<Population*> population;
    
    //! Mapping of year to index in the population vector. The years are stored
    //! as strings to work around a limitation in the XML parsing helper
    //! routines.
    std::map<std::string,int> yearToMapIndex;

    typedef std::map<std::string,int>::const_iterator CYearMapIterator;
    typedef std::vector<Population*>::iterator PopulationIterator;
    typedef std::vector<Population*>::const_iterator CPopulationIterator;

    //! A binary functor to compare Population years so that they can be sorted.
    struct YearComparator : public std::binary_function<Population*, Population*, bool> {
        /*!
        * \brief Determine if the left hand operand is less than the right.
        * \param lhs The left hand side operand for the sort comparison.
        * \param rhs The right hand side operand for the sort comparison.
        * \return True if lhs is less than rhs, false otherwise.
        */
        bool operator()( const Population* lhs, const Population* rhs ){
            return lhs->getYear() < rhs->getYear();
        }
    };

};

#endif // _DEMOGRAPHIC_H_

