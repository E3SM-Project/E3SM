#ifndef _AGE_COHORT_H_
#define _AGE_COHORT_H_
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
* \file age_cohort.h
* \ingroup Objects-SGM
* \brief The AgeCohort class header file.
* \author Sonny Kim
* \author Katherine Chung
*/

#include <memory>
#include "util/base/include/iround_trippable.h"
#include "util/base/include/ivisitable.h"

// forward declare headers
class Male;
class Female;
class Gender;
class Tabs;

/*! 
* \ingroup Objects
* \brief A class which stores population information for a more specific age range
*/

class AgeCohort: public IRoundTrippable, IVisitable {
public:
    AgeCohort();
    ~AgeCohort();
    void XMLParse( const xercesc::DOMNode* node );
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( std::ostream& out, Tabs* tabs ) const;
    void completeInit();
    void initCalc();
    
    static const std::string& getXMLNameStatic();

    double calcMaleBirth();
    double calcFemaleBirth();
    double calcSurvMalePop();
    double calcSurvFemalePop();

    void setMalePop( double aMalePopulation );
    void setFemalePop( double aFemalePopulation );
    double getMalePop() const;
    double getFemalePop() const;
    const std::string& getAgeGroup() const;
    int getLowerAgeBound() const;
    int getUpperAgeBound() const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;
protected:
    std::string ageGroup; //!< age group name
    std::auto_ptr<Male> male; //!< male gender object
    std::auto_ptr<Female> female; //!< male gender object
private:
    bool parseGender( xercesc::DOMNode* aNode );
    mutable int mLowerAgeBound;
    mutable int mUpperAgeBound;
};

#endif // _AGE_COHORT_H_
