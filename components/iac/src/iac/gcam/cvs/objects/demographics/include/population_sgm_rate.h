#ifndef _POPULATION_SGM_RATE_H_
#define _POPULATION_SGM_RATE_H_
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
* \file population_sgm_rate.h
* \ingroup Objects-SGM
* \brief The PopulationSGMRate class header file.
* \author Sonny Kim
* \author Katherine Chung
*/

#include <vector>

#include "demographics/include/population.h"

// forward class declaration
class AgeCohort;

/*! 
* \ingroup Objects-SGM
* \brief An object which contains the PopulationSGMRate information for a region.
* \details PopulationSGMRate reads in population information for the first year
*  then uses survival rates, etc, to calculate the populations for future years.
*/

class PopulationSGMRate : public Population
{
public:
    PopulationSGMRate();
    virtual ~PopulationSGMRate();
    
    virtual Population* cloneAndInterpolate( const int aNewYear, const Population* aNextPopulation ) const;

    virtual void completeInit( const std::vector<double>& femalePopFromPrev = std::vector<double>(),
        const std::vector<double>& malePopFromPrev = std::vector<double>() );
    
    virtual void initCalc();
    static const std::string& getXMLNameStatic();
    
    const std::vector<double> getSurvMalePop() const;
    const std::vector<double> getSurvFemalePop() const;
    
    double getWorkingAgePop() const;
    double getWorkingAgePopMale() const;
    double getWorkingAgePopFemale() const;

    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
protected:
    std::vector<AgeCohort*> ageCohort; //!< array of age cohorts
    typedef std::vector<AgeCohort*>::iterator AgeCohortIterator;
    typedef std::vector<AgeCohort*>::const_iterator CAgeCohortIterator;
    static const std::string XML_NAME; //!< node name for toXML methods

    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string &nodeName, const xercesc::DOMNode* curr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( std::ostream& out, Tabs* tabs ) const;
private:
    void fillMalePop( const std::vector<double>& aMalePop );
    void fillFemalePop( const std::vector<double>& aFemalePop );
    void calcPop();
    void clear();
};

#endif // _POPULATION_SGM_RATE_H_

