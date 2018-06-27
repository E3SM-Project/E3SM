#ifndef _MODEL_TIME_H_
#define _MODEL_TIME_H_
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
* \file model_time.h
* \ingroup Objects
* \brief The Modeltime class header file.
* \author Sonny Kim
*/

#include <vector>
#include <map>
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"
/*! 
* \ingroup Objects
* \brief A class which defines the time information necessary for the model to run.
* \todo This class needs to be cleaned up and documented. 
* \author Sonny Kim
*/

class Modeltime: public IParsable, public IRoundTrippable
{
private:
    //! Model start year (read-in).
    int mStartYear;

    //! Model end year (read-in).
    int mEndYear;
    
    //! The final year in which calibration occurs.
    int mFinalCalibrationYear;

    //! Maximum number of model periods (calculated).
    int mMaxPeriod;

    //! Index of time steps.
    std::vector<int> mPeriodToTimeStep;

    //! Model period to year.
    std::vector<int> mPeriodToYear;

    //! Year to model period map object.
    std::map<int,int> mYearToPeriod;
    
    //! Debugging flag to make sure the modeltime params have been
    //! set before attempting to get the modeltime attributes.
    bool mIsInitialized;

    // member functions
    void initMembers( const std::map<int, int>& aYearToTimeStep );
    
    //! Private constructor to prevent creating a Modeltime.
    Modeltime();
    //! Private undefined copy constructor to prevent creating another Modeltime.
    Modeltime( const Modeltime& aModelTime );
    //! Private undefined assign operator to prevent creating another Modeltime.
    Modeltime& operator=( const Modeltime& aModelTime );
public:
    static const Modeltime* getInstance();
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
    // IRoundTrippable methods
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    
    void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;

    static const std::string& getXMLNameStatic();
    
    int getBasePeriod() const;
    int getStartYear() const;
    int getEndYear() const;
    int gettimestep( const int aPeriod ) const { return mPeriodToTimeStep[ aPeriod ]; } // years from last to current per
    int getmaxper() const { return mMaxPeriod; }  // max modeling periods

    int getper_to_yr( const int aPeriod ) const;
    int getyr_to_per( const int aYear ) const;

    bool isModelYear( const int aYear ) const;

    int getFinalCalibrationPeriod() const;
};

#endif // _MODEL_TIME_H_
