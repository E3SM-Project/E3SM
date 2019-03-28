#ifndef _ASIMPLE_CARBON_CALC_H_
#define _ASIMPLE_CARBON_CALC_H_
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
 * \file asimple_carbon_calc.h
 * \ingroup Objects
 * \brief The ASimpleCarbonCalc class header file.
 * \author James Blackwood
 */
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/time_vector.h"
#include "ccarbon_model/include/icarbon_calc.h"

class LandUseHistory;

/*!
 * \brief The abstract base class of all simple carbon content calculators.
 * \details This class serves to share code between simple carbon content
 *          calculators, and to differentiate them from more computationally
 *          intensive carbon content calculators.
 */
class ASimpleCarbonCalc : public ICarbonCalc {
public:
    ASimpleCarbonCalc();
    virtual ~ASimpleCarbonCalc();

    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const = 0;
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const = 0;

    virtual void completeInit() = 0;

    virtual void initLandUseHistory( const LandUseHistory* aHistory );

    virtual void calc( const int aPeriod, const int aEndYear );

    virtual double getNetLandUseChangeEmission( const int aYear ) const;

    virtual void setTotalLandUse( const double aLandUse,
                                  const int aPeriod );

    virtual double getActualAboveGroundCarbonDensity( const int aYear ) const = 0;
    
    virtual void setActualAboveGroundCarbonDensity( const double aAboveGroundCarbonDensity,
                                           const int aPeriod ) = 0;

    virtual double getActualBelowGroundCarbonDensity( const int aYear ) const = 0;

    virtual void setActualBelowGroundCarbonDensity( const double aBelowGroundCarbonDensity,
                                           const int aPeriod ) = 0;
	
	virtual int getMatureAge( ) const = 0;
    
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;

    virtual void acceptDerived( IVisitor* aVisitor, const int aPeriod ) const;

    virtual double getAboveGroundCarbonSubsidyDiscountFactor( );

    virtual double getBelowGroundCarbonSubsidyDiscountFactor( );
	
	virtual double getAboveGroundCarbonStock( const int aYear ) const;
	
    virtual double getBelowGroundCarbonStock( const int aYear ) const;

    virtual void setSoilTimeScale( const int aTimeScale );

protected:
    //! Total land used by period.
    objects::PeriodVector<double> mLandUse;

    //! Stored emissions which are necessary to clear the total emissions when
    //! recalculating a period.
    objects::PeriodVector<objects::YearVector<double>*> mStoredEmissions;

    //! Total emissions by year.
    objects::YearVector<double> mTotalEmissions;

    //! Time scale for soil carbon emissions
    int mSoilTimeScale;

    /*! 
     * \brief The land use history for the land leaf or it's parent land node.
     * \details Weak pointer to the land use history either for this leaf
     *          or the parent land type. The historical land share will be set to
     *          1 if this is the land use history for the leaf, and the historical share
     *          if it is for the parent land type.
     */
    const LandUseHistory* mLandUseHistory;
    
    //! The difference in the sigmoid curve by year offset + 1 - year offset.
    //! This value get precomputed during initcalc to avoid doing the computationally
    //! expensive operations during calc.
    std::vector<double> precalc_sigmoid_diff;
    
    //! Flag to ensure historical emissions are only calculated a single time
    //! since they can not be reset.
    bool mHasCalculatedHistoricEmiss;

    void calcAboveGroundCarbonEmission( const double aCarbonDiff,
                                        const int aYear,
                                        const int aEndYear,
                                        objects::YearVector<double>& aEmissVector);

    void calcBelowGroundCarbonEmission( const double aCarbonDiff,
                                        const int aYear,
                                        const int aEndYear,
                                        objects::YearVector<double>& aEmissVector);
private:
    void calcSigmoidCurve( const double aCarbonDiff,
                           const int aYear,
                           const int aEndYear,
                           objects::YearVector<double>& aEmissVector);
};

#endif // _ASIMPLE_CARBON_CALC_H_
