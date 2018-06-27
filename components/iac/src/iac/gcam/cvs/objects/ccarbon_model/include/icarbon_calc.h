#ifndef _ICARBON_CALC_H_
#define _ICARBON_CALC_H_
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
 * \file icarbon_calc.h
 * \ingroup Objects
 * \brief The ICarbonCalc interface file.
 * \author James Blackwood, Josh Lurz
 */

#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/ivisitable.h"
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"
#include "ccarbon_model/include/carbon_model_utils.h"
// Forward declarations
class IInfo;
class Tabs;
class LandUseHistory;

/*!
 * \brief An interface to an object responsible for determining the annual
 *        carbon content of a land area.
 * \details Carbon content calculators are used to determine how much carbon is
 *          contained in a land leaf over time, and any net carbon emissions
 *          This calculation can be simple (net carbon accounting) or complex 
 *          (full carbon cycle with feedbacks). The carbon content may change as 
 *          a result of internal processes (re-growth or feedbacks), as a result
 *          of land-use changes, or due changes in the technologies used to grow
 *          and harvest crops. Changes in the carbon content of the land will
 *          result in net carbon emissions. The calculated carbon content of the
 *          land is also used to determine the carbon subsidy. Carbon content
 *          calculators do not calculate land use change, this is given to them
 *          by the land leaf.
 */
class ICarbonCalc: public IVisitable,
                   public IRoundTrippable,
                   public IParsable 
{
public:
    //! Constructor
    inline ICarbonCalc();

    //! Destructor.
    inline virtual ~ICarbonCalc();

    virtual ICarbonCalc* clone() const = 0;

    // Documentation is inherited.
    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;

    // Documentation is inherited.
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const = 0;
    
    // Documentation is inherited.
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const = 0;
    
    /*!
     * \brief Complete the initialization of the carbon calculator.
     */
    virtual void completeInit() = 0;

    /*!
     * \brief Initialize the historical land use for the carbon calculation.
     * \details Initializes the carbon calculation with values for the historical
     *          land use of the leaf containing the carbon calculator.
     * \param aHistory Historical land use object which may be null.
     * \param aShare Estimated share of total land used by the containing leaf.
     */
    virtual void initLandUseHistory( const LandUseHistory* aHistory ) = 0;

    /*!
     * \brief Conduct carbon calculations for a period.
     * \details Conduct carbon calculations for a period once all information is
     *          set. This should be called once per iteration, but there may be
     *          multiple iterations per period. This may perform historical
     *          calculations up to the current period if the historical years
     *          have not been calculated yet.  Future emissions will only be
     *          calculated to the given end year.
     * \param aPeriod The current model period that is being calculated.
     * \param aEndYear The year to calculate future emissions to.
     */
    virtual void calc( const int aPeriod, const int aEndYear ) = 0;

    /*!
     * \brief Get the net land use change emissions for a given year.
     * \details Returns the net land use change emission at the year. This is an
     *          annual emission. The value is a sum of the above and below
     *          ground uptake and emissions. A positive return value represents
     *          an emission, a negative represents net uptake.
     * \param aYear Year for which to get the land use change emission.
     * \return Annual net land use change emission for the year.
     */
    virtual double getNetLandUseChangeEmission( const int aYear ) const = 0;

    /*!
     * \brief Set the total land use for a period.
     * \details Sets the total land area used by the containing land leaf for a
     *          period.
     * \param aLandUse Amount of land used.
     * \param aPeriod Model period.
     */
    virtual void setTotalLandUse( const double aLandUse,
                                  const int aPeriod ) = 0;

    virtual double getActualAboveGroundCarbonDensity( const int aYear ) const = 0;
    
    virtual void setActualAboveGroundCarbonDensity( const double aAboveGroundCarbonDensity,
                                           const int aPeriod ) = 0;

    virtual double getActualBelowGroundCarbonDensity( const int aYear ) const = 0;

    virtual void setActualBelowGroundCarbonDensity( const double aBelowGroundCarbonDensity,
                                           const int aPeriod ) = 0;
    
		virtual int getMatureAge( ) const = 0;

    virtual double getAboveGroundCarbonSubsidyDiscountFactor( ) = 0;

    virtual double getBelowGroundCarbonSubsidyDiscountFactor( ) = 0;

	virtual double getAboveGroundCarbonStock( const int aYear ) const = 0;
	
    virtual double getBelowGroundCarbonStock( const int aYear ) const = 0;
	
    virtual void setSoilTimeScale( const int aTimeScale ) = 0;
    /*!
     * \brief Set the mature age of this land cover type; used by simple carbon calculator.
     * \details Sets the mature age.  For crops this is 1 (the default), but forests
     *            grow slowly, and their emission curve is calculated using this.
     * \param aMatureAge The age (in years) this land cover type reaches full biomass.
     * \param aPeriod Model period.
     */
    virtual void setMatureAge( const int aMatureAge ) = 0;    

    // Documentation is inherited.
    virtual void accept( IVisitor* aVisitor,
                         const int aPeriod ) const = 0;

    virtual void acceptDerived( IVisitor* aVisitor,
                         const int aPeriod ) const = 0;
};

// Inline function definitions.
ICarbonCalc::ICarbonCalc(){
}

ICarbonCalc::~ICarbonCalc(){
}

#endif // _ICARBON_CALC_H_
