#ifndef _LAND_CARBON_DENSITIES_H_
#define _LAND_CARBON_DENSITIES_H_
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
 * \file land_carbon_densities.h
 * \ingroup Objects
 * \brief The LandCarbonDensities class header file.
 * \author James Blackwood
 */

#include <xercesc/dom/DOMNode.hpp>
#include "ccarbon_model/include/asimple_carbon_calc.h"

/*!
 * \brief A simple carbon content calculator used for unmanaged land leaves.
 * \details The unmanaged land carbon calculator reads-in value for above and
 *          below ground carbon for all years. These read-in values are then
 *          used to calculate the carbon value of the land and the net land use
 *          change emissions.
 *
 *          <b>XML specification for LandCarbonDensities</b>
 *          - XML name: \c unmanaged-carbon-calc
 *          - Contained by: UnmanagedLandLeaf
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements:
 *              - \c above-ground-carbon LandCarbonDensities::mAboveGroundCarbon
 *              - \c below-ground-carbon LandCarbonDensities::mBelowGroundCarbon
 */
class LandCarbonDensities : public ASimpleCarbonCalc {
public:
    LandCarbonDensities();
    virtual ~LandCarbonDensities();
    virtual LandCarbonDensities* clone() const;

    static const std::string& getXMLNameStatic();

    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;

    virtual void completeInit();

    virtual double getActualAboveGroundCarbonDensity( const int aYear ) const;
    
    virtual void setActualAboveGroundCarbonDensity( const double aAboveGroundCarbonDensity,
                                           const int aYear );

    virtual double getActualBelowGroundCarbonDensity( const int aYear ) const;

    virtual void setActualBelowGroundCarbonDensity( const double aBelowGroundCarbonDensity,
                                           const int aYear );

    virtual void setMatureAge( const int aMatureAge );    
	
	virtual int getMatureAge( ) const;

protected:
    //! Actual above ground carbon content by year.  Note that these are for
    //! GCAM model years only since historical years have a seperate value.
    objects::YearVector<double> mAboveGroundCarbon;

    //! Actual below ground carbon content by year.  Note that these are for
    //! GCAM model years only since historical years have a seperate value.
    objects::YearVector<double> mBelowGroundCarbon;

    //! Average above ground carbon content (read in).
    double mAvgAboveGroundCarbon;

    //! Average below ground carbon content (read in).
    double mAvgBelowGroundCarbon;
	
    //! Age at maturity.  This is used to grow forests slowly.
    int mMatureAge;      
    
    //! Above ground carbon adjustments for feedbacks
    objects::YearVector<double> mAboveGroundAdj;

    //! Below ground carbon adjustments for feedbacks
    objects::YearVector<double> mBelowGroundAdj;
};

#endif // _LAND_CARBON_DENSITIES_H_
