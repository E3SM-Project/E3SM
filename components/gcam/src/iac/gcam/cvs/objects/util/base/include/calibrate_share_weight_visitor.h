#ifndef _CALIBRATE_SHARE_WEIGHT_VISITOR_H_
#define _CALIBRATE_SHARE_WEIGHT_VISITOR_H_
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
 * \file calibrate_share_weight_visitor.h
 * \ingroup Objects
 * \brief CalibrateShareWeightVisitor class header file.
 * \author Pralit Patel
 */

#include "util/base/include/default_visitor.h"
#include <string>

class GDP;

/*! 
 * \ingroup Objects
 * \brief A visitor which determines if a logit nest has calibrated values
 *        and will adjust share weights algebraiclly to reproduce those values
 *        at the current set of prices.
 * \details For both subectors and technologies we first sum all of the calibration
 *          values and find the largest child in terms of largest calibration value
 *          so that we may make the share weights relative to that subsector/technology.
 *          We then back out the appropriate share weight to reproduce the calibration shares
 *          using the following equation: (TODO: figure out how to format this)
 *          
 *          s_i = ( cal_i / cal_total ) / ( scaled_gdp_percapita ^ fuel_pref_elasticity_i )
 *          sw_i = ( s_i / s_r ) * ( price_r / price_i ) ^ logit_exp
 *
 *          Where _i indicates the current child, and _r indicates the child we are
 *          making the share weights relative to.
 * \author Pralit Patel
 * \note Calibration of detailed buildings are also handled by this class.
 * \warning This class never actually checks whether calibration is active.
 * \warning This methodology will not work if in a given nest there is a mix
 *          of subsectors/technologies with and without calibrated values.
 */
class CalibrateShareWeightVisitor : public DefaultVisitor {
public:

    CalibrateShareWeightVisitor( const std::string& aRegionName, const GDP* aGDP );

    // Documentation for visitor methods is inherited.
    virtual void startVisitSector( const Sector* aSector,
                                   const int aPeriod );

    virtual void endVisitSector( const Sector* aSector,
                                const int aPeriod );

    virtual void startVisitSubsector( const Subsector* aSubsector,
                                     const int aPeriod );

    // need to visit BuildingGenericDmdTechnology to handle the specialized
    // calibration for detailed buildings
    virtual void startVisitBuildingGenericDmdTechnology( const BuildingGenericDmdTechnology* aBuildingTech,
                                                        const int aPeriod );

private:
    //! Name of the Region the for which we are calibrating
    std::string mCurrentRegionName;

    //! Name of the sector currently being tabulated.
    std::string mCurrentSectorName;

    //! The GDP for the current region. TODO: is this really necessary
    const GDP* mGDP;
};

#endif // _CALIBRATE_SHARE_WEIGHT_VISITOR_H_
