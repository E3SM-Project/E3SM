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
 * \file calibrate_resource_visitor.cpp
 * \ingroup Objects
 * \brief The CalibrateResourceVisitor class source file.
 * \author Kate Calvin
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "util/base/include/calibrate_resource_visitor.h"
#include "resources/include/subresource.h"
#include "resources/include/renewable_subresource.h"
#include "resources/include/resource.h"
#include "resources/include/grade.h"

#include "marketplace/include/marketplace.h"
#include "containers/include/scenario.h"

using namespace std;

extern Scenario* scenario;

/*!
 * \brief Constructor
 * \param aRegionName Name of the region if starting the visiting below the
 *        region level.
 */
CalibrateResourceVisitor::CalibrateResourceVisitor( const string& aRegionName )
:mCurrentRegionName( aRegionName )
{
}

void CalibrateResourceVisitor::startVisitResource( const AResource* aResource,
                                                   const int aPeriod )
{
    mCurrentResourceName = aResource->getName();
}

void CalibrateResourceVisitor::endVisitResource( const AResource* aResource,
                                                 const int aPeriod )
{
    mCurrentResourceName.clear();
}


void CalibrateResourceVisitor::startVisitSubResource( const SubResource* aSubResource, const int aPeriod ) {
    // If calibration is active and a calibrated production quantity was read in, 
    // then we need to calculate the price adder needed to produce that quantity
    if( aSubResource->mCalProduction[ aPeriod ] != -1.0 && aPeriod > 0 ){
        // First, calculate the cumulative production.  This is equal to 
        // cumulative production in the previous period plus the calibrated
        // production times the timestep.  Note: this assumes constant production
        // in all years with in a timestep. This is necessary to prevent erratic
        // production in the future.

        double tempCumulProd = aSubResource->mCumulProd[ aPeriod - 1 ] + 
            aSubResource->mCalProduction[ aPeriod ] * scenario->getModeltime()->gettimestep( aPeriod );

        // Next, determine which grade of resource is produced to get to the
        // cumulative production needed.
        double temp_cumulative = 0.0;
        int gr_ind = 0;
        double gr_avail = 0.0;
        while ( temp_cumulative < tempCumulProd ) {
            gr_avail = aSubResource->mGrade[ gr_ind ]->getAvail();
            temp_cumulative += gr_avail;
            gr_ind++;
        }

        // Then, calculate the fraction of the next grade that will be produced
        double fractGrade = 0.0;
        if ( gr_avail > 0.0 ) {
            fractGrade = ( tempCumulProd - ( temp_cumulative - gr_avail ) )
                / gr_avail;
        }

        // Next, calculate the effective price.  This is the price needed
        // to produce the calibrated production quantity in this period.
        // mEffectivePrice = cost of the next highest grade - 
        // ( 1 - fractGrade )*( cost of higher grade - cost of lower grade )
        double low_cost = 0.0;
        if ( gr_ind > 0 ) {
            low_cost = aSubResource->mGrade[ gr_ind - 1 ]->getCost( aPeriod );
        }
        double tempEffectivePrice = aSubResource->mGrade[ gr_ind ]->getCost( aPeriod ) - 
            ( 1 - fractGrade ) * ( aSubResource->mGrade[ gr_ind ]->getCost( aPeriod ) - low_cost );

        double mktPrice = scenario->getMarketplace()->getPrice( mCurrentResourceName, 
                                                                mCurrentRegionName, 
                                                                aPeriod );

        // Finally, calculate the price adder. This is the difference between the
        // effective price and the global price
        const_cast<SubResource*>( aSubResource )->mPriceAdder[ aPeriod ] = tempEffectivePrice - mktPrice;
    }
}

void CalibrateResourceVisitor::startVisitSubRenewableResource( const SubRenewableResource* aSubResource, 
                                                                                const int aPeriod ) {
    // If calibration is active and a calibrated production quantity was read in, 
    // then we need to calculate the price adder needed to produce that quantity
    if( aSubResource->mCalProduction[ aPeriod ] != -1.0 && aPeriod > 0 ){
        // First, calculate the annual production
        double tempAnnualProd = aSubResource->mCalProduction[ aPeriod ];

        // Next, determine which grade of resource is produced to get to the
        // annual production needed.
        double temp_annual = 0.0;
        int gr_ind = 0;
        double gr_avail = 0.0;
        double prev_avail = 0.0;
        while ( temp_annual < tempAnnualProd ) {
            // For renewable resources, getAvail() returns the fraction of the maxSubResource
            // that can be produced.  So, total production is getAvail() * maxSubResource
            gr_avail = ( aSubResource->mGrade[ gr_ind ]->getAvail() - prev_avail ) 
                                            * aSubResource->getMaxSubResource();
            prev_avail = aSubResource->mGrade[ gr_ind ]->getAvail();
            temp_annual += gr_avail;
            gr_ind++;
        }

        // Then, calculate the fraction of the next grade that will be produced
        double fractGrade = 0.0;
        if ( gr_avail > 0.0 ) {
            fractGrade = ( tempAnnualProd - ( temp_annual - gr_avail ) )
                / gr_avail;
        }

        // Next, calculate the effective price.  This is the price needed
        // to produce the calibrated production quantity in this period.
        // mEffectivePrice = cost of the next highest grade - 
        // ( 1 - fractGrade )*( cost of higher grade - cost of lower grade )
        double low_cost = 0.0;
        if ( gr_ind > 0 ) {
            low_cost = aSubResource->mGrade[ gr_ind - 1 ]->getCost( aPeriod );
        }
        
        double tempEffectivePrice = aSubResource->mGrade[ gr_ind ]->getCost( aPeriod ) - 
            ( 1 - fractGrade ) * ( aSubResource->mGrade[ gr_ind ]->getCost( aPeriod ) - low_cost );


        double mktPrice = scenario->getMarketplace()->getPrice( mCurrentResourceName, 
                                                                mCurrentRegionName, 
                                                                aPeriod );

        // Finally, calculate the price adder. This is the difference between the
        // effective price and the global price
        const_cast<SubRenewableResource*>( aSubResource )->mPriceAdder[ aPeriod ] = tempEffectivePrice - mktPrice;
    }
}
