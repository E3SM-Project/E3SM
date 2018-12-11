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
 * \file carbon_model_utils.cpp
 * \ingroup Objects
 * \brief CarbonModelUtils class source file.
 * \author Jim Naslund and Ming Chang
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include <cfloat>

#include "ccarbon_model/include/carbon_model_utils.h"
#include "util/base/include/util.h"
#include "land_allocator/include/land_use_history.h"
#include "climate/include/iclimate_model.h"

using namespace std;

/*!
 * \brief Constructor.
 */
CarbonModelUtils::CarbonModelUtils(){
}

//! Default destructor
CarbonModelUtils::~CarbonModelUtils() {
}

/*!
 * \brief Get the land usage for a year.
 * \details Returns the land usage for a year. An appropriate historical or
 *          calculated value will be returned depending on the year.
 *          Land-use history values always have precedence
 * \param aYear the year.
 * \param aLandUseHistory the land use history object.
 * \param aHistoricalShare the historical share.
 * \param aLandUse the land use object.
 * \return Land usage for the year.
 * \warning This function does not check if the model projection data is valid.
 *          For example, its potentially possible that land-use data for 2020
 *          is requested before that period is calculated.
 * \todo Address the above warning.
 */
double CarbonModelUtils::getLandUse( const unsigned int aYear,
                                     const LandUseHistory* aLandUseHistory,
                                     const objects::PeriodVector<double>& aLandUse ){
    // If the year is within the range of the history use the historical
    // allocation. The land use history may be null if none was read-in.
    unsigned int maxHistoryYear = 0;
    if( aLandUseHistory ){
        maxHistoryYear = aLandUseHistory->getMaxYear();
    }

    unsigned int basePeriod =
        max( scenario->getModeltime()->getyr_to_per( max( static_cast<unsigned int>( 1975 ), maxHistoryYear ) ), 1 );
    basePeriod = min( basePeriod, static_cast<unsigned int>(scenario->getModeltime()->getmaxper() - 1) ); 
    // Store the first calculated year to save time.
    const unsigned int baseYear
        = static_cast<unsigned int>( scenario->getModeltime()->getper_to_yr( basePeriod ) );

    double landUse;

    if( aYear <= maxHistoryYear ){
        landUse = aLandUseHistory->getAllocation( aYear );
    }
    // If the year is between the last historical year and the first
    // calculated year interpolate between the two.
    else if( aYear <= baseYear && maxHistoryYear != 0 ){
        landUse = util::linearInterpolateY( aYear, maxHistoryYear, baseYear,
                                            aLandUseHistory->getAllocation( maxHistoryYear ),
                                            aLandUse[ basePeriod ] );
    }
    // Otherwise use data interpolated from the current model projection.
    else {
        landUse = interpYearHelper( aLandUse, aYear );
    }
    return landUse;
}

/*!
 * \brief A static function to return the starting year to index the arrays.
 * \todo Use a read-in value for the start year.
 * \return The start year.
 */
int CarbonModelUtils::getStartYear(){
    const static int START_YEAR = scenario->getClimateModel()
        ? scenario->getClimateModel()->getCarbonModelStartYear() : scenario->getModeltime()->getStartYear();
    return START_YEAR;
}

/*!
 * \brief Return the last year of the climate calculation.
 * \author Jim Naslund
 * \return The last year of the climate calculation.
 */
int CarbonModelUtils::getEndYear(){
    const static int END_YEAR = scenario->getModeltime()->getEndYear();
    return END_YEAR;
}

/*
 * \brief Returns a parameter which defines the time scale for the soil
 *        emissions decay function.
 * \return Soil decay function time scale parameter.
 * \todo This should be dynamic by land type.
 */
double CarbonModelUtils::getSoilTimeScale(){
    const static double SOIL_TIME_SCALE = 40;
    return SOIL_TIME_SCALE;
}

/*!
 * \brief Helper function to interpolate a value for a year from a PeriodVector.
 * \details Calculates a linearly interpolated value for the year. If the year
 *          is before the first period of the vector, the first value is
 *          returned. If the year is after the last period of the vector, the
 *          last value is used. Otherwise a value is linearly interpolated
 *          between the nearest two periods.
 * \param aPeriodVector Vector from which to interpolate the value.
 * \param aYear Year for which to interpolate a value.
 * \return Interpolated value for the year.
 */
double CarbonModelUtils::interpYearHelper( const objects::PeriodVector<double>& aPeriodVector,
                                           const unsigned int aYear ){
    // If the year is before the first period of the model use the value
    // in the base period.
    const Modeltime* modeltime = scenario->getModeltime();
    if( aYear <= static_cast<unsigned int>( modeltime->getStartYear() ) ){
        return *aPeriodVector.begin();
    }

    // If the year is after the end period use the value in the last period.
    if( aYear > static_cast<unsigned int>( modeltime->getEndYear() ) ){
        return *aPeriodVector.last();
    }
    
    // Find the period containing aYear. This cannot be zero because the year
    // was already checked against the start year.
    int currPeriod = modeltime->getyr_to_per( aYear );

    // Find the last year of the current period.
    int lastYear = modeltime->getper_to_yr( currPeriod );

    // Find the first year of the current period.
    int firstYear = modeltime->getper_to_yr( currPeriod - 1 );

    // Interpolate the result.
    return util::linearInterpolateY( aYear, firstYear, lastYear,
                                     aPeriodVector[ currPeriod - 1 ],
                                     aPeriodVector[ currPeriod ] );
}

/*!
 * \brief Helper function to get a value for a year from a YearVector.
 * \details Determines a value for a given year from a YearVector. If the year
 *          is within the range of the year vector the value for the year will
 *          be returned, otherwise a value will be extrapolated. The
 *          extrapolation considers all values before the first year equal to
 *          the first year, and all values after the end year to be equal to the
 *          end year.
 * \param aYearVector Vector from which to interpolate the value.
 * \param aStartYear First year of the climate model.
 * \param aEndYear Last year of the climate model.
 * \param aYear Year for which to interpolate a value.
 * \return Interpolated value for the year.
 */
double CarbonModelUtils::interpYearHelper( const objects::YearVector<double>& aYearVector,
                                           const unsigned int aStartYear,
                                           const unsigned int aEndYear,
                                           const unsigned int aYear ){
    // If the year is before the first year of the carbon cycle use the value in the
    // first year.
    if( aYear < aStartYear ){
        return *aYearVector.begin();
    }

    // If the year is after the last year of the carbon cycle use the value in the last year.
    if( aYear > aEndYear ){
        return *aYearVector.last();
    }

    // Return the value from inside the range of the vector.
    return aYearVector[ aYear ];
}
