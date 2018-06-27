#ifndef _ITARGET_H_
#define _ITARGET_H_
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
 * \file itarget.h
 * \ingroup Objects
 * \brief The ITarget interface file.
 * \author Josh Lurz
 */

/*!
 * \brief Interface to represent a target.
 */
class ITarget {
public:
    /*!
     * \brief A year flag to indicate that getStatus should use getYearOfMaxTargetValue
     *        to figure out the appropriate year to get the status in.
     */
    static int getUseMaxTargetYearFlag() {
        static const int MAX_TARGET_YEAR = -1;
        return MAX_TARGET_YEAR;
    }

    /*!
     * \brief Get the status of the last trial for a given period.
     * \details Returns whether the last trial was over the target, under the
     *          target, or solved the target for a given period.
     * \param aTolerance Solution tolerance.
     * \param aYear Year in which to check the target.
     * \return The status of the last trial.
     */
    virtual double getStatus( const int aYear ) const = 0;
    
    /*!
     * \brief Find the year in which the target is at it's maximum value over all
     *        valid model years.
     * \details Should getStatus be called with year equal to getUseMaxTargetYearFlag()
     *          it will use this method to find the appropriate year to check.
     * \return The year in which the target is at it's maximum.
     */
    virtual int getYearOfMaxTargetValue() const = 0;
};

#endif // _ITARGET_H_
