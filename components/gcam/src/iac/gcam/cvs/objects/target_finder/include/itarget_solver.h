#ifndef _ITARGET_SOLVER_H_
#define _ITARGET_SOLVER_H_
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
 * \file itarget_solver.h
 * \ingroup Objects
 * \brief The ITargetSolver interface file.
 * \author Pralit Patel
 */

/*!
 * \brief Interface to represent a target finder solution algorithm.
 */
class ITargetSolver {
public:
    /*!
     * \brief A flag to indicate an undefined parameter.
     */
    static int undefined() {
        static const int UNDEFINED = -1;
        return UNDEFINED;
    }
    
    /*!
     * \brief Get the next trial price to use.
     * \details The first value of the return is the next trial price to use,
     *          the second is a boolean to indicate if we have already solved
     *          the target.
     * \return A pair of values including the next trial value to use and if we
     *         have arrived at the target.
     */
    virtual std::pair<double, bool> getNextValue() = 0;
    
    /*!
     * \brief Get the number of attempts this algorithm has taken at finding a
     *        solution.
     * \return The number of iterations by this solver.
     */
    virtual unsigned int getIterations() const = 0;
};

#endif // _ITARGET_SOLVER_H_
