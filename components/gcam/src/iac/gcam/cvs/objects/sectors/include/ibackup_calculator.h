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

#ifndef _IBACKUP_CALCULATOR_H_
#define _IBACKUP_CALCULATOR_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*!
 * \file ibackup_calculator.h
 * \ingroup Objects
 * \brief The IBackupCalculator interface header file.
 * \author Marshall Wise, Josh Lurz
 */

#include "util/base/include/istandard_component.h"

// Forward declaration
class IInfo;

/*!
 * \ingroup Objects
 * \brief Interface which defines methods for calculating an average and
 *        marginal amount of backup capacity required per unit of output.
 * \details Defines an interface to an object which determines the backup
 *          capacity required for a Sector. The backup capacity is determined
 *          per unit of output, but may use trial values to allow computation
 *          based on the total output. Backup requirements are specific to
 *          sectors that produce electricity.
 * \author Josh Lurz
 */
class IBackupCalculator : public IParsedComponent {
public:
    // Clone operator must be declared explicitly even though it is inherited
    // from IStandardComponent so that the return type can be changed. Since
    // this class is a subtype of IStandardComponent, this is legal and referred
    // to as a covariant return type.
    virtual IBackupCalculator* clone() const = 0;

    /*!
     * \brief Pass parameter information into the backup calculator from the technology.
     * \param aTechInfo An info object containing information to be passed to the backup component.
     * \author Steve Smith
     */
    virtual void initCalc( const IInfo* aTechInfo ) = 0;

    /*!
     * \brief Compute backup required for the marginal unit of energy output.
     * \details Compute backup required per resource energy output on the margin
     *          (since energy output is what the modeled market is based on).
     * \param aSector The name of the sector which requires backup capacity.
     * \param aElectricSector The name of the electricity sector into which the
     *        sector having a backup amount calculated for will feed.
     * \param aResource The name of the resource the sector consumes.
     * \param aRegion Name of the containing region.
     * \param aReserveMargin Reserve margin for the electricity sector.
     * \param aAverageGridCapacityFactor The average electricity grid capacity
     *        factor.
     * \param aPeriod Model period.
     * \return Reserve capacity per marginal intermittent electricity resource
     *         output.
     */
    virtual double getMarginalBackupCapacity( const std::string& aSector,
                                              const std::string& aElectricSector,
                                              const std::string& aResource,
                                              const std::string& aRegion,
                                              const double aReserveMargin,
                                              const double aAverageGridCapacityFactor,
                                              const int aPeriod ) const = 0;

    /*!
     * \brief Compute the average backup required per unit for the intermittent
     *        subsector.
     * \details Computes the average quantity of backup capacity required per
     *          unit of energy output.
     * \param aSector The name of the sector which requires backup capacity.
     * \param aElectricSector The name of the electricity sector into which the
     *        sector having a backup amount calculated for will feed.
     * \param aResource The name of the resource the sector consumes.
     * \param aRegion Name of the containing region.
     * \param aReserveMargin Reserve margin for the electricity sector.
     * \param aAverageGridCapacityFactor The average electricity grid capacity
     *        factor.
     * \param aPeriod Model period.
     * \return The average backup capacity required per unit of output.
     */
    virtual double getAverageBackupCapacity( const std::string& aSector,
                                             const std::string& aElectricSector,
                                             const std::string& aResource,
                                             const std::string& aRegion,
                                             const double aReserveMargin,
                                             const double aAverageGridCapacityFactor,
                                             const int aPeriod ) const = 0;
};

#endif // _IBACKUP_CALCULATOR_H_
