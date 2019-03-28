#ifndef _INDIRECT_EMISSIONS_CALCULATOR_H_
#define _INDIRECT_EMISSIONS_CALCULATOR_H_
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
 * \file indirect_emissions_calculator.h
 * \ingroup Objects
 * \brief IndirectEmissionsCalculator class header file.
 * \author Josh Lurz
 */
#include "util/base/include/default_visitor.h"
#include "util/base/include/time_vector.h"
#include "util/base/include/value.h"
#include <map>
#include <string>

/*! 
 * \ingroup Objects
 * \brief Object which calculates total indirect emissions for all sectors and
 *        upstream emissions coefficients.
 * \details Calculates the indirect emissions and upstream coefficients for all
 *          sectors within a region. Indirect emissions are all upstream
 *          emissions for the fuels used by a sector. They do not include direct
 *          emissions from the sector. The upstream emissions coefficient for a
 *          sector includes upstream and direct emissions since it is used by
 *          other sectors to calculate their indirect emissions. Indirect
 *          emissions are equal to the sum of all inputs multiplied by their
 *          upstream coefficients.
 * \warning This class should only be used to generate indirect emissions and
 *          upstream coefficients for one region at a time.
 * \todo Once the Access database is removed this can be simplified to only work
 *       for one period at a time and not keep track of total sector level
 *       indirect emissions.
 * \author Josh Lurz
 */
class IndirectEmissionsCalculator : public DefaultVisitor {
public:
    IndirectEmissionsCalculator();

    double getUpstreamEmissionsCoefficient( const std::string& aSector,
                                            const int aPeriod ) const;

    double getIndirectEmissions( const std::string& aSector,
                                 const int aPeriod ) const;

    // IVisitor interface methods.
    void startVisitSector( const Sector* aSector, const int aPeriod );
    
    void endVisitSector( const Sector* aSector, const int aPeriod );
    
    void startVisitTechnology( const Technology* aTechnology,
                               const int aPeriod );
private:
    /*!
     * \brief Storage type for the coefficients that contains one value per
     *        sector and period.
     * \todo Convert to a one period storage structure once the Access database
     *       is removed.
     */
    typedef std::map<std::string, objects::PeriodVector<Value> > DoubleMap;
    
    //! Map of sector name to upstream emissions coefficient.
    DoubleMap mUpstreamEmissionsCoefficients;

    //! Map of sector name to indirect emissions.
    DoubleMap mIndirectEmissions;

    //! Current sector name.
    std::string mCurrSectorName;

    //! Total emissions from the current sector.
    double mCurrTotalEmissions;

    //! Total indirect emissions from the current sector.
    double mCurrIndirectEmissions;
    
    //! Total output for the current sector.
    double mCurrOutput;
};

#endif // _INDIRECT_EMISSIONS_CALCULATOR_H_
