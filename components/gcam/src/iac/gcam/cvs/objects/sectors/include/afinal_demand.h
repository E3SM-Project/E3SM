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

#ifndef _AFINAL_DEMAND_H_
#define _AFINAL_DEMAND_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file afinal_demand.h
 * \ingroup Objects
 * \brief The AFinalDemand abstract base class header file.
 * \author Josh Lurz
 */
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/ivisitable.h"
#include "util/base/include/iround_trippable.h"
#include "util/base/include/iparsable.h"
// Forward declarations
class GDP;
class Demographic;
class IInfo;
class Tabs;

/*! 
 * \ingroup Objects
 * \brief A description of a final demand, the end users in the economy.
 * \details Final demand do not produce outputs, but only consume goods which
 *          can include services. They are not linked to each other directly,
 *          but are in aggregate the consumers of the economy.
 */

class AFinalDemand: public IParsable,
                    public IRoundTrippable,
                    public IVisitable
{
public:
    /*!
     * \brief Destructor.
     */
    virtual ~AFinalDemand();

    // Documentation is inherited
    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const = 0;
    
    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const = 0;
    
    /*!
     * \brief Get the name of the final demand.
     * \return The name of the final demand.
     */
    virtual const std::string& getName() const = 0;
    
    /*!
     * \brief Complete the initialization of the final demand.
     * \details This method is called after all data is parsed so that the
     *          object may perform any final initializations.
     * \param aRegionName Region name.
     * \param aRegionInfo Regional information container.
     */
    virtual void completeInit( const std::string& aRegionName,
                               const IInfo* aRegionInfo ) = 0;

    /*!
     * \brief Initialize the final demand for a given period.
     * \details This method is called at the beginning of each period so that
     *          the object may perform any initializations for the model period.
     * \param aRegionName Region name.
     * \param aGDP Regional GDP.
     * \param aDemograhics Region demographics.
     * \param aPeriod Model period.
     */
    virtual void initCalc( const std::string& aRegionName,
                           const GDP* aGDP,
                           const Demographic* aDemographics,
                           const int aPeriod ) = 0;
    
    /*!
     * \brief Calculate the quantity of the demand and add it to the
     *        marketplace.
     * \details Calculates the final demand for the input good based on a set of
     *          drivers, such as GDP, population and price. The demand for the
     *          final good is then added to the marketplace.
     * \param aRegionName Region name.
     * \param aDemographics Regional demographics.
     * \param aGDP Regional GDP container.
     * \param aPeriod Model period.
     */
    virtual void setFinalDemand( const std::string& aRegionName,
                                 const Demographic* aDemographics,
                                 const GDP* aGDP,
                                 const int aPeriod ) = 0;

    /*!
     * \brief Get the market price of the service weighted by the quantity of
     *          service demanded.
     * \todo Should this be weighted on final energy? Service is arbitrary.
     * \param aRegionName Region name.
     * \param aPeriod Model period.
     */
    virtual double getWeightedEnergyPrice( const std::string& aRegionName,
                                           const int aPeriod ) const = 0;

    /*!
     * \brief Tabulate all known fixed demands.
     * \details Adds to a central location the values of any known fixed
     *          demands.
     * \param aRegionName Region name.
     * \param aPeriod Model period.
     */
    virtual void tabulateFixedDemands( const std::string& aRegionName,
                                       const Demographic* aDemographic,
                                       const GDP* aGDP,
                                       const int aPeriod ) = 0;

    /*!
     * \brief Write output to the CSV file.
     * \param aRegionName Region name.
     */
    virtual void csvOutputFile( const std::string& aRegionName ) const = 0;

    /*!
     * \brief Write output to the database.
     * \param aRegionName Region name.
     */
    virtual void dbOutput( const std::string& aRegionName ) const = 0;

    // Documentation is inherited.
    virtual void accept( IVisitor* aVisitor,
                         const int aPeriod ) const = 0;
};

// Inline function definitions.
inline AFinalDemand::~AFinalDemand(){}

#endif // _AFINAL_DEMAND_H_

