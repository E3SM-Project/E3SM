#ifndef _IINVESTOR_H_
#define _IINVESTOR_H_
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
 * \file iinvestor.h
 * \ingroup Objects
 * \brief The IInvestor interface file.
 * \author Josh Lurz
 */
#include <string>

class Tabs;
class Demographic;
class IInvestable;
class NationalAccount;

/*! 
 * \ingroup Objects
 * \brief An interface to a class which controls SGM sector level investment in
 *        new technologies.
 * \details The IInvestor interface represents a sector level decision maker
 *          which controls calculating the total investment for a
 *          ProductionSector, and distributing it to the subsectors. The
 *          IInvestor will make this decision once per iteration. Each
 *          ProductionSector has it's own independent IInvestor, although the
 *          actions of the IInvestors are linked through the investment market.
 * \author Josh Lurz
 */
class IInvestor
{
public:
    IInvestor();
    virtual ~IInvestor();

    /*!
     * \brief Complete the initialization of an IInvestor.
     * \details Finishes initializes an IInvestor before it is used.
     * \param aRegionName Region containing the investor.
     * \param aSectorName Sector for which the IInvestor is determining
     *        investment.
     */
    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName ) = 0;

    /*!
     * \brief Initialization of an IInvestor for each period.
     * \param aRegionName Region containing the investor.
     * \param aSectorName Sector for which the IInvestor is determining
     *        investment.
     */
    virtual void initCalc( std::vector<IInvestable*>& aInvestables,
                   NationalAccount& aNationalAccount, 
                   const Demographic* aDemographic,
                   const int aPeriod ) = 0;

    // TODO: Inherit and make documentation inherited.
    virtual void XMLParse( const xercesc::DOMNode* node ) = 0; 
    virtual void toDebugXML( const int period, std::ostream& out,
                             Tabs* tabs ) const = 0;
    virtual void toInputXML( std::ostream& out, Tabs* tabs ) const = 0;
    
    /*!
     * \brief Calculate a total investment level for the sector and distribute
     *        it to the subsectors of the sector.
     * \details This is the main investment method. It is called by the
     *          ProductionSector to determine the total quantity of investment
     *          for the sector and to distribute the investment to the
     *          subsectors of the ProductionSector.
     * \param aInvestables The subsectors of the ProductionSector as an
     *        IInvestable vector.
     * \param aNationalAccount Regional national accounts container.
     * \param aDemographic Regional demographics container.
     * \param aPeriod Model period for which to calculate investment.
     */
    virtual double calcAndDistributeInvestment( std::vector<IInvestable*>& aInvestables,
                                                NationalAccount& aNationalAccount, 
                                                const Demographic* aDemographic,
                                                const int aPeriod ) = 0;
    virtual void setEfficiencyConditions( std::vector<IInvestable*>& aInvestables,
                                          NationalAccount& aNationalAccount, 
                                          const Demographic* aDemographic,
                                          const int aPeriod ) const = 0;
};

//! Constructor
inline IInvestor::IInvestor(){
}

//! Destructor
inline IInvestor::~IInvestor(){
}

#endif // _IINVESTOR_H_
