#ifndef _TRANSECTOR_H_
#define _TRANSECTOR_H_
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
* \file tran_sector.h
* \ingroup Objects
* \brief The Transportation Sector header file. 
* \author Marshall Wise, Sonny Kim, Josh Lurz
* \date $Date: 2006/02/21 15:32:03 $
* \version $Revision: 1.13.2.1 $
*/

#include <vector>
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include "sectors/include/energy_final_demand.h"

// Forward declarations
class GDP;

/*! 
* \ingroup Objects
* \brief This class represents a specialized transporation demand sector.
* \details The Transportation sector is a specialized type of energy demand.
* \author Marshall Wise, Sonny Kim, Josh Lurz
*/
class TranSector : public EnergyFinalDemand
{
public:
    TranSector();
    
    virtual void initCalc( const std::string& aRegionName,
                           const GDP* aGDP,
                           const Demographic* aDemographics,
                           const int aPeriod );

    virtual void completeInit( const std::string& aRegionName,
                               const IInfo* aRegionInfo );

    
    virtual void setFinalDemand( const std::string& aRegionName,
                                const Demographic* aDemographics,
                                const GDP* aGDP,
                                const int aPeriod );

    static const std::string& getXMLNameStatic();
protected:
    std::vector<double> mPercentLicensed; //!< Percent of population licensed
    
    //! Cummulative technical change.
    objects::PeriodVector<double> mTechChange;

    double mBaseScaler; //!< Calculated constant scaler to scale base output
    double mBaseScalerNotLic; //!< Calculated constant scaler to scale base unlicensed output
    
    void calcBaseScalers( const GDP* aGDP, const int aPeriod );
    double calcTechChange( const int aPeriod ) const;
    
    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ); 
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
};

#endif // _TRANSSECTOR_H_

