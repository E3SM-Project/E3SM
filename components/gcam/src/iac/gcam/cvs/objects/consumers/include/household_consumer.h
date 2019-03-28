#ifndef _HOUSEHOLD_CONSUMER_H_
#define _HOUSEHOLD_CONSUMER_H_
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
* \file household_consumer.h
* \ingroup Objects
* \brief HouseholdConsumer class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <xercesc/dom/DOMNode.hpp>

#include "consumers/include/consumer.h"

class NationalAccount;
class Demographic;
class Tabs;
class MoreSectorInfo;
class IVisitor;

/*! 
* \ingroup Objects
* \brief An object representing a set of households.
* \details TODO
* \author Pralit Patel, Sonny Kim
*/

class HouseholdConsumer : public Consumer
{
    friend class SocialAccountingMatrix;
    friend class DemandComponentsTable;
    friend class SectorReport;
    friend class SGMGenTable;
    friend class GovtResults;
    friend class XMLDBOutputter;
public:
    HouseholdConsumer();
    virtual HouseholdConsumer* clone() const;

    void copyParam( const BaseTechnology* baseTech,
                    const int aPeriod );

    void copyParamsInto( HouseholdConsumer& householdConsumerIn,
                         const int aPeriod ) const;

    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName );
    
    virtual void initCalc( const MoreSectorInfo* aMoreSectorInfo,
                           const std::string& aRegionName, 
                           const std::string& aSectorName,
                           NationalAccount& nationalAccount, 
                           const Demographic* aDemographics,
                           const double aCapitalStock,
                           const int aPeriod );

    virtual void operate( NationalAccount& aNationalAccount, const Demographic* aDemographics, 
        const MoreSectorInfo* aMoreSectorInfo, const std::string& aRegionName, 
        const std::string& aSectorName, const bool aIsNewVintageMode, const int aPeriod );

    static const std::string& getXMLNameStatic();

    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;
    double getLaborSupply() const;
protected:
    bool isCoefBased() const { return true; }
    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
private:
    void calcBaseScalerLand( const std::string& regionName, const int period );
    void calcBaseScalerLaborMale( const std::string& regionName, const int&period );
    void calcBaseScalerLaborFemale( const std::string& regionName, const int period );
    void calcBaseScalerSavings( const std::string& regionName, const int period );
    void calcBaseLaborDemandPerHH( const std::string& regionName, const int period );
    void calcBaseLandDemandPerHH( const std::string& regionName, const int period );
    void calcFactorDemand( const std::string& regionName, int period );
    void calcFactorSupply( const std::string& regionName, int period );
    double calcSocialSecurityTax( NationalAccount& nationalAccount, const std::string& regionName, int period );
    void calcSavings( double disposableIncome, const std::string& regionName, int period );
    void calcLandSupply( const std::string& regionName, int period );
    void calcLaborSupply( const std::string& regionName, int period );
    void calcIncome( NationalAccount& nationalAccount, const std::string& regionName, int period );
    void calcNoHouseholds( const Demographic* aDemographics, int aPeriod );
    void calcBudget( const std::string& aRegionName, const int aPeriod );
    const std::string getBudgetMarketName() const;
    const std::string getPriceIndexMarketName() const;

    // corresponds to S2, R1, L1
    double baseScalerLand; //!< base land scaler
    double baseScalerLaborMale; //!< base labor scaler
    double baseScalerLaborFemale; //!< base labor scaler
    double baseScalerSavings; //!< base savings scaler

    // corresponds to S0, R0, L0
    double maxLandSupplyFrac; //!< max land supply fraction
    double maxLaborSupplyFracMale; //!< max labor supply fraction
    double maxLaborSupplyFracFemale; //!< max labor supply fraction
    double fixedLaborSupplyFracUnSkLab; //!< fixed unskilled labor supply fraction for both gender
    double fixedLaborSupplyFracSkLab; //!< fixed skilled labor supply fraction for both gender
    double maxSavingsSupplyFrac; //!< maximum savings supply fraction

    double baseLandDemandPerHH; //!< base land demand per household
    double baseLaborDemandPerHH; //!< base labor demand per household
    double landDemand; //!< land demand
    double laborDemand; //!< labor demand
    double householdLandDemand; //!< household land demand
    double householdLaborDemand; //!< household labor demand

    double socialSecurityTaxRate; //!< social security tax rate
    //! land income tax rate
    Value landIncomeTaxRate;

    //! labor income tax rate
    Value laborIncomeTaxRate;

    //! dividends income tax rate
    Value dividendsIncomeTaxRate;

    double personsPerHousehold; //!< people per household
    double numberOfHouseholds; //!< number of households
    double totalLandArea; //!< total land area
    
    double baseTransfer; //!< base year government transfer
    double transfer; //!< government transfer

    double baseLandSupply; //!< base year land supply, in physical units

    double landSupply; //!< land supply
    double laborSupplyMaleUnSkLab; //!< male unskilled labor supply
    double laborSupplyFemaleUnSkLab; //!< female unskilled labor supply
    double laborSupplyUnSkLab; //!< total unskilled labor supply
    double laborSupplyMaleSkLab; //!< male skilled labor supply
    double laborSupplyFemaleSkLab; //!< female skilled labor supply
    double laborSupplySkLab; //!< total skilled labor supply
    
    double mInitialSavings;

    double workingAgePopMale; //!< population of working age males(from Demographics)
    double workingAgePopFemale; //!< population of working age females(from Demographics)
    double workingAgePop; //!< population of working age (from Demographics)

    void copy( const HouseholdConsumer& householdConsumerIn );
};

#endif // _HOUSEHOLD_CONSUMER_H_
