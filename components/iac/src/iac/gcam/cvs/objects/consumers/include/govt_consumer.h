#ifndef _GOVT_CONSUMER_H_
#define _GOVT_CONSUMER_H_
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
* \file govt_consumer.h
* \ingroup Objects
* \brief GovtConsumer class header file.
*
*  Detailed description.
*
* \author Katherine Chung
* \author Sonny Kim
*/

#include <string>
#include <xercesc/dom/DOMNode.hpp>

#include "consumers/include/consumer.h"
#include "util/base/include/value.h"

class NationalAccount;
class Demographic;
class Tabs;
class MoreSectorInfo;

/*! 
 * \brief A consumer representing government purchases.
 */
class GovtConsumer : public Consumer
{
    friend class SocialAccountingMatrix;
    friend class DemandComponentsTable;
    friend class SectorReport;
    friend class SGMGenTable;
    friend class GovtResults;
    friend class XMLDBOutputter;
public:
    GovtConsumer();
    GovtConsumer* clone() const;

    void copyParam( const BaseTechnology* baseTech,
                    const int aPeriod );

    void copyParamsInto( GovtConsumer& govtConsumerIn,
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

    void operate( NationalAccount& aNationalAccount, const Demographic* aDemographics, 
        const MoreSectorInfo* aMoreSectorInfo, const std::string& aRegionName,
        const std::string& aSectorName, const bool isNewVintageMode, const int aPeriod );

    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;
    static const std::string& getXMLNameStatic();
protected:
    virtual bool isCoefBased() const { return true; }
    void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    const std::string& getXMLName() const;
private:
    void calcTotalTax( NationalAccount& aNationalAccount, const std::string& aRegionName,
                       const int aPeriod );
    void calcTransfer( NationalAccount& nationalAccount, const Demographic* demographics, 
        const std::string& regionName, int period );
    void calcBaseCoef( NationalAccount& nationalAccount, const Demographic* aDemographics );
    void calcGovtTaxOrSubsidy( const std::string& regionName, int period );
    void calcGovtCapitalDemand( const std::string& regionName, int period ); // not currently used
    void calcSubsidy( NationalAccount& nationalAccount, const std::string& regionName, int period );
    void calcDeficit( const std::string& regionName, int period );
    void calcIncome( NationalAccount& nationalAccount, const Demographic* demographics, 
        const std::string& regionName, int period );
    
    void calcBudget();
    
    // Read in parameters.
    Value mBaseDeficit; //!< Base Real deficit
    Value mBaseTransfer; //!< Read in base year real transfer   
    Value mTaxProportional; //!< Proportional Tax
    Value mTaxAdditive; //!< Additive Tax
    Value mTaxCorporate; //!< Corporate income tax
    Value mTaxIBT; //!< Indirect business tax
    Value mRho;

    // Calculated parameters.
    Value mBaseTransferPopCoef;
    Value mSigma;
};

#endif // _GOVT_CONSUMER_H_
