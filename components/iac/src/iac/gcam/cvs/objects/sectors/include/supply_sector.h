#ifndef _SUPPLY_SECTOR_H_
#define _SUPPLY_SECTOR_H_
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
* \file supply_sector.h
* \ingroup Objects
* \brief The SupplySector class header file.
* \author James Blackwood
*/
#include <string>
#include "sectors/include/sector.h"
class NationalAccount;
class IInfo;
class DependencyFinder;
/*! 
* \ingroup Objects
* \brief This class represents a single supply sector.
* \author James Blackwood, Sonny Kim
*/
class SupplySector: public Sector
{
public:
    explicit SupplySector( const std::string& aRegionName );
    virtual ~SupplySector(){};
    static const std::string& getXMLNameStatic();
    
    virtual void completeInit( const IInfo* aRegionInfo,
                               DependencyFinder* aDepFinder,
                               ILandAllocator* aLandAllocator );

    
    virtual void initCalc( NationalAccount* aNationalAccount,
                           const Demographic* aDemographics,
                           const int aPeriod );

    virtual void calcFinalSupplyPrice( const GDP* aGDP, const int aPeriod );
    
    virtual void supply( const GDP* aGDP,
                         const int aPeriod );

    virtual void operate( NationalAccount& aNationalAccount, const Demographic* aDemographic,
                          const int aPeriod ){};

    virtual void postCalc( const int aPeriod );

    virtual void dbOutput( const GDP* aGDP,
                           const IndirectEmissionsCalculator* aIndEmissCalc ) const;
protected:
    
    /*!
     * \brief Class responsible for setting final energy into the calibration
     *        market.
     */
    class FinalEnergySupplier {
    public:
        FinalEnergySupplier( const std::string& aSectorName );

        void setFinalEnergy( const std::string& aRegionName,
                             const double aFinalEnergy,
                             const int aPeriod );
    private:
        //! The cached name of the TFE market.
        std::string mTFEMarketName;
    };

    virtual double getEnergyInput( const int aPeriod ) const;
    virtual double getOutput( const int aPeriod ) const;
    virtual double getPrice( const GDP* aGDP, const int aPeriod ) const;
    virtual void setMarket();
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ); 

    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;

    virtual const std::string& getXMLName() const;

    //! Temporary hack for CCTP. Biomass adder.
    std::vector<double> mBiomassAdder;

    //! The object responsible for setting final energy supply into the
    //! calibration market.
    std::auto_ptr<FinalEnergySupplier> mFinalEnergySupplier;

    //! Whether the sector has a trial supply market.
    bool mHasTrialSupplyMarket;
    //! Trial supply market prices
    std::vector<double> mPriceTrialSupplyMarket;

private:
    const static std::string XML_NAME; //!< node name for toXML methods 
};

#endif // _SUPPLY_SECTOR_H_
