#ifndef _REGION_CGE_H_
#define _REGION_CGE_H_
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
* \file region_cge.h
* \ingroup Objects-SGM
* \brief The RegionCGE class header file.
* \author Sonny Kim
*/

#include <map>
#include <iosfwd>

// User headers
#include "containers/include/region.h"

// Forward declare headers
class ProductionSector;
class FinalDemandSector;
class FactorSupply;
class IVisitor;
class NationalAccount;
// TEMP
class SGMGenTable;
class Tabs;

/*!
 * \ingroup Objects-SGM
 * \brief This derived Region class contains SGM specific information.
 * \details A RegionCGE basically consists of ProductionSectors, FinalDemandSectors,
 *          FactorSupplies, and Resources.  They also contain NationalAccounts for regional
 *          accounting of value flows and various reporting structures.  For SGM these
 *          steps are completed in the following order to operate a region:
 *           - Set the price of the numeraire good as the price of the price-index.  This
 *             is nescessary because there is no direct link between the two markets and 
 *             it must be taken care of before any of the prices are taken into consideration.
 *             Note that there may be a better way to ensure that the price of the numeraire is
 *             equal to the price-index however this was the easiest solution.
 *           - Calculate the price of the Capital good.  This would be the CES aggregate price
 *             of the goods that make up the Capital good which is contained in the investment
 *             consumer.  This must be taken care of before new investment demands are calculated
 *             by production sectors.
 *           - Operate all production sectors.  This includes domestic and import sectors.  It will
 *             operate all old vintage technologies, distrubute the level of new investment, and operate
 *             the new vintage technologies.  This will calculate the supply of all the regular goods.
 *           - Operate all consumers which includes government, housholds, investment, and trade.  The
 *             consumers will provide the supply for the factor goods besides resources.
 *           - The last step is to calculate the supply for all resources.
 *
 *          Note that at the moment FactorSupply is barely used and could possibly be removed or moved
 *          into the appropriate consumer.
 *
 *          <b>XML specification for RegionCGE</b>
 *          - XML name: \c regionCGE
 *          - Contained by: World
 *          - Parsing inherited from class: Region
 *          - Attributes:
 *          - Elements:
 *              - \c FinalDemandSector::getXMLNameStatic() RegionCGE::finalDemandSector
 *              - \c FactorSupply::getXMLNameStatic() RegionCGE::factorSupply
 *              - \c ProductionSector::getXMLNameStatic() Region::supplySector
 *              - \c NationalAccount::getXMLNameStatic() RegionCGE::mNationalAccounts
 *
 * \author Sonny Kim
 */
class RegionCGE : public Region
{
    friend class SocialAccountingMatrix;
    friend class DemandComponentsTable;
    friend class SectorReport;
    friend class SGMGenTable;
    friend class InputOutputTable;
    friend class XMLDBOutputter;
public:
    RegionCGE();
    ~RegionCGE(); 
    static const std::string& getXMLNameStatic();
    virtual void completeInit();
    virtual void initCalc( const int period);
    virtual void postCalc( const int aPeriod );
    virtual void calc( const int period );
    virtual void updateMarketplace( const int period );
    virtual void updateAllOutputContainers( const int period );
    virtual void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
    virtual void csvSGMGenFile( std::ostream& aFile ) const;

protected:
    const static std::string XML_NAME; //!< node name for toXML method.
    std::vector<FinalDemandSector*> finalDemandSector; //!< vector of pointers to final demand sector objects
    std::vector<FactorSupply*> factorSupply; //!< vector of pointers to factor supply objects
    std::vector<SGMGenTable*> mOutputContainers; //!< vector of output containers
    std::vector<NationalAccount*> mNationalAccounts; //!< vector of NationalAccounts, one for each period.
    IVisitor* mCalcCapitalGoodPriceVisitor; //!< vistor to calculate the price of the capital good

    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    void operate( const int period );
private:
    void createSGMGenTables();
    void clear();
};

#endif // _REGION_CGE_H_
