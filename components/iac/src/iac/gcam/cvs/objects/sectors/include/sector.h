#ifndef _SECTOR_H_
#define _SECTOR_H_
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
* \file sector.h
* \ingroup Objects
* \brief The Sector class header file.
* \author Sonny Kim
*/

#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include <map>
#include <memory>
#include <list>

#include "containers/include/national_account.h" // lets use an auto_ptr instead.
#include "util/base/include/ivisitable.h"
#include "util/base/include/iround_trippable.h"
#include "util/base/include/inamed.h"
#include "util/base/include/object_meta_info.h"
#include "util/base/include/time_vector.h"

// Forward declarations
class Subsector;
class Summary;
class ILogger;
class GDP;
class Tabs;
class IInfo;
class DependencyFinder;
class Demographic;
class NationalAccount;
class MoreSectorInfo;
class SocialAccountingMatrix;
class ILandAllocator;
class IndirectEmissionsCalculator;
class AGHG;

/*! 
* \ingroup Objects
* \brief This class represents a single good that is produced, transformed, or consumed.
*
* All production, consumption, and transformation (other than resource extraction)
* is contained within the Sector class. Each Sector represents a distinct good that
* can either be supplied or demanded. The demand Sector derived from this class contains
* a few classes where changes are necessary, although most of the basic mechanisms are unchanged.
*
* \author Sonny Kim, Steve Smith, Josh Lurz
*/

class Sector: public IVisitable,
              public IRoundTrippable,
              public INamed
{
    // TODO: Remove the need for these.
    friend class SocialAccountingMatrix;
    friend class DemandComponentsTable;
    friend class SectorReport;
    friend class SGMGenTable;
    friend class XMLDBOutputter;
    // TODO: Replace with visitor.
    friend class TotalSectorEmissions;
    friend class CalQuantityTabulator;
    friend class CalibrateShareWeightVisitor;
protected:
    std::string name; //!< Sector name
    std::string mOutputUnit; //!< unit of good or service produced by sector
    std::string mInputUnit; //!< unit of input demanded by sector
    std::string mPriceUnit; //!< price unit of good or service produced by sector
    std::string regionName; //!< region name

    //! Type of the sector.
    std::string mSectorType;

    double mBaseOutput; //!< Read in base year output. TODO: Move to demand sector.
    std::auto_ptr<IInfo> mSectorInfo; //!< Pointer to the sector's information store.
    std::vector<Subsector*> subsec; //!< subsector objects
    typedef std::vector<Subsector*>::iterator SubsectorIterator;
    typedef std::vector<Subsector*>::const_iterator CSubsectorIterator;
    
    //! Sector price in $/service. TODO: Move to supply and production sector.
    double mBasePrice;
    //! Sector price by period updated with solution prices.
    std::vector<double> mPrice;

    std::vector<Summary> summary; //!< summary for reporting
    std::map<std::string,int> subSectorNameMap; //!< Map of subSector name to integer position in vector.
    std::auto_ptr<MoreSectorInfo> moreSectorInfo; //! Additional sector information needed below sector

    typedef ObjECTS::TObjectMetaInfo<> object_meta_info_type;
    typedef std::vector<object_meta_info_type> object_meta_info_vector_type;
    object_meta_info_vector_type mObjectMetaInfo; //!< Vector of object meta info to pass to mSectorInfo

    //! A map of a keyword to its keyword group
    std::map<std::string, std::string> mKeywordMap;

    //! Logit exponential used for the subsector competition.
    objects::PeriodVector<double> mSubsectorLogitExp;

    static double getDefaultSubsectorLogitExp();

    virtual void toInputXMLDerived( std::ostream& aOut, Tabs* aTabs ) const = 0;
    virtual void toDebugXMLDerived( const int period, std::ostream& aOut, Tabs* aTabs ) const = 0;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ) = 0;
    virtual const std::string& getXMLName() const = 0;
    
    double getFixedOutput( const int period ) const;
    const std::vector<double> calcSubsectorShares( const GDP* aGDP, const int aPeriod ) const;
    static const std::string& getDefaultSectorType();
    const std::string& getSectorType() const;

    bool outputsAllFixed( const int period ) const;
    
    double getCalOutput( const int period ) const;

    virtual double getPrice( const GDP* aGDP, const int aPeriod ) const;

public:
    explicit Sector( const std::string& aRegionName );
    virtual ~Sector();
    virtual const std::string& getName() const;

    virtual void XMLParse( const xercesc::DOMNode* node );
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;

    virtual void completeInit( const IInfo* aRegionInfo,
                               DependencyFinder* aDepFinder,
                               ILandAllocator* aLandAllocator ) = 0;

    virtual void initCalc( NationalAccount* aNationalAccount,
                           const Demographic* aDemographics,
                           const int aPeriod ) = 0;

    bool isAllCalibrated( const int period, double calAccuracy, const bool printWarnings ) const;

    virtual void supply( const GDP* aGDP,
                         const int aPeriod ) = 0;

    void calcCosts( const int aPeriod );

    virtual double getOutput( const int period ) const = 0;

    virtual void calcFinalSupplyPrice( const GDP* aGDP, const int aPeriod ) = 0;

    void emission( const int period );
    
    virtual void csvOutputFile( const GDP* aGDP,
                                const IndirectEmissionsCalculator* aIndirectEmissCalc ) const;

    virtual void dbOutput( const GDP* aGDP,
                           const IndirectEmissionsCalculator* aIndEmissCalc ) const = 0;

    std::map<std::string, double> getfuelcons( const int period ) const;
    double getConsByFuel( const int period, const std::string& key) const;
    std::map<std::string, double> getemission( const int period ) const;
    std::map<std::string, double> getemfuelmap( const int period ) const;
    void updateSummary( const std::list<std::string>& aPrimaryFuelList, const int period );

    virtual void operate( NationalAccount& nationalAccount, const Demographic* aDemographic, const int period ) = 0;
    void updateMarketplace( const int period );
    virtual void postCalc( const int aPeriod );

    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
private:
    void clear();
};

#endif // _SECTOR_H_
