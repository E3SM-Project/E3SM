#ifndef _SUBRSRC_H_
#define _SUBRSRC_H_
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
* \file subresource.h
* \ingroup Objects
* \brief The SubResource class header file.
* \author Sonny Kim
*/
#include <memory>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/ivisitable.h"
#include "util/base/include/value.h"
#include "util/base/include/time_vector.h"

// Forward declarations.
class Grade;
class GDP;
class IInfo;
class IVisitor;
class Tabs;

/*! 
* \ingroup Objects
* \brief SubResource is a class that contains grades.
* \author Sonny Kim
*/

class SubResource: public IVisitable
{
	friend class XMLDBOutputter;
    friend class CalibrateResourceVisitor;    
    friend class EnergyBalanceTable;
public:
    SubResource();
    virtual ~SubResource();
    std::string getName() const;
    void XMLParse( const xercesc::DOMNode* tempnode );
    virtual void completeInit( const IInfo* aResourceInfo );
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    static const std::string& getXMLNameStatic();
    virtual void cumulsupply( double aPrice, int aPeriod );
    virtual void initCalc( const std::string& aRegionName, const std::string& aResourceName, const int aPeriod );
    virtual void postCalc( const std::string& aRegionName, const std::string& aResourceName, const int aPeriod );
    double getCumulProd( const int aPeriod ) const;
    virtual void annualsupply( int aPeriod, const GDP* aGdp, double aPrice, double aPrevPrice );
    double getAnnualProd( int aPeriod ) const;
    double getAvailable( int aPeriod ) const;
    void dbOutput( const std::string& regname, const std::string& secname); 
    void csvOutputFile(const std::string &regname, const std::string& sname); 
    void updateAvailable( const int period );
    virtual double getVariance() const;
    virtual double getAverageCapacityFactor() const;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
protected:
    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* node ) = 0;
    virtual void toXMLforDerivedClass( std::ostream& out, Tabs* tabs ) const;

    std::string mName; //!< SubResource name
    std::auto_ptr<IInfo> mSubresourceInfo; //!< The subsector's information store.
    std::vector<double> mAvailable; //!< total available resource
    std::vector<double> mAnnualProd; //!< annual production of SubResource
    std::vector<double> mCumulProd; //!< cumulative production of SubResource
    std::vector<double> mCumulativeTechChange; //!< Cumulative Technical Change for this subsector
    std::vector<double> mEffectivePrice; //!< effective price (global price + price adder)
    std::vector<double> mCalProduction; //!< calibrated production 

    objects::PeriodVector<Value> mTechChange; //!< technical change
    objects::PeriodVector<Value> mEnvironCost; //!< Environmental costs
    objects::PeriodVector<Value> mSeveranceTax; //!< Severance Tax (exogenous)
    objects::PeriodVector<Value> mPriceAdder; //!< price adder used for calibration purposes
    
    // Cumulative technical change needs to be in sub-resource sector 
    std::vector<Grade*> mGrade; //!< amount of SubResource for each grade
    std::map<std::string,int> mGradeNameMap; //!< Map of grade name to integer position in vector. 

private:
    static const std::string XML_NAME; //!< node name for toXML methods
};


/*! 
* \ingroup Objects
* \brief A class which defines a SubDepletableResource object, which is a container for multiple grade objects.
* \author Steve Smith
*/
class SubDepletableResource: public SubResource {
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* node );
};

/*! 
* \ingroup Objects
* \brief A class which defines a SubFixedResource object, which is a container for multiple grade objects.
* \author Steve Smith
*/
class SubFixedResource: public SubResource {
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* node );
};

#endif // _SUBRSRC_H_
