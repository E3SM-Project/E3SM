#ifndef _FACTOR_SUPPLY_H_
#define _FACTOR_SUPPLY_H_
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
* \file factor_supply.h
* \ingroup Objects
* \brief FactorSupply class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <vector>
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <memory>
#include "util/base/include/ivisitable.h"
#include "util/base/include/iround_trippable.h"

class MoreSectorInfo;

/*! 
* \ingroup Objects
* \brief The supplier of the initial goods such as land, labor and capital used
*        to create intermediate products,
* \details
* \todo This object is barely used, it may be possible to get rid of it and move
*       it's current functionality into the consumers, or move the object into 
*       the consumers where it could take over some of the factory supply duties
*       the consumers currently take care of.
* \author Pralit Patel, Sonny Kim
*/

class FactorSupply: public IVisitable, public IRoundTrippable
{
    friend class SGMGenTable;
    friend class XMLDBOutputter;
public:
    FactorSupply();
    const std::string& getName() const;
    bool XMLParse( const xercesc::DOMNode* node );
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;

    static const std::string& getXMLNameStatic();
    void completeInit( const std::string& aRegionName );
    void initCalc( const std::string& aRegionName, const int period );

    double getSupply( const std::string& aRegionName, const int period ) const;
    void calcPricePaid( const std::string& aRegionName, const int period );
    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;

protected:
    void setMarket( const std::string& aRegionName );
    const std::string& getXMLName() const;
    void setMarket();
    std::string name; //!< CHANGE
    std::string marketName; //!< regional market. this isn't read in right now.
    double mBasePrice; //!< Price in the base year. only applied to capital, should be in derived class.
    double mBaseSupply; //!< Supply in the base year, only applies to land and labor.
private:
    static const std::string XML_NAME; //!< node name for toXML methods
    std::auto_ptr<MoreSectorInfo> moreSectorInfo; //! Additional sector information needed below sector
};

#endif // _FACTOR_SUPPLY_H_
