#ifndef _MORE_SECTOR_INFO_H_
#define _MORE_SECTOR_INFO_H_
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
* \file more_sector_info.h
* \ingroup Objects
* \brief MoreSectorInfo class header file.
* \author Sonny Kim
*/

#include <string>
#include <vector>
#include <map>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/iround_trippable.h"

/*! 
* \ingroup Objects
* \brief A container which contains additional information about an SGM sector.
* \details TODO
* \author Sonny Kim
*/

class MoreSectorInfo: public IRoundTrippable
{
public:
    enum MoreSectorInfoType {
        ENERGY_CURRENCY_CONVERSION,
        INVEST_TAX_CREDIT_RATE,
        CORP_INCOME_TAX_RATE,
        IND_BUS_TAX_RATE,
        MAX_CORP_RET_EARNINGS_RATE,
        CORP_RET_EARNINGS_RATE,
        HH_RET_EARNINGS_RATE,
        RET_EARNINGS_PARAM,
		TRANSPORTATION_COST,
		TRAN_COST_MULT,
		PRICE_ADJUST_MULT,
		PROPORTIONAL_TAX_RATE,
		ADDITIVE_TAX
    };

	MoreSectorInfo();
    virtual ~MoreSectorInfo();
	void XMLParse( const xercesc::DOMNode* node );
	void toInputXML( std::ostream& out, Tabs* tabs ) const;
	void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
	static const std::string& getXMLNameStatic();
	void reset();
	void setType( const MoreSectorInfoType aType, const double aValue );
    double getValue( const MoreSectorInfoType aType ) const;

private:
    const std::string& getXMLName() const;
    const std::string enumToName( const MoreSectorInfoType aType ) const;
    const std::string enumToXMLName( const MoreSectorInfoType aType ) const;
	static const std::string XML_NAME; //!< node name for toXML methods
	std::map<MoreSectorInfoType, double> mSectorInfoMap; //!< Map relating additional sector info by type to value
    typedef std::map<MoreSectorInfoType, double>::const_iterator CInfoTypeIterator;
};

#endif // _MORE_SECTOR_INFO_H_
