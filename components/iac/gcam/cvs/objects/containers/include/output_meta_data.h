#ifndef _OUTPUT_META_DATA_H_
#define _OUTPUT_META_DATA_H_
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
* \file output_meta_data.h
* \ingroup Objects
* \brief OutputMetaData header file.
* \author Josh Lurz
*/
#include <list>
#include <string>

#include "util/base/include/ivisitable.h"
#include "util/base/include/iround_trippable.h"

/*! 
* \ingroup Objects
* \brief A container of read-in metadata used for outputting information and
*        passed to the dataviewer.
* \author Josh Lurz
*/

class OutputMetaData: public IVisitable, public IRoundTrippable
{
    friend class IVisitor;
public:
    OutputMetaData();
    static const std::string& getXMLNameStatic();
    void XMLParse( const xercesc::DOMNode* aNode );
    void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    const std::list<std::string>& getPrimaryFuelList() const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;
private:
    //! List of names of primary fuels.     
    std::list<std::string> mPrimaryFuels;

    //! List of variables which are summable.
    std::list<std::string> mSummableVariables;

    //! List of variables which have year attribute.
    std::list<std::string> mHasYearVariables;

    //! Scenario summary.
    std::string mScenarioSummary;
};

#endif // _OUTPUT_META_DATA_H_

