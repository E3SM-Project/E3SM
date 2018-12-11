#ifndef _GENDER_H_
#define _GENDER_H_
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
* \file gender.h
* \ingroup Objects
* \brief The Gender class header file.
* \author Sonny Kim
* \author Katherine Chung
*/

#include <xercesc/dom/DOMNode.hpp>
#include <string>
#include "util/base/include/iround_trippable.h"
#include "util/base/include/ivisitable.h"

/*! 
* \ingroup Objects
* \brief The base class for Male and Female objects. Both have populations and survival rates.
*/

class Gender: public IRoundTrippable, IVisitable {
    friend class XMLDBOutputter; // For getXMLName()
public:
	Gender();
	void XMLParse( const xercesc::DOMNode* node );
	void toInputXML( std::ostream& out, Tabs* tabs ) const;
	void toDebugXML( std::ostream& out, Tabs* tabs ) const;
	static const std::string& getXMLNameStatic();
	double calcSurvivingPop();
	double getPopulation() const;
	void setPopulation( double aPopulation );
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const = 0;
protected:
    virtual void toDebugXMLDerived( std::ostream& out, Tabs* tabs ) const = 0;
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const = 0;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ) = 0;
    virtual const std::string& getXMLName() const = 0;
    static const std::string XML_NAME; //!< node name for toXML methods
    double mPopulation; //!< population for this gender
    double mSurvivalRate; //!< survival rate
    double mSurvivingPop; //!< surviving Population, calculated
};

#endif // _GENDER_H_
