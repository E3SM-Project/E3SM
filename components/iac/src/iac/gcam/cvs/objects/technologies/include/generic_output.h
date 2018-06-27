#ifndef _GENERIC_OUTPUT_H_
#define _GENERIC_OUTPUT_H_
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
 * \file generic_output.h
 * \ingroup Objects
 * \brief GenericOutput class header file.
 * \author Kate Calvin
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>

class Tabs;

#include "technologies/include/primary_output.h"
#include "util/base/include/value.h"

/*! 
 * \ingroup Objects
 * \brief An output used only to drive non-CO2 emissions.
 * \details Generic outputs are used to drive non-CO2 emissions. The 
 *          structure is similar to the PrimaryOutput, but these 
 *          outputs are not added to the market.
 *          
 *          <b>XML specification for GenericOutput</b>
 *          - XML name: None. The generic output is generated automatically
                        and cannot be parsed.
 *
 * \author Kate Calvin
 */
class GenericOutput : public PrimaryOutput {
public:
    /*
     * \brief Constructor
     * \param aSectorName Name of the sector and primary output.
     */
    GenericOutput( const std::string& aSectorName );
    
    virtual ~GenericOutput();

    virtual GenericOutput* clone() const;

    virtual bool isSameType( const std::string& aType ) const;

    virtual const std::string& getXMLReportingName() const;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const int aPeriod );

    virtual void setPhysicalOutput( const double aPrimaryOutput,
                                    const std::string& aRegionName,
                                    ICaptureComponent* aCaptureComponent,
                                    const int aPeriod );

private:
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 
};

#endif // _GENERIC_OUTPUT_H_
