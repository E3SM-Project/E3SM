#ifndef _RES_SECONDARY_OUTPUT_H_
#define _RES_SECONDARY_OUTPUT_H_
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
 * \file secondary_output.h
 * \ingroup Objects
 * \brief SecondaryOutput class header file.
 * \author Josh Lurz, Patrick Luckow
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>

#include "technologies/include/secondary_output.h"

/*! 
 * \ingroup Objects
 * \brief An output representing a secondary output of the Technology.
 * \details Secondary outputs represent a product produced in conjunction with a
 *          primary product which may have a positive value. Output of the
 *          secondary output is not optimized, it is always produced at a fixed
 *          ratio to the primary output, and may not be discarded. The revenue
 *          from the secondary output must not exceed the cost of the primary
 *          output. The secondary output must be the primary output of another
 *          sector, and the total secondary output of the good must not exceed
 *          demand.
 *          
 *          <b>XML specification for SecondaryOutput</b>
 *          - XML name: \c secondary-output
 *          - Contained by: Technology
 *          - Parsing inherited from class: None.
 *          - Attributes:
 *              - \c name SecondaryOutput::mName
 *          - Elements:
 *              - \c output-ratio SecondaryOutput::mOutputRatio
 *          
 * \author Josh Lurz, Patrick Luckow
 */
class RESSecondaryOutput: public SecondaryOutput
{
    friend class OutputFactory;
public:
    /*!
     * \brief Get the XML name for the class.
     * \return The XML name for the class.
     */
    static const std::string& getXMLNameStatic();

    virtual RESSecondaryOutput* clone() const;

    virtual bool isSameType( const std::string& aType ) const;


    virtual const std::string& getXMLReportingName() const;

    
    virtual void setPhysicalOutput( const double aPrimaryOutput,
                                   const std::string& aRegionName,
                                   ICaptureComponent* aCaptureComponent,
                                   const int aPeriod );

protected:
    /*!
     * \brief Protected constructor so the class can only be created by the
     *        OutputFactory.
     */
    RESSecondaryOutput();
    

private:
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 
};

#endif // _RES_SECONDARY_OUTPUT_H_
