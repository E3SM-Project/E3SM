#ifndef _SO2_EMISSIONS_H_
#define _SO2_EMISSIONS_H_
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
* \file so2_emissions.h
* \ingroup Objects
* \brief SO2Emissions class header file.
* \author
*/

#include "emissions/include/acomplex_emissions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>

/*! 
 * \ingroup Objects
 * \brief A class that represents SO2 emissions.
 * \warning SO2Emissions will not work with an emissions coefficient of a type
 *          other than ReadEmissionsCoef.
 * \todo Currently this method overrides initCalc and sets the emissions coefficient
 *       delegate to 1.  The emissions coefficient is contained in the emissionsDriver
 *       function.  This is a hack and should be refactored.
 * \author Jim Naslund
 */
class SO2Emissions: public AComplexEmissions { 
public:
    SO2Emissions();

    SO2Emissions* clone() const;
    static const std::string& getXMLNameStatic();
    virtual const std::string& getName() const;

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aFuelName,
                           const IInfo* aLocalInfo,
                           const int aPeriod );

protected:
    virtual const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void parseName( const std::string& aNameAttr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    void setAdjMaxCntrl();
    virtual void adjustMaxCntrl(const double GDPcap);
    virtual double emissionsDriver( const double inputIn, const double outputIn ) const;
    
    double ashRetention; //!< percentage of output that remains in ash form
    double percentSulfur; //!< sulfur content of input (percentage)  
    double gjPerTonne; //!< fuel energy(in GJ) per metric ton
    double finalSulfur; // Asymptotic final sulfur content (percentage) 
};

#endif // _SO2_EMISSIONS_H_

