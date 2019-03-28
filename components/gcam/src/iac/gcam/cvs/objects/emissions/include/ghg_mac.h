#ifndef _GHG_MAC_H_
#define _GHG_MAC_H_
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
* \file ghg_mac.h
* \ingroup objects
* \brief The ghg_mac class header file.
* \author Nick Fernandez
*/

#include <xercesc/dom/DOMNode.hpp>
#include <memory>
#include "emissions/include/aghg.h"

class Curve;
class Tabs;


/*! 
* \ingroup objects
* \brief This class represents a Marginal Abatement Curve for a GHG.
* \author Nick Fernandez
*/
class GhgMAC {
public:
    GhgMAC();
    ~GhgMAC();
    GhgMAC( const GhgMAC& other );
    GhgMAC& operator=( const GhgMAC& other );
    GhgMAC* clone() const;

    static const std::string& getXMLNameStatic();
    void XMLParse( const xercesc::DOMNode* node );
    void initCalc( const std::string& ghgName );
    
    double findReduction( const std::string& region, const int period ) const;
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;

protected:
    std::string name; //!< name of mac

    const std::string& getXMLName() const;
    double shiftNatGas( const int period, const std::string& regionName, const double carbonPrice) const;
    double adjustPhaseIn( const int period ) const;
    double adjustTechCh( const int period, const int finalReductionPeriod, const double maxReduction ) const;
    double shiftCostReduction( const int period, const double costReductionRate ) const;
    double getMACValue( const double carbonPrice ) const;

    //! number of periods over which phase in occurs. can be a non-integer
    double phaseIn;
    
    //! the initial range over which carbon price changes due to the standard range of Nat. Gas price changes
    double fuelShiftRange;
    
    //! the annual rate at which carbon price is shifted due to technological change (Direct number. i.e. 1% = 0.01)
    double costReductionRate; 
    
    //! Increase maximum reduction to this value (due to tech change)  (Direct number. i.e. 85% max reduction = 0.85)
    double finalReduction; 
    
    //! Base year from which to start decreasing costs
    int baseCostYear; 
    
    //! Year in which maximum reduction should be implemented
    int finalReductionYear; 
    
    //! turns off reductions if carbon Price is less than 0;
    bool noBelowZero;
    
    //! Name of fuel who's price changes cause a shift in the curve
    std::string curveShiftFuelName; 
    
    //! The underlying Curve
    std::auto_ptr<Curve> macCurve; 
    
private:
    static const std::string XML_NAME; //!< node name for toXML methods
    void copy( const GhgMAC& other );
};

#endif // _GHG_MAC_H_
