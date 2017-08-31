#ifndef _AEMISSIONS_COEF_H_
#define _AEMISSIONS_COEF_H_
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
 * \file aemissions_coef.h
 * \ingroup Objects
 * \brief AEmissionsCoef header file.
 * \author Jim Naslund
 */

class IInfo;
class Tabs;

/*! 
 * \ingroup Objects
 * \brief An abstract class that encapsulates the emissions coefficient.
 * \details This class encapsulates the emissions coefficient.  Its subclasses
 *          determine how the emissions coefficient is calculated.
 *          This class provides empty implementations of setCoef, updateCoef,
 *          and initCalc because derived classes may not need to use them.
 * \author Jim Naslund
 */
class AEmissionsCoef{

public:
    AEmissionsCoef( const double aEmissionsCoef = 0 );
    virtual AEmissionsCoef* clone() const = 0;

    double getCoef() const;
    virtual void setCoef( const double aEmissionsCoef );
    virtual void updateCoef( const double aOutput );
    /*!
     * \brief Performs any initializations that only need to be done once per period.
     * \param aSubSectorInfo Pointer to the IInfo object
     * \param aName The name of the region
     */
    virtual void initCalc( const IInfo* aSubsectorInfo, const std::string& aName, const int aPeriod ) = 0;
    virtual double getEmissions( const double aOutput ) const;
    virtual double getInputEmissions() const;
    virtual double calcMaxCntrl( const double aFinalEmissCoef, const double aB,
                                 const double aMultiplier ) const;
    virtual bool needsCalcForAdjustment() const;
    bool getOverride() const;
    virtual void overrideCoef( const double aEmissionsCoef );

    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( std::ostream& out, Tabs* tabs ) const;
    /*!
     * \brief Get the XML node name for output to XML.
     * \details This public function accesses the XML_NAME of the coefficient.
     *          This way the tag is always consistent for both read-in and output and can be easily changed.
     *          This function may be virtual to be overridden by derived class pointers.
     * \author Jim Naslund
     * \return The constant XML_NAME.
     */
    virtual const std::string& getXMLName() const = 0;
    virtual double getXMLValue() const;

protected:
    double mEmissionsCoef;
    bool mOverridesFuturePeriods; //!< Whether this period's emissCoef overrides future period's emissCoef
  
};


#endif // _AEMISSIONS_COEF_H_
