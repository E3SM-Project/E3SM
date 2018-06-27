#ifndef _ICAPTURE_COMPONENT_H_
#define _ICAPTURE_COMPONENT_H_
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
 * \file icapture_component.h
 * \ingroup Objects
 * \brief ICaptureComponent interface header file.
 * \author Josh Lurz
 */
#include <vector>
#include "util/base/include/istandard_component.h"
class DependencyFinder;
class IInput;

/*! 
 * \ingroup Objects
 * \brief This object is responsible for determining the quantity and cost of
 *        sequestered emissions for a Technology.
 * \details This object can be added on to Technologies so that they can
 *          sequester their emissions instead of emitting them. The capture
 *          component is responsible for determining the fraction of emissions
 *          captured, the cost and efficiency loss for capturing the emissions,
 *          and how the emissions are disposed.
 * \author Josh Lurz
*/
class ICaptureComponent : public IParsedComponent { 
public:
    // Clone operator must be declared explicitly even though it is inherited
    // from IStandardComponent so that the return type can be changed. Since
    // this class is a subtype of IStandardComponent, this is legal and referred
    // to as a covariant return type.
    virtual ICaptureComponent* clone() const = 0;

    /*! \brief Returns whether the type of the object is the same as the passed
    *          in type.
    * \param aType Type to check the object's type against.
    * \return Whether the type of the object is the same as the passed in type.
    */
    virtual bool isSameType( const std::string& aType ) const = 0;
    
    /*! \brief Parse the data for this object starting at a given node.
    * \param aNode Root node from which to parse data.
    */
    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;
    
    /*! \brief Write data from this object in an XML format so that it can be
    *          read back in later as input.
    * \param aOut Filestream to which to write.
    * \param aTabs Object responsible for writing the correct number of tabs. 
    */
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const = 0;
    
    /*! \brief Write data from this object in an XML format for debugging.
    * \param aPeriod Period for which to write data.
    * \param aOut Filestream to which to write.
    * \param aTabs Object responsible for writing the correct number of tabs. 
    */
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const = 0;
    
    /*! \brief Complete the initialization of the capture component.
    * \param aRegionName Region name.
    * \param aSectorName Sector name.
    * \param aDependencyFinder Regional dependency finder.
    */
    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               DependencyFinder* aDependencyFinder ) = 0;
    
    /*!
     * \brief Initialize the capture component for a given period.
     * \param aRegionName Region name.
     * \param aSectorName Sector name.
     * \param aFuelName Name of the fuel being consumed.
     * \param aPeriod Model period.
     */
    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const std::string& aFuelName,
                           const int aPeriod ) = 0;

    /*!
     * \brief Get the cost of storing one unit of carbon.
     * \details Calculates and returns the cost of storing a unit of carbon.
     *          This may either be based on a storage market or a read-in value.
     *          There is no switch to turn off sequestration if there is no tax,
     *          as this causes the calculation of the total social cost of a
     *          policy to be incorrect.
     * \param aRegionName Name of the region containing the capture component.
     * \param aPeriod Period in which the carbon is being stored.
     * \return The cost of storing one unit of carbon.
     */
    virtual double getStorageCost( const std::string& aRegionName,
                                   const std::string& aGHGName,
                                   const int aPeriod ) const = 0;
    
    /*! 
     * \brief Get the fraction of emissions captured by the capture component.
     * \return The fraction of emissions captured.
     */
    virtual double getRemoveFraction( const std::string& aGHGName ) const = 0;
    
    /*!
     * \brief Calculate  the amount of emissions that are sequestered.
    * \param aRegionName Name of the region.
    * \param aGHGName Name of the GHG stored.
    * \param aTotalEmissions Total technology emissions before capture occurs.
    * \param aPeriod Model period.
    */
    virtual double calcSequesteredAmount( const std::string& aRegionName,
                                          const std::string& aGHGName,
                                          const double aTotalEmissions,
                                          const int aPeriod ) = 0;
    
    /*!
     * \brief Get the amount of the emissions sequestered.
     * \param aGHGName Name of the GHG to capture.
     * \param aGetGeologic Whether to get geologically sequestered emissions.
     * \param aPeriod Period for which to get emissions.
     * \pre The sequestered amount has been calculated.
     * \pre calcSequesteredAmount
     */
    virtual double getSequesteredAmount( const std::string& aGHGName, 
                                         const bool aGetGeologic,
                                         const int aPeriod ) const = 0;
    
    /*! \brief Adjust the technology's inputs for the effects of the capture
    *          component.
    * \details Capture components cause technologies to incur additional capital
    *          costs and decrease the efficiency of fuel usage.
    * \param aInputs List of inputs to adjust.
    * \param aPeriod Model period.
    */
    virtual void adjustInputs( const std::string& aRegionName,
                               std::vector<IInput*>& aInputs,
                               const int aPeriod ) const = 0;
};

#endif // _ICAPTURE_COMPONENT_H_
