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
 * wind_technology.h
 * Created: 03/20/2007
 * Version: 03/20/2007
 *
 * This software, which is provided in confidence, was prepared by employees
 * of Pacific Northwest National Laboratory operated by Battelle Memorial
 * Institute. Battelle has certain unperfected rights in the software
 * which should not be copied or otherwise disseminated outside your
 * organization without the express written authorization from Battelle.
 * All rights to the software are reserved by Battelle.   Battelle makes no
 * warranty, express or implied, and assumes no liability or responsibility
 * for the use of this software.
 */

#if !defined( __WIND_TECHNOLOGY_H )
#define __WIND_TECHNOLOGY_H      // prevent multiple includes

// include files ***********************************************************

#include "technologies/include/intermittent_technology.h"

// namespaces **************************************************************

// class: WindTechnology ***************************************************

/*!
 * \ingroup objects::biomass
 * \brief A technology class for wind power
 * \details This class contains a set of routines that implement
 *          technology for wind power
 *
 *   <b>XML specification for WindTechnology</b>
 *   - XML name: \c wind-technology
 *   - Contained by: Technology
 *   - Parsing inherited from class: None.
 *   - Attributes: none
 *   - Elements:
 *   - \c  WindTechnology::
 *
 * \author Kevin Walker
 * \date $ Date $
 * \version $ Revision $
 */
class WindTechnology : public IntermittentTechnology
{
public :

   typedef IntermittentTechnology   parent;

   /*! Default constructor
    *  \param aName the name of the technology
    *  \param aYear the year
    */
   WindTechnology(
      const std::string& aName = std::string(),
      const int          aYear = -1 );
   /*! Copy constructor
    *  \param other the instance to copy
    */
   WindTechnology(const WindTechnology& other);

   //! Destructor
   virtual ~WindTechnology(void);

   // Documentation is inherited
   virtual void calcCost(
      const std::string& aRegionName,
      const std::string& aSectorName,
      const int          aPeriod );

   // Documentation is inherited
   virtual WindTechnology* clone( void ) const;

   // Documentation is inherited
   virtual void completeInit(
      const std::string&              aRegionName,
      const std::string&              aSectorName,
      const std::string&              aSubsectorName,
      DependencyFinder*               aDepFinder,
      const IInfo*                    aSubsectorIInfo,
      ILandAllocator*                 aLandAllocator );

   // Documentation is inherited
   static const std::string& getXMLNameStatic1D( void );

   // Documentation is inherited
   virtual void initCalc(
      const std::string& aRegionName,
      const std::string& aSectorName,
      const IInfo*       aSubsectorIInfo,
      const Demographic* aDemographics,
      PreviousPeriodInfo& aPrevPeriodInfo,
      const int          aPeriod );

    virtual double getCalibrationOutput( const bool aHasRequiredInput,
                                         const std::string& aRequiredInput, 
                                         const int aPeriod ) const;

protected :

   //! XML tag name indices
   enum
   {
      AIR_DENSITY_KEY,
      AVERAGE_WIND_SPEED_KEY,
      CAPITAL_COST_KEY,
      CUTOUT_SPEED_KEY,
      FCR_KEY,
      GRID_CONNECTION_COST_KEY,
      OM_KEY,
      REFERENCE_HEIGHT_KEY,
      ROTOR_DIAMETER_KEY,
      TURBINE_DENSITY_KEY,
      TURBINE_DERATING_KEY,
      TURBINE_HUB_HEIGHT_KEY,
      TURBINE_RATING_KEY,
      WIND_FARM_LOSS_KEY,
      WIND_VELOCITY_EXPONENT_KEY,

      //! Number of keys
      NUM_KEYS
   };

    virtual double getResourceToEnergyRatio( const std::string& aRegionName,
                                             const std::string& aSectorName,
                                             const int aPeriod );

   static std::string sXMLTagNames[NUM_KEYS];

   /*! Compute the ideal turbine output
    *  \param aAveWindSpeed the average Wind Speed
    *  \param aDiameter the turbine blade diameter (in meters)
    *  \param aAirDensity the average Air density (in g/m^3)
    */
   static double calcIdealTurbineOutput(
      double aAveWindSpeed,
       double aDiameter,
       double aAirDensity );

   /*! Compute the realized turbine output
    *  \param apInfo pointer to the market info
    */
   virtual double calcRealizedTurbineOutput( const IInfo* apInfo ) const;

   /*! Calculate the resource area in km^2
    *  \param aRegionName the region name
    *  \param aSectorName the sector name
    *  \param aVariableDemand the variable demand
    *  \param aPeriod the period
    *  \return the resource area
    */
   virtual double calcResourceArea(
      const std::string& aRegionName,
      const std::string& aSectorName,
      double             aVariableDemand,
      const int          aPeriod );

   /*! Compute the capture coefficient for a turbine with a finite power rating.
    *  \param aAveWindSpeed the average Wind Speed
    *  \param aRating the turbine Rating (in MW)
    *  \param aDiameter the turbine Blade Diameter (in meters)
    *  \param aAirDensity the average Air density (in g/m^3)
    *  \param aCutoutSpeed the cut-out Speed (m/s)
    */
   static double calcTurbineCoefficient(
      double aAveWindSpeed,
       double aRating,
       double aDiameter,
       double aAirDensity,
       double aCutoutSpeed );

   double computeWindPowerVariance(
      double aAveWindSpeed,
      double aRating,
      double aDiameter,
      double aAirDensity,
      double aCutoutSpeed ) const;

   // Documentation is inherited
   virtual const std::string& getXMLName1D( void ) const;

   // Documentation is inherited
   virtual void toDebugXMLDerived(
      const int     period,
      std::ostream& out,
      Tabs*         tabs ) const;

   // Documentation is inherited
   virtual void toInputXMLDerived(
      std::ostream& out,
      Tabs*         tabs ) const;

   // Documentation is inherited
   virtual bool XMLDerivedClassParse(
      const std::string&      nodeName,
      const xercesc::DOMNode* curr );

private :

   static const double kWhrtoGJ;

   //! Capital cost [2004 $/kW] based on rated capacity
   double mCapitalCost;

   //! Connection cost
   double mConnectCost;

   //! Turbine cutout speed [m/s]
   double mCutOutSpeed;

   //! fixed charge rate [unitless]
   double mFCR;

   //! Generation cost
   double mGenerationCost;

   //! Grid connection cost [$/km/MW capacity] based on rated capacity
   double mGridConnectionCost;

   //! [2004 $/kW/Yr] based on rated capacity
   double mOM;

   //! Realized turbine output
   double mRealizedTurbineOutput;

   //! Turbine rotor diameter [m]
   double mRotorDiameter;

   //! Land per turbine [MW/km2]
   double mTurbineDensity;

   //! Percentage of loss from ideal operation [unitless]
   double mTurbineDerating;

   //! Turbine hub height [m]
   double mTurbineHubHeight;

   //! Turbine output rating [unitless]
   double mTurbineRating;

   //! Wind capacity factor
   double mWindCapacityFactor;

   //! Percentage of loss from turbine to grid [unitless]
   double mWindFarmLoss;

   //! Wind Power Variance
   mutable double mWindPowerVariance;

   //! The maximum energy loss fraction
   double mMaxLoss;

   //! Exponent for energy loss as a function of penetration
   double mEfficiencyLossExponent;

   virtual const std::string& getTechCostName( ) const;
};

#endif   // __WIND_TECHNOLOGY_H

// end of wind_technology.h ************************************************

