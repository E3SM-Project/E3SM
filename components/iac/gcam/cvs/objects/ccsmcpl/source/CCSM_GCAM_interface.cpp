/*! 
* \file CCSM_GCAM_interface.cpp
* \ingroup Objects
* \brief CCSM gcam driver source file.
* \author John Truesdale
*/


#include "util/base/include/auto_file.h"
#include "../../ccsmcpl/include/CCSM_GCAM_interface.h"
#include "containers/include/world.h"

  ofstream outFile;

  LoggerFactoryWrapper loggerFactoryWrapper;

  // Pointer for a scenario
  Scenario* scenario; // model scenario info
  
  std::auto_ptr<CcsmRunner> runner;

  vector<string> CCSM_GCAM_interface::regionName;
  vector<string> CCSM_GCAM_interface::landType;
  vector<string> CCSM_GCAM_interface::cropName;

/*! \brief Constructor
* \details This is the constructor for the CCSM_driver class.
*/


CCSM_GCAM_interface::CCSM_GCAM_interface(){
}

//! Destructor. 
CCSM_GCAM_interface::~CCSM_GCAM_interface(){
}

/*! \brief Initializer for GCAM.
* \details
*  Initialize gcam log files and read in configuration
*  and base model information.  Create and setup scenario
*/

void CCSM_GCAM_interface::initGCAM(void)
{

  // identify default file names for control input and logging controls
  string configurationArg = "configuration.xml";
  string loggerFactoryArg = "log_conf.xml";

  // Add OS dependent prefixes to the arguments.
  const string configurationFileName = configurationArg;
  const string loggerFileName = loggerFactoryArg;

  // Initialize the LoggerFactory
  bool success = XMLHelper<void>::parseXML( loggerFileName, &loggerFactoryWrapper );

  // Get the main log file.
  ILogger& mainLog = ILogger::getLogger( "main_log" );
  mainLog.setLevel( ILogger::WARNING );

 // print disclaimer
  mainLog << "LEGAL NOTICE" << endl;
  mainLog << "This computer software was prepared by Battelle Memorial Institute," << endl;
  mainLog << "hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830" << endl;
  mainLog << "with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE" << endl;
  mainLog << "CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY" << endl;
  mainLog << "LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this" << endl;
  mainLog << "sentence must appear on any copies of this computer software." << endl;

  // print export control notice
  mainLog << "EXPORT CONTROL" << endl;
  mainLog << "User agrees that the Software will not be shipped, transferred or" << endl;
  mainLog << "exported into any country or used in any manner prohibited by the" << endl;
  mainLog << "United States Export Administration Act or any other applicable" << endl;
  mainLog << "export laws, restrictions or regulations (collectively the 'Export Laws')." << endl;
  mainLog << "Export of the Software may require some form of license or other" << endl;
  mainLog << "authority from the U.S. Government, and failure to obtain such" << endl;
  mainLog << "export control license may result in criminal liability under" << endl;
  mainLog << "U.S. laws. In addition, if the Software is identified as export controlled" << endl;
  mainLog << "items under the Export Laws, User represents and warrants that User" << endl;
  mainLog << "is not a citizen, or otherwise located within, an embargoed nation" << endl;
  mainLog << "(including without limitation Iran, Syria, Sudan, Cuba, and North Korea)" << endl;
  mainLog << "    and that User is not otherwise prohibited" << endl;
  mainLog << "under the Export Laws from receiving the Software." << endl;
  mainLog << "" << endl;
  mainLog << "Copyright 2011 Battelle Memorial Institute.  All Rights Reserved." << endl;
  mainLog << "Distributed as open-source under the terms of the Educational Community " << endl;
  mainLog << "License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php" << endl;
  mainLog << "" << endl;
  mainLog << "For further details, see: http://www.globalchange.umd.edu/models/gcam/" << endl;
  
  mainLog << "Running GCAM model code base version " << __ObjECTS_VER__ << " revision "
	  << __REVISION_NUMBER__ << endl << endl;
  
  // Parse configuration file.
  mainLog.setLevel( ILogger::NOTICE );
  mainLog << "Parsing input files..." << endl;
  Configuration* conf = Configuration::getInstance();
  success = XMLHelper<void>::parseXML( configurationFileName, conf );
  // Check if parsing succeeded. Non-zero return codes from main indicate


  // Initialize the timer.  Create an object of the Timer class.
  Timer timer;
  timer.start();

  // ******  packup gcamo  hard coded regions and landtypes for now
  regionName.push_back("USA");
  regionName.push_back("Canada");
  regionName.push_back("Western Europe");
  regionName.push_back("Japan");
  regionName.push_back("Australia_NZ");
  regionName.push_back("Former Soviet Union");
  regionName.push_back("China");
  regionName.push_back("Middle East");
  regionName.push_back("Africa");
  regionName.push_back("Latin America");
  regionName.push_back("Southeast Asia");
  regionName.push_back("Eastern Europe");
  regionName.push_back("Korea");
  regionName.push_back("India");
    
    
  landType.push_back("Build-up");
  landType.push_back("Cropland");
  landType.push_back("Pasture");
  landType.push_back("ForestHarvested");
  landType.push_back("Forest");
  landType.push_back("Other Land");
    
  cropName.push_back("biomass");
  cropName.push_back("Corn");
  cropName.push_back("eucalyptus");
  cropName.push_back("FiberCrop");
  cropName.push_back("FodderGrass");
  cropName.push_back("FodderHerb");
  cropName.push_back("Forest");
  cropName.push_back("Grassland");
  cropName.push_back("Jatropha");
  cropName.push_back("miscanthus");
  cropName.push_back("MiscCrop");
  cropName.push_back("OilCrop");
  cropName.push_back("OtherArableLand");
  cropName.push_back("OtherGrain");
  cropName.push_back("PalmFruit");
  cropName.push_back("Pasture");
  cropName.push_back("Rice");
  cropName.push_back("RockIceDesert");
  cropName.push_back("Root_Tuber");
  cropName.push_back("Shrubland");
  cropName.push_back("SugarCrop");
  cropName.push_back("Tundra");
  cropName.push_back("UnmanagedForest");
  cropName.push_back("UnmanagedPasture");
  cropName.push_back("UrbanLand");
  cropName.push_back("Wheat");
  cropName.push_back("willow");

  // Create an empty exclusion list so that any type of IScenarioRunner can be
  // created.

  //  create the scenario runner
  runner = auto_ptr<CcsmRunner>( new CcsmRunner );

  // Setup the scenario.
  success = runner->setupScenarios( timer );

  const Modeltime* modeltime = runner->getInternalScenario()->getModeltime();

  gcamStartYear = modeltime->getStartYear();
  gcamEndYear = modeltime->getEndYear();

  timer.stop();

}
void CCSM_GCAM_interface::setDensityGCAM(int *yyyymmdd, int *tod, double *gcami, int *gcami_fdim1_nflds, int *gcami_fdim2_datasize)
{

  bool success = false;
  ILogger& mainLog = ILogger::getLogger( "main_log" );
  mainLog.setLevel( ILogger::NOTICE );
  
  if (*yyyymmdd/10000 < gcamStartYear ) {
    mainLog << "returning from setDensityGCAM: current date: " <<*yyyymmdd<<" is before GCAM starting date: " << gcamStartYear << endl;
    return;
  }
  if (*yyyymmdd/10000 > gcamEndYear ) {
    mainLog << "returning from setDensityGCAM: current date: " <<*yyyymmdd<<" is after GCAM ending date: " << gcamEndYear << endl;
    return;
  }
  
  const Modeltime* modeltime = runner->getInternalScenario()->getModeltime();
  int curryear = *yyyymmdd/10000;
  int finalCalibrationYear=modeltime->getper_to_yr(modeltime->getFinalCalibrationPeriod());
  int period = modeltime->getyr_to_per( curryear );
  int modelyear= modeltime->getper_to_yr(period);

  mainLog << "curryear is " << curryear << endl;
  mainLog << "period is " << period << endl;
  mainLog << "modelyear is " << modelyear << endl;
  mainLog << "finalCalibrationYear is " << finalCalibrationYear << endl;

  // if we are past final calibration year then push carbon densities for this year
  if (curryear > finalCalibrationYear) {
    mainLog << "GCAM is setting Carbon data for year" << curryear << endl;
    
    //unpack gcami and set above and below ground carbon fields
    SetCarbonDensity& setDensities = runner->getInternalScenario()->getWorld()->getAdjustCarbonDensityVisitor();

    int gcami_nflds=2;
    int nreg = 14;
    int naez = 18;
    int ncrops = 27;
    int i=-1*gcami_nflds;
    for (int r=0;r<nreg;r++){
      for (int aez=1;aez<=naez;aez++){
	for (int ic=0;ic<ncrops;ic++){
	    mainLog << "inside for loop r="<<r<<" aez= "<<aez<<" crop "<<ic<<endl;
	  i+=gcami_nflds;
	  if (gcami[i] > 0.) {
	    setDensities.setCarbonDensityToPush( regionName[r], aez, cropName[ic], gcami[i], gcami[i+1],  curryear );
	    mainLog << "Setting above/below ground carbon for region "<<r<<" aez "<<aez<<" crop "<<ic<<" year "<<curryear<<"agc="<<gcami[i]<<"bgc="<<gcami[i+1] << endl;
	  }
	}
      }
    }
	    mainLog << "pushing carbon into the model" << endl;

    // push set carbon densities into the model
    runner->getInternalScenario()->accept( &setDensities, period );

    // for now just write a restart dataset every step
    //    runner->writeRestart( period, curryear );
    return;
  }
}
void CCSM_GCAM_interface::runGCAM(int *yyyymmdd, int *tod, double *gcami, int *gcami_fdim1_nflds, int *gcami_fdim2_datasize, double *gcamo,int *gcamo_fdim1_nflds,int *gcamo_fdim2_datasize, double *gcamoemis,int *gcamoemis_fdim1_nflds,int *gcamoemis_fdim2_datasize,int *yr1,int *yr2,int *sneakermode,int *write_rest)
{
  
  int curryear = *yyyymmdd/10000;
  bool success = false;
  ILogger& mainLog = ILogger::getLogger( "main_log" );
  mainLog.setLevel( ILogger::NOTICE );
  
  if (curryear < gcamStartYear ) {
    mainLog << "returning from runGCAM: current date: " <<*yyyymmdd<<" is before GCAM starting date: " << gcamStartYear << endl;
    return;
  }
  if (curryear > 2104 ) {
    mainLog << "returning from runGCAM: current date: " <<*yyyymmdd<<" is after GCAM ending date: " << gcamEndYear << endl;
    return;
  }
  
  // Extrapolate GCAM data for year 2100 fill gcamo gcamoemis
  if (curryear >= 2096 ) {
    const Modeltime* modeltime = runner->getInternalScenario()->getModeltime();
    int period2095 = modeltime->getyr_to_per( 2095 );
    int period2090 = modeltime->getyr_to_per( 2090 );
    
    mainLog << "curryear is " << curryear << endl;
    
    // calculate bounding years to pass to gcam2glm
    
    double tmp=0.;
    double tmp2090=0.;
    double tmp2095=0.;
    int in = -1;
    int gcamoyear[2];
    if (curryear< 2100) {
      gcamoyear[0]=2090;
      gcamoyear[1]=2095;
    }else{
      mainLog << "PAST GCAM ending year of "<< gcamEndYear << endl;
      mainLog << "Extrapolating GCAM data to 2100 using 2090 and 2095" << endl;
      gcamoyear[0]=2095;
      gcamoyear[1]=2100;
    }
    GetGLMData glmData2090;
    runner->getInternalScenario()->accept(&glmData2090, period2090 );
    GetGLMData glmData2095;
    runner->getInternalScenario()->accept(&glmData2095, period2095 );

    for (int yr = 0;yr<2;yr++) {
      mainLog << "gcamoyear[" << yr <<"]=" << gcamoyear[yr] << endl;
      for (int i = 0;i<14;i++) {
	for (int aez = 1;aez<19;aez++) {
	  for (int n = 0;n<*gcamo_fdim1_nflds;n++) {
	    in++;
	    if (landType[n] != "ForestHarvested") { 
	      if (gcamoyear[yr]==2090) {
		tmp =  glmData2090.getLandCover( regionName[i], aez, landType[n], gcamoyear[yr] );	
	      }
	      if (gcamoyear[yr]==2095) {
		tmp =  glmData2095.getLandCover( regionName[i], aez, landType[n], gcamoyear[yr] );	
	      }
	      if (gcamoyear[yr]==2100) {
		tmp2090 =  glmData2090.getLandCover( regionName[i], aez, landType[n], 2090 );	
		tmp2095 =  glmData2095.getLandCover( regionName[i], aez, landType[n], 2095 );	
		tmp=(tmp2095-tmp2090)+tmp2095;
   	        mainLog << "values for 2090 and 2095 are" << tmp2090 <<"/"<< tmp2095 << endl;
	      }
	      gcamo[in] = tmp<0.?0.:tmp;
	      mainLog << "packing gcamo" << regionName[i] << ":" << aez << ":" << landType[n] <<" gcamo[" << in << "]=" << gcamo[in] << endl;
	      
	    }else{
	      if (gcamoyear[yr]==2090) {
		tmp = glmData2090.getProductionInCarbon( regionName[i], aez, "Forest", gcamoyear[yr] );	 
	      }
	      if (gcamoyear[yr]==2095) {
		tmp = glmData2095.getProductionInCarbon( regionName[i], aez, "Forest", gcamoyear[yr] );	 
	      }
	      if (gcamoyear[yr]==2100) {
		tmp2090 = glmData2090.getProductionInCarbon( regionName[i], aez, "Forest", 2090 );	 
		tmp2095 = glmData2095.getProductionInCarbon( regionName[i], aez, "Forest", 2095 );	 
		tmp=(tmp2095-tmp2090)+tmp2095;
   	        mainLog << "values for 2090 and 2095 are" << tmp2090 <<"/"<< tmp2095 << endl;
	      }
	      gcamo[in] = tmp<0.?0.:tmp;
	      mainLog << "packing gcamo" << regionName[i] <<":" << aez << ":" << landType[n] <<" gcamo[" << in << "]=" << gcamo[in] << endl;
	    }  //landtype
	  }  //for n
	}  //for aez
      }  //for i
    }  // for yr
    *yr1=gcamoyear[0];
    *yr2=gcamoyear[1];
    
    // packup gcamoemiss array with just co2 emissions for now will add others later
    in=-1;
    const int numGases = 1;
    const string gases[] = { "CO2" };
    //    const int numGases = 8;
    //    const string gases[] = { "BC", "CH4", "CO", "NH3", "NMVOC", "NOx", "OC", "SO2" };
    const int numCategories = 12;
    const string categories[] = { "AGR", "AIR", "AWB", "DOM", "ENE", "IND", "LCF", "SAV", "SHIP", "SLV", "TRA", "WST" };
    
    RCPEmissionsVisitor emissionsVisitor2090;
    RCPEmissionsVisitor emissionsVisitor2095;
    //    const Modeltime* modeltime = scenario->getModeltime();
    //    const int year = modeltime->getper_to_yr( aPeriod );
    runner->getInternalScenario()->accept( &emissionsVisitor2090, period2090 );
    runner->getInternalScenario()->accept( &emissionsVisitor2095, period2095 );
    std::vector<Region*> regions=runner->getInternalScenario()->getWorld()->getRegions();
    for( RegionIterator region = regions.begin(); region != regions.end(); ++region ){
      for( int catI = 0; catI < numCategories; ++catI ) {
	for( int gasI = 0; gasI < numGases; ++gasI ) {
	  in++;
	  if (gcamoyear[1]==2100) {
	    tmp2090 = emissionsVisitor2090.getEmissions( (*region)->getName(), categories[ catI ], gases[ gasI ], 2090 );
	    tmp2095 = emissionsVisitor2095.getEmissions( (*region)->getName(), categories[ catI ], gases[ gasI ], 2095 );
	    tmp=(tmp2095-tmp2090)+tmp2095;
	    gcamoemis[in] = tmp;
	    mainLog << "values for 2090 and 2095 are" << tmp2090 <<"/"<< tmp2095 << endl;
	    mainLog << "packing gcamoemis" << (*region)->getName() << ":" << gases[ gasI ] << ":" << categories[ catI ] << gcamoyear[1] <<" gcamoemis[" << in << "]=" << gcamoemis[in] << endl;
	  }else{
	    tmp = emissionsVisitor2095.getEmissions( (*region)->getName(), categories[ catI ], gases[ gasI ], gcamoyear[1] );
	  }
	}
      }
    }
    
  }else{   //Run GCAM and fill gcamo and gcamoemis

    const Modeltime* modeltime = runner->getInternalScenario()->getModeltime();
    int finalCalibrationYear=modeltime->getper_to_yr(modeltime->getFinalCalibrationPeriod());
    int period = modeltime->getyr_to_per( curryear );
    int modelyear= modeltime->getper_to_yr(period);
    
    mainLog << "curryear is " << curryear << endl;
    mainLog << "period is " << period << endl;
    mainLog << "modelyear is " << modelyear << endl;
    mainLog << "finalCalibrationYear is " << finalCalibrationYear << endl;
    
    
    if(modeltime->isModelYear( curryear)) {
      
      Timer timer;
      
      // TODO: is this necessary, it will be the same as currYear
      mainLog << "Running GCAM for year" << modelyear << endl;
      mainLog << "calculating period=" << period << endl;
      
      mainLog.precision(20);
      
      // Initialize the timer.  Create an object of the Timer class.
      timer.start();
      
      success = runner->runScenarios( period, true, timer );
      
      timer.stop();
      
      // write restarts if needed.
      if(write_rest) {
	mainLog << "write_rest: " << *write_rest << endl;
	runner->writeRestart( period, curryear );
      }
    }
    
    //only worry about packing gcamo past final calibration year
    
    if (curryear>=finalCalibrationYear) {
      
      // calculate bounding years to pass to gcam2glm
      int bndyr1,bndyr2;
      mainLog << "sneakermode: " << *sneakermode << endl;
      if (*sneakermode) {
	bndyr1 =  period<4?gcamStartYear:modeltime->getper_to_yr(period-3);
	bndyr2 =  modeltime->getper_to_yr(period);
      }else{
	if(modeltime->isModelYear( curryear)) {
	  bndyr1 =  period==0?gcamStartYear:modeltime->getper_to_yr(period-1);
	  bndyr2 =  modeltime->getper_to_yr(period);
	}else{
	  bndyr1 =  period==0?gcamStartYear:modeltime->getper_to_yr(period-2);
	  bndyr2 =  modeltime->getper_to_yr(period-1);
	}
      }
      
      //    mainLog << "prev year: " << prevyear << endl;
      
      // //     // ****** This is an example of how to query GCAM data for GLM ******
      // //     // A Period value of -1 indicates gather data for all model periods.
      
      // mainLog << "calling runner accept with period: " << period << endl;
      // GetGLMData glmData;
      // runner->getInternalScenario()->accept(&glmData, period );
      
      // // Region and land categories are GLM categories, the years are GCAM years
      // // for example 1990 to 2095 in 15 curryear time steps.
      // int aez=7 ;
      // mainLog << "Year: " << modelyear << endl;
      // mainLog << "USA Cropland: " << glmData.getLandCover( "USA", aez, "Cropland", modelyear ) << endl;
      // mainLog << "USA Forest: " << glmData.getLandCover( "USA", aez, "Forest", modelyear ) << endl;
      // mainLog << "USA ForestHarvested: " << glmData.getLandCover( "USA", aez, "ForestHarvested", modelyear ) << endl;
      // mainLog << "USA Grassland: " << glmData.getLandCover( "USA", aez, "Grassland", modelyear ) << endl;
      // mainLog << "USA Pasture: " << glmData.getLandCover( "USA", aez, "Pasture", modelyear ) << endl;
      // mainLog << "USA Build-up: " << glmData.getLandCover( "USA", aez, "Build-up", modelyear ) << endl;
      // mainLog << "USA Other Land: " << glmData.getLandCover( "USA", aez, "Other Land", modelyear ) << endl;
      // mainLog << "USA Forest production in tC: " << glmData.getProductionInCarbon( "USA", aez, "Forest", modelyear ) << endl;
      // pair<double, double> carbonDensity = glmData.getCarbonDensity( "USA", aez, "Corn", modelyear );
      // mainLog << "USA Corn carbon densities: " << carbonDensity.first <<  ", " << carbonDensity.second << endl;
      
      double tmp=0.;
      int in = -1;
      int glmdataPeriod;
      int gcamoyear[2];
      for (int yr = 0;yr<2;yr++) {
	GetGLMData glmData;
	gcamoyear[yr] = yr==0?bndyr1:bndyr2;
	glmdataPeriod = modeltime->getyr_to_per( gcamoyear[yr] );
	runner->getInternalScenario()->accept(&glmData, glmdataPeriod );
	mainLog << "gcamoyear[" << yr <<"]=" << gcamoyear[yr] << endl;
	mainLog << "gcamoyear period" << glmdataPeriod << endl;
	for (int i = 0;i<14;i++) {
	  for (int aez = 1;aez<19;aez++) {
	    for (int n = 0;n<*gcamo_fdim1_nflds;n++) {
	      in++;
	      if (landType[n] != "ForestHarvested") { 
		tmp =  glmData.getLandCover( regionName[i], aez, landType[n], gcamoyear[yr] );	
		gcamo[in] = tmp<0.?0.:tmp;
		mainLog << "packing gcamo" << regionName[i] << ":" << aez << ":" << landType[n] <<" gcamo[" << in << "]=" << gcamo[in] << endl;
	      }else{
		tmp = glmData.getProductionInCarbon( regionName[i], aez, "Forest", gcamoyear[yr] );	 
		gcamo[in] = tmp<0.?0.:tmp;
		mainLog << "packing gcamo" << regionName[i] <<":" << aez << ":" << landType[n] <<" gcamo[" << in << "]=" << gcamo[in] << endl;
	      } //fi landtype
	    } //for fieldnum
	  } //for aez
	} //for region
      } //for year 0,1
      *yr1=gcamoyear[0];
      *yr2=gcamoyear[1];

      // packup gcamoemiss array with just co2 emissions for now will add others later
      in=-1;      
      const int numGases = 1;
      const string gases[] = { "CO2" };
      //    const int numGases = 8;
      //    const string gases[] = { "BC", "CH4", "CO", "NH3", "NMVOC", "NOx", "OC", "SO2" };
      const int numCategories = 12;
      const string categories[] = { "AGR", "AIR", "AWB", "DOM", "ENE", "IND", "LCF", "SAV", "SHIP", "SLV", "TRA", "WST" };
      
      RCPEmissionsVisitor emissionsVisitor;
      //    const Modeltime* modeltime = scenario->getModeltime();
      //    const int year = modeltime->getper_to_yr( aPeriod );
      runner->getInternalScenario()->accept( &emissionsVisitor, period );
      std::vector<Region*> regions=runner->getInternalScenario()->getWorld()->getRegions();
      for( RegionIterator region = regions.begin(); region != regions.end(); ++region ){
	for( int catI = 0; catI < numCategories; ++catI ) {
	  for( int gasI = 0; gasI < numGases; ++gasI ) {
	    in++;
	    tmp = emissionsVisitor.getEmissions( (*region)->getName(), categories[ catI ], gases[ gasI ], gcamoyear[1] );
	    //		rcpLog << ( region - regions.begin() + 1 ) << ", " << gases[ gasI ] << ", " << categories[ catI ] << ", "
	    //		       << year << ", "
	    //		       << tmp << endl;
	    gcamoemis[in] = tmp<0.?0.:tmp;
	    mainLog << "packing gcamoemis" << (*region)->getName() << ":" << gases[ gasI ] << ":" << categories[ catI ] << gcamoyear[1] <<" gcamoemis[" << in << "]=" << tmp << endl;
	  } //gasI
	} //catI
      } //region
    } //if final calibration year
  }  //end curryr==2100
}
void CCSM_GCAM_interface::finalizeGCAM()
{
  Timer timer;
  // Initialize the timer.  Create an object of the Timer class.
  timer.start();
  ILogger& mainLog = ILogger::getLogger( "main_log" );
  mainLog.setLevel( ILogger::NOTICE );
  mainLog << "calling finalize" << endl;
  runner->printOutput(timer);
  timer.stop();
}
