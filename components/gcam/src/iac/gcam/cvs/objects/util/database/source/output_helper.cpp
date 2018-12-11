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
 * \file output_helper.cpp
 * \ingroup Objects
 * \brief Contains the functions to write out single records to output file or
 *        table in the database. includes functions to open, create, and close
 *        database and tables.
 * \author Sonny Kim
 */

#include "util/base/include/definitions.h"
#define DUPLICATE_CHECKING 0
#if(__HAVE_DB__)
#include <afxdisp.h>
#include <dbdao.h>
#include <comdef.h>
#endif

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <ctime>
#include "containers/include/scenario.h" 
#include "containers/include/world.h"
#include "util/base/include/model_time.h"
#include "util/base/include/configuration.h"
#include "util/base/include/util.h"

#if( DUPLICATE_CHECKING )
#include "util/base/include/util.h"
#include <algorithm>
#endif

using namespace std;

extern time_t gGlobalTime;
extern ofstream outFile;
extern Scenario* scenario;

#if( DUPLICATE_CHECKING )
vector<string> gHashes;
#endif

#if(__HAVE_DB__)
// define global DB engine and database
CdbDBEngine dben;
CdbDatabase db;
CdbTableDef DBoutTD; // table definitions for creating tables in the database
CdbField tfield; // temporary field for creating fields in tables
CdbRecordset DBoutrst; // record set for writing results
const char *DBout = "DBout"; // name of the table for outputs compatible with dataviewer
#endif

/*! Output single records to file.
*
* Names of categories, subcategories, and variables are written
* as strings in the argument.  Values are also passed as arguments.
*
*/

void fileoutput3( string var1name,string var2name,string var3name,
              string var4name,string var5name,string uname,vector<double> dout) {
    outFile <<var1name<<","<<var2name<<","<<var3name<<","
            <<var4name<<","<<var5name<<","<<uname<<",";
    for (int i=0;i< static_cast<int>( dout.size() );i++) {
        outFile << dout[i]<<",";
    }
    outFile << endl;
}

/*! Output single records MiniCAM style to the database. 
*
* Names of categories, subcategories, variables, and units are written
* as strings passed in by the argument.  Values are also passed as arguments.
* \warning Variable id field has been taken out of the dbout database table
and this function no longer compensate for this field.  Ensure that the 
revised dbout table is used.
*/
void dboutput4(string var1name,string var2name,string var3name,string var4name,
               string uname,vector<double> dout)
{
#if( DUPLICATE_CHECKING )
    // Create a hash of the values.
    const string newHash = var1name + "," + var2name + "," + var3name + "," + var4name;
    // See if its already been added to the map.
    if( find( gHashes.begin(), gHashes.end(), newHash ) != gHashes.end() ){
        // Print an error since we found a duplicate.
        cout << "A value has already been added for: " << var1name << ", " << var2name
            << ", " << var3name << ", " << var4name << endl;
    }
    // Otherwise add it.
    else {
        gHashes.push_back( newHash );
    }
#endif // DUPLICATE_CHECKING

#if(!__HAVE_DB__)
    fileoutput3( var1name,var2name,var3name,var4name,"",uname,dout);
#else
    // MS Access can not handle blank strings. This is a programming error
    // and would ideally be handled with an assert. However, running in debug
    // mode until the database writes out can take a very long time and so
    // these errors would not be found often enough.
    if ( var1name.empty() || var2name.empty() || var3name.empty() || var4name.empty() || 
        uname.empty() ) {
            cout << "Variable name is blank! Check data and output routines." << endl;
            cout << "V1: " << var1name << " V2: " << var2name << " V3: " << var3name
                << " V4: " << var4name << " unit: " << uname << endl;
            return;
        }

        const map<string,int> regionMap = scenario->getWorld()->getOutputRegionMap();

        // only output regions already defined in regionMap
        map<string,int>::const_iterator iter = regionMap.find( var1name );
        if( iter != regionMap.end()) {
            int i=0,j;
            // COleVariant does not take int, must be short or long
            // now write to the database
            // MUST Update before adding another
            try {
                DBoutrst.AddNew(); // now the current record is this empty new one

                const long runID = util::createMinicamRunID( gGlobalTime );

                DBoutrst.SetField(0L, COleVariant(runID,VT_I4)); // run id
                DBoutrst.SetField(1L, COleVariant(short(iter->second))); // region id
                DBoutrst.SetField(2L, COleVariant(var2name.c_str(), VT_BSTRT)); // category
                DBoutrst.SetField(3L, COleVariant(var3name.c_str(), VT_BSTRT)); // subcategory
                DBoutrst.SetField(4L, COleVariant(var4name.c_str(), VT_BSTRT)); // variable
                DBoutrst.SetField(5L, COleVariant(uname.c_str(), VT_BSTRT)); // units
                // get scenario name and output to dbout
                string scenarioName = scenario->getName();
                DBoutrst.SetField(6L, COleVariant(scenarioName.c_str(), VT_BSTRT));
                for (i=0;i< static_cast<int>( dout.size() );i++) {
                    j = 7+i;
                    DBoutrst.SetField(j, COleVariant(dout[i]));
                }
                DBoutrst.Update(); // save and write the record        
            } 
            catch( _com_error e ){
                printf("Error:*************************************************\n");
                printf("Code = %08lx\n", e.Error());
                printf("Message = %s\n", e.ErrorMessage());
                printf("Source = %s\n", (LPCSTR) e.Source());
                printf("Description = %s\n", (LPCSTR) e.Description());
            }
        }
        // if region not found in regionMap, do nothing
#endif
}

#if(__HAVE_DB__)
//! Open connection to the database
void openDB(void) {
    Configuration* conf = Configuration::getInstance();
    string dbFile = conf->getFile( "dbFileName" );
    // Open a global Jet database in exclusive, read/write mode.
    try {
        db = dben.OpenDatabase( dbFile.c_str()); 
    }
    catch(...) {
        cout<<"Error opening database: "<< dbFile << endl; 
    }
}

//! Close connection to the database
void closeDB() {
    db.Close();
}

/*! \brief Create and open the main database output table.
* 
* Record set for dbout table is opened here for writing to by output_helper.cpp 
* The code to delete and create new dbout table is commented out.
* Results are continually appended to the dbout table.  This function no longer
* clears out and recreates the dbout table.
* \warning The output database has been greatly simplified and fewer tables exists in
the database.  Results are continually appended to the dbout table, which is now the
only table of results.  RunLabel and RegionInfo talbles exists in the database and 
are recreated in createMCvarid().
*/

void createDBout() {
    // open DBout table as a record set (DBoutrst) for writing
    DBoutrst = db.OpenRecordset(DBout,dbOpenDynaset);
    // **** End DBout table *****   
}

//! Create run label and region info tables that are necessary for the dataviewer.xls.
void createMCvarid() {
    string sqltemp; // temporary string for sql

    // *** Delete and Create DbRunLabels table ***
    const char *dbtrun = "DBRunLabels"; // database table names
    // first delete existing table
    try { db.TableDefs.Delete(dbtrun); } 
    catch (...) {cout<<"\nError deleting "<<dbtrun<<" table\n";}
    // Reusing RunLabelTD and tfield to add fields to the table
    CdbField tfield; // temporary field for creating fields in tables
    CdbTableDef RunLabelTD = db.CreateTableDef(dbtrun); // temporary RunLabel table definition
    // create run id field
    tfield = RunLabelTD.CreateField("RunID",dbLong);
    RunLabelTD.Fields.Append(tfield);
    // create run label field
    tfield = RunLabelTD.CreateField("RunLabel",dbText);
    RunLabelTD.Fields.Append(tfield);
    // create run comments field
    tfield = RunLabelTD.CreateField("RunComments",dbText);
    RunLabelTD.Fields.Append(tfield);
    // create the region info table
    try {db.TableDefs.Append(RunLabelTD);}
    catch(...) { cout<<"\nError appending "<<dbtrun<<" table to database\n";}

    // insert into the new run labels table the list of runs
    sqltemp = "INSERT INTO ";
    sqltemp += dbtrun;
    sqltemp += " SELECT DISTINCT RunID,RunLabel";
    sqltemp = sqltemp + " FROM " + DBout;
    try {
        db.Execute(sqltemp.c_str()); }
    catch(...) {
        cout<<"\nError executing sql: \""<<sqltemp<<"\"\n"; }
    // *** end DbRunLabels table ***

    // *** Delete and Create RegionInfo table ***
    const char *dbtregion = "RegionInfo"; // database table names
    // first delete existing table
    try { db.TableDefs.Delete(dbtregion); } 
    catch (...) {cout<<"\nError deleting "<<dbtregion<<" table\n";}
    // create new region info table definition
    CdbTableDef RegionInfoTD = db.CreateTableDef(dbtregion);
    // create region name field
    tfield = RegionInfoTD.CreateField("RegionLabel",dbText);
    RegionInfoTD.Fields.Append(tfield);
    // create region id field
    tfield = RegionInfoTD.CreateField("RegionID",dbLong);
    RegionInfoTD.Fields.Append(tfield);
    // create the region info table
    try {db.TableDefs.Append(RegionInfoTD);}
    catch(...) { cout<<"\nError appending "<<dbtregion<<" table to database\n";}

    // create a record set for the region info table
    CdbRecordset RegionInfoRst = db.OpenRecordset(dbtregion,dbOpenDynaset);
    typedef map<string,int>:: const_iterator CI;
    string regstr; // temporary string for region names
    map<string,int> regionMap = scenario->getWorld()->getOutputRegionMap();
    for (CI rmap=regionMap.begin(); rmap!=regionMap.end(); ++rmap) {
        RegionInfoRst.AddNew(); // now the current record is this empty new one
        regstr = rmap->first; // region name
        RegionInfoRst.SetField(_T("RegionLabel"), COleVariant(regstr.c_str(), VT_BSTRT));
        RegionInfoRst.SetField(_T("RegionID"), COleVariant(long(rmap->second), VT_I4));
        RegionInfoRst.Update(); // save and write the record
        // add a Global region that includes all regions to sum all regional outputs
        if(rmap->second != 0) {
            RegionInfoRst.AddNew(); // now the current record is this empty new one
            RegionInfoRst.SetField(_T("RegionLabel"), COleVariant("zGlobal", VT_BSTRT));
            RegionInfoRst.SetField(_T("RegionID"), COleVariant(long(rmap->second), VT_I4));
            RegionInfoRst.Update(); // save and write the record
        }
    }
    RegionInfoRst.Close();
    // *** end RegionInfo table ***
}
// if not working with database
#else
void closeDB(void) {
    // do nothing if not utilizing database
}
void openDB(void) {
    // do nothing if not utilizing database
}
void createDBout(void) {
    // do nothing if not utilizing database
}
void createMCvarid(void) {
    // do nothing if not utilizing database
}
#endif
