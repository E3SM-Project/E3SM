//===-- Test driver for OMEGA Metadata -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Metadata
///
/// This driver tests the Metadata capabilities for the OMEGA
/// model
//
//===-----------------------------------------------------------------------===/

#include <iostream>
#include <limits>

#include "DataTypes.h"
#include "Logging.h"
#include "MetaData.h"

using namespace OMEGA;

const I4 FILL_VALUE = 0;

int printResult(bool result, bool expected, std::string msg_pass,
                std::string msg_fail) {

   int RetVal = 0;

   if (result == expected) {
      std::cout << msg_pass << ": PASS" << std::endl;
   } else {
      RetVal = 1;
      std::cout << msg_fail << ": FAIL" << std::endl;
   }

   return RetVal;
}

int testMetaDim() {

   int RetVal = 0;

   const std::string DimName{"MyDim"};
   const I4 DimValue{1};

   // test has()
   RetVal += printResult(MetaDim::has(DimName), false,
                         "'" + DimName + "' is not created",
                         "'" + DimName + "' should not exist");

   // test create()
   auto Dim1 = MetaDim::create(DimName, DimValue);

   RetVal +=
       printResult(MetaDim::has(DimName), true, "'" + DimName + "' is created",
                   "'" + DimName + "' should exist");

   // test get()
   auto Dim2 = MetaDim::get(DimName);

   RetVal += printResult(Dim1 == Dim2, true, "get() returns correct instance.",
                         "get() returns incorrect instance.");

   // test getLength()
   I4 Length;
   Dim1->getLength(Length);

   RetVal += printResult(DimValue == Length, true,
                         "getLength() returns correct length.",
                         "getLength() returns incorrect length.");

   // test destroy()
   MetaDim::destroy(DimName);

   RetVal += printResult(MetaDim::has(DimName), false,
                         "'" + DimName + "' is destroyed correctly",
                         "'" + DimName + "' is not destroyed");

   return RetVal;
}

int testMetaData() {

   int RetVal = 0;

   const std::string ArrayName{"MyArray"};
   const std::string DimName{"MyDim"};
   const I4 DimValue{1};
   int ret;

   // test has()
   RetVal += printResult(MetaData::has(CodeMeta), false,
                         "'" + CodeMeta + "' is not created",
                         "'" + CodeMeta + "' should not exist");

   // test create() 1
   auto Data1 = MetaData::create(CodeMeta);

   RetVal += printResult(MetaData::has(CodeMeta), true,
                         "'" + CodeMeta + "' is created",
                         "'" + CodeMeta + "' should exist");

   // test create() 2
   auto Data2 = MetaData::create(SimMeta, {std::make_pair("Meta1", 1),
                                           std::make_pair("Meta2", 2),
                                           std::make_pair("Meta3", 3)});

   RetVal +=
       printResult(MetaData::has(SimMeta), true, "'" + SimMeta + "' is created",
                   "'" + SimMeta + "' should exist");

   std::map<std::string, std::any> *VarMeta = Data2->getAllEntries();
   std::string MetaName;
   int MetaVal, Count = 1;

   // use std::map functionality to loop through all metadata
   for (const auto &Pair : *VarMeta) {
      MetaName = Pair.first; // retrieve name part of meta pair
      MetaVal  = std::any_cast<int>(Pair.second); // retrieve value part

      RetVal += printResult(MetaName == "Meta" + std::to_string(MetaVal), true,
                            "'" + SimMeta + "' has correct MetaName",
                            "'" + SimMeta + "' has wrong MetaName");

      RetVal += printResult(MetaVal == Count, true,
                            "'" + SimMeta + "' has correct MetaVal",
                            "'" + SimMeta + "' has wrong MetaVal");

      // do whatever needed with metadata
      // metaVal can be cast into the appropriate type using
      //   std::any_cast<type>(metaVal)
      Count++;
   }

   // test create() 3
   std::vector<std::shared_ptr<MetaDim>> Dimensions;
   Dimensions.push_back(MetaDim::create(DimName, DimValue));

   auto Data3 = ArrayMetaData::create(
       ArrayName,
       "Description",                   /// long Name or description
       "Units",                         /// units
       "StdName",                       /// CF standard Name
       std::numeric_limits<int>::min(), /// min valid value
       std::numeric_limits<int>::max(), /// max valid value
       FILL_VALUE,                      /// scalar used for undefined entries
       1,                               /// number of dimensions
       Dimensions                       /// dim pointers
   );

   RetVal += printResult(MetaData::has(ArrayName), true,
                         "'" + ArrayName + "' is created",
                         "'" + ArrayName + "' should exist");

   // test get()
   auto Data4 = MetaData::get(ArrayName);

   RetVal +=
       printResult(Data3 == Data4, true, "get() returns correct instance.",
                   "get() returns incorrect instance.");

   // test hasEntry()
   RetVal += printResult(Data4->hasEntry("FillValue"), true,
                         "'" + ArrayName + "' has a fill value.",
                         "'" + ArrayName + "' does not have a fill value");

   // test getEntry()
   I4 FillValue;
   ret = Data3->getEntry("FillValue", FillValue);

   RetVal += printResult(ret == 0, true, "getEntry() returns zero",
                         "getEntry() returns non-zero");

   RetVal += printResult(FILL_VALUE == FillValue, true,
                         "'" + ArrayName + "' has a correct fill value",
                         "'" + ArrayName + "' has an incorrect fill value");

   // test addEntry()
   const R8 NewValue = 2.0;
   ret               = Data3->addEntry("NewMeta", NewValue);

   RetVal += printResult(ret == 0, true, "addEntry() returns zero",
                         "addEntry() returns non-zero");

   R8 R8Value;
   ret = Data3->getEntry("NewMeta", R8Value);

   RetVal += printResult(NewValue == R8Value, true,
                         "getEntry() returns correct value",
                         "getEntry() returns incorrect value");

   // test removeEntry()
   ret = Data3->removeEntry("NewMeta");

   RetVal += printResult(ret == 0, true, "removeEntry() returns zero",
                         "removeEntry() returns non-zero");

   RetVal += printResult(Data3->hasEntry("NewMeta"), false,
                         "'NewMeta' is removed correctly",
                         "'NewMeta' is not removed");

   // test destroy()
   ret = MetaData::destroy(SimMeta);

   RetVal += printResult(ret == 0, true, "destroy() returns zero",
                         "destroy() returns non-zero");

   RetVal += printResult(MetaData::has(SimMeta), false,
                         "'" + SimMeta + "' is correclty removed",
                         "'" + SimMeta + "' is not removed.");

   return RetVal;
}

int testMetaGroup() {

   int RetVal = 0;

   const std::string GroupName{"MyGroup"};
   const std::string FieldName{"MyField"};
   const std::string DimName{"MyDim"};
   const I4 DimValue{1};
   int ret;

   // test has()
   RetVal += printResult(MetaGroup::has(GroupName), false,
                         "'" + GroupName + "' is not created",
                         "'" + GroupName + "' should not exist");

   // test create()
   auto Group1 = MetaGroup::create(GroupName);

   RetVal += printResult(MetaGroup::has(GroupName), true,
                         "'" + GroupName + "' is created",
                         "'" + GroupName + "' should exist");

   // test get()
   auto Group2 = MetaGroup::get(GroupName);

   RetVal +=
       printResult(Group1 == Group2, true, "get() returns correct instance.",
                   "get() returns incorrect instance.");

   // test hasField()
   RetVal += printResult(Group1->hasField(FieldName), false,
                         "'" + FieldName + "' is not in a group",
                         "'" + FieldName + "' is in a group");

   // test addField()
   auto Data1 = MetaData::create(FieldName);

   ret = Group1->addField(FieldName);

   RetVal += printResult(ret == 0, true, "addField() returns zero.",
                         "addField() returns non-zero.");

   RetVal += printResult(Group1->hasField(FieldName), true,
                         "'" + FieldName + "' is in a group",
                         "'" + FieldName + "' is not in a group");

   std::map<std::string, std::shared_ptr<MetaData>> *Fields;
   Fields = Group1->getAllFields();

   std::string FieldNamePart;
   std::shared_ptr<MetaData> FieldPart;

   // use std::map functionality to loop through all fields
   for (const auto &Pair : *Fields) {
      FieldNamePart = Pair.first;  // retrieve field name part
      FieldPart     = Pair.second; // retrieve field part

      RetVal += printResult(FieldNamePart == FieldName, true,
                            "Correct FieldName is returned",
                            "Incorrect FieldName is returned");
   }

   // test getField()
   auto Data2 = Group1->getField(FieldName);

   RetVal +=
       printResult(Data1 == Data2, true, "getField() returns correct instance.",
                   "getField() returns incorrect instance.");

   // test removeField()
   ret = Group1->removeField(FieldName);

   RetVal += printResult(ret == 0, true, "removeField() returns zero.",
                         "removeField() returns non-zero.");

   RetVal += printResult(Group1->hasField(FieldName), false,
                         "'" + FieldName + "' is not in a group",
                         "'" + FieldName + "' is in a group");

   // test destroy()
   MetaGroup::destroy(GroupName);

   RetVal += printResult(MetaGroup::has(GroupName), false,
                         "'" + GroupName + "' is destroyed correctly",
                         "'" + GroupName + "' is not destroyed");

   return RetVal;
}

std::vector<std::shared_ptr<MetaDim>> initMetaDim(const std::string DimName,
                                                  const I4 DimValue) {

   std::vector<std::shared_ptr<MetaDim>> Dimensions;
   Dimensions.push_back(MetaDim::create(DimName, DimValue));

   return Dimensions;
}

void initMetaData(const std::string FieldName,
                  const std::vector<std::shared_ptr<MetaDim>> Dimensions) {

   auto Data = ArrayMetaData::create(
       FieldName,
       "Description",                   /// long Name or description
       "Units",                         /// units
       "StdName",                       /// CF standard Name
       std::numeric_limits<int>::min(), /// min valid value
       std::numeric_limits<int>::max(), /// max valid value
       FILL_VALUE,                      /// scalar used for undefined entries
       1,                               /// number of dimensions
       Dimensions                       /// dim pointers
   );
}

void initMetaGroup(const std::string GroupName) {

   auto Group1 = MetaGroup::create(GroupName);
}

int testMetaInit() {

   int RetVal = 0;

   const std::string GroupName{"MyInitGroup"};
   const std::string FieldName{"MyInitField"};
   const std::string DimName{"MyInitDim"};
   const I4 DimValue{1};
   int ret;

   auto dimensions = initMetaDim(DimName, DimValue);

   initMetaData(FieldName, dimensions);

   initMetaGroup(GroupName);

   RetVal += printResult(MetaGroup::has(GroupName), true,
                         "'" + GroupName + "' is created",
                         "'" + GroupName + "' should exist");

   // test get()
   auto Group1 = MetaGroup::get(GroupName);

   // test hasField()
   RetVal += printResult(Group1->hasField(FieldName), false,
                         "'" + FieldName + "' is not in a group",
                         "'" + FieldName + "' is in a group");

   ret = Group1->addField(FieldName);

   RetVal += printResult(ret == 0, true, "addField() returns zero.",
                         "addField() returns non-zero.");

   RetVal += printResult(Group1->hasField(FieldName), true,
                         "'" + FieldName + "' is in a group",
                         "'" + FieldName + "' is not in a group");

   // test getField()
   auto Data1 = MetaData::get(FieldName);
   auto Data2 = Group1->getField(FieldName);

   RetVal +=
       printResult(Data1 == Data2, true, "getField() returns correct instance.",
                   "getField() returns incorrect instance.");

   // test removeField()
   ret = Group1->removeField(FieldName);

   RetVal += printResult(ret == 0, true, "removeField() returns zero.",
                         "removeField() returns non-zero.");

   RetVal += printResult(Group1->hasField(FieldName), false,
                         "'" + FieldName + "' is not in a group",
                         "'" + FieldName + "' is in a group");

   // test MetaGroup::destroy()
   MetaGroup::destroy(GroupName);

   RetVal += printResult(MetaGroup::has(GroupName), false,
                         "'" + GroupName + "' is destroyed correctly",
                         "'" + GroupName + "' is not destroyed");

   // test MetaData::destroy()
   MetaData::destroy(FieldName);

   RetVal += printResult(MetaData::has(FieldName), false,
                         "'" + FieldName + "' is destroyed correctly",
                         "'" + FieldName + "' is not destroyed");

   // test MetaDim::destroy()
   MetaDim::destroy(DimName);

   RetVal += printResult(MetaDim::has(DimName), false,
                         "'" + DimName + "' is destroyed correctly",
                         "'" + DimName + "' is not destroyed");

   return RetVal;
}

int main(int argc, char **argv) {

   int RetVal = 0;

   try {

      RetVal += testMetaDim();

      RetVal += testMetaData();

      RetVal += testMetaGroup();

      RetVal += testMetaInit();

   } catch (const std::exception &Ex) {
      std::cout << Ex.what() << ": FAIL" << std::endl;
      RetVal += 1;

   } catch (...) {
      std::cout << "Unknown: FAIL" << std::endl;
      RetVal += 1;
   }

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;
}
