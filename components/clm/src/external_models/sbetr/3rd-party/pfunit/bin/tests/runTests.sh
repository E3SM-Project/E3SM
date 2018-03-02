python testParser.py

mkdir -p outputs
for file in simple beforeAfter TestA TestCaseA MpiTestCaseB ParameterizedTestCaseB MpiParameterizedTestCaseC
do
   ../pFUnitParser.py inputs/${file}.pf outputs/${file}.F90
   diff outputs/${file}.F90 expectedOutputs/${file}.F90
done
