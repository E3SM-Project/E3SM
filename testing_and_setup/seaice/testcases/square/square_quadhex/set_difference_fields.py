from netCDF4 import Dataset

#files = {"./output_hex_wachspress_0082x0094_120/output.2000.nc":
#         ["./output_hex_pwl_0082x0094_120/output.2000.nc",
#          "./output_hex_weak_0082x0094_120/output.2000.nc"],
#         "./output_quad_wachspress_0080x0080_120/output.2000.nc":
#         ["./output_quad_pwl_0080x0080_120/output.2000.nc",
#          "./output_quad_weak_0080x0080_120/output.2000.nc"]}

files = {"./output_hex_wachspress_0082x0094_120/output.2000.nc":
         ["./output/output.2000.nc"]}

#fieldnames = ["uVelocity","stressDivergenceU"]
#fieldnames = ["stressDivergenceU"]
fieldnames = ["uVelocity","vVelocity","stressDivergenceU","stressDivergenceV","strain11var","strain22var","strain12var"]


for filenameBase in files:
    for filenameDiff in files[filenameBase]:
        for fieldname in fieldnames:

            print(fieldname)

            filein = Dataset(filenameBase,"r")

            field1 = filein.variables[fieldname][:]

            dimensionsBase = filein.variables[fieldname].dimensions

            filein.close()


            fileDiff = Dataset(filenameDiff,"a")

            field2 = fileDiff.variables[fieldname][:]

            fieldDiff = fileDiff.createVariable(fieldname+"Diff", "d", dimensions=dimensionsBase)

            fieldDiff[:] = field2[:] - field1[:]

            fileDiff.close()
