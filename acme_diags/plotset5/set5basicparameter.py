#reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'

test_name = '20160520.A_WCYCL1850.ne30'
results_dir = 'testresults'
regrid_tool = 'esmf'
regrid_method = 'linear'

#backend = 'cartopy'
backend = 'vcs'

results_dir='./numbers/'
diff_title = 'test - reference'
diff_colormap = 'bl_to_darkred'

# output_format = ['png', 'pdf']

## To add custmer derived variable:
## 1. create a function for calculating dervied variable
#
#def albedo(solin, fsntoa):
#    """TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension"""
#    var = (solin - fsntoa) / solin
#    var.units = "dimensionless"
#    var.long_name = "TOA albedo"
#    return var
#
## 2. generate a dictionary which can be appended to the default derived_variables dictionary.
# 
##derived_variables = {
##    'derived_var_key': [
##        (['output_var1','output_var2'], lambda var1,var2: function(var1,var2))
##    ]
##}
#
#derived_variables = {
#    'ALBEDO': [
#        (['SOLIN', 'FSNTOA'], lambda solin, fsntoa: albedo(solin, fsntoa))
#    ]
#}
#
## 3. In mydiags.json: generate a new entry for 'derived_var_key'
##For example:
##{
##    "set5": [
##          {
##      "case_id": "set5_CERES-EBAF",
##      "variables": "ALBEDO",
##      "ref_name": "CERES-EBAF",
##      "reference_name": "CERES-EBAF March 2000-Feb 2013",
##      "season": ["ANN"],
##                }
##              ]
##} 


