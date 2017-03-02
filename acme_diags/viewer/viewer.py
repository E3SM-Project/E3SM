import os
import stat
from output_viewer.build import build_viewer
from output_viewer.utils import rechmod
from output_viewer.index import OutputIndex, OutputPage, OutputFile, OutputRow, OutputGroup
from acme_diags.acme_parameter import ACMEParameter

parameter = ACMEParameter()
parameter.case_id = 'set5_PRECT_GPCP'
parameter.reference_name = 'GPCP (yrs1979-2009)'
parameter.var = 'PRECT'

# Use the commented code when this is called by the driver
# path = os.path.abspath(parameter.case_id)
path = os.path.join('/Users/shaheen2/github/acme_diags/acme_diags/plotset5', parameter.case_id)

index = OutputIndex(parameter.reference_name)

page1 = OutputPage("Set 5", ['Description', 'ANN', 'DJF', 'JJA', 'MAM', 'SON'])

cols = []
cols.append('Some description for %s' % parameter.case_id)

files_in_dir = [fnm for fnm in os.listdir(path) if '.png' in fnm]
print files_in_dir

for f in files_in_dir:
    cols.append(OutputFile(f))

r1 = OutputRow(parameter.var, cols)

g1 = OutputGroup('Variables for this obs')
#g2 = OutputGroup('Some Other Group')

page1.addGroup(g1)
#page1.addGroup(g2)
page1.addRow(r1, 0)


index.addPage(page1)
#index.addPage(page2)
index.toJSON(os.path.join(path, "index.json"))






default_mask = stat.S_IMODE(os.stat(path).st_mode)
rechmod(path, default_mask)

if os.access(path, os.W_OK):
    default_mask = stat.S_IMODE(os.stat(path).st_mode)  # mode of files to be included
    build_viewer(os.path.join(path, "index.json"), diag_name="ACME Diagnostics", default_mask=default_mask)


if os.path.exists(os.path.join(path, "index.html")):
    should_open = raw_input("Viewer HTML generated at %s/index.html. Would you like to open in a browser? y/[n]: " % path)
    if should_open and should_open.lower()[0] == "y":
        import webbrowser
        webbrowser.open("file://" + os.path.join(path, "index.html"))
else:
    print "Failed to generate the viewer."
    sys.exit(1)
