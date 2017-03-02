import os
import stat
from output_viewer.build import build_viewer
from output_viewer.utils import rechmod
from output_viewer.index import OutputIndex, OutputPage, OutputFile, OutputRow, OutputGroup

path = '/Users/shaheen2/github/acme_diags/acme_diags/plotset5'

ind = OutputIndex("Name of the model")

page1 = OutputPage("Set 5", ['Description', 'ANN', 'DJF', 'JJA', 'MAM', 'SON'])

cols = []
cols.append('some description')
f1 = OutputFile(os.path.join(path, 'set5_PRECT_GPCP/GPCP_ANN.png'))
cols.append(f1)
f1 = OutputFile(os.path.join(path, 'set5_PRECT_GPCP/GPCP_DJF.png'))
cols.append(f1)
f1 = OutputFile(os.path.join(path, 'set5_PRECT_GPCP/GPCP_JJA.png'))
cols.append(f1)
f1 = OutputFile(os.path.join(path, 'set5_PRECT_GPCP/GPCP_MAM.png'))
cols.append(f1)
f1 = OutputFile(os.path.join(path, 'set5_PRECT_GPCP/GPCP_SON.png'))
cols.append(f1)

r1 = OutputRow('PRECT', cols)

g1 = OutputGroup('Variables for this obs')
g2 = OutputGroup('Some Other Group')

page1.addGroup(g1)
page1.addGroup(g2)
page1.addRow(r1, 0)


ind.addPage(page1)
#ind.addPage(page2)
ind.toJSON(os.path.join(path, "index.json"))






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
