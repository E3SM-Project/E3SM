import os
import stat
from output_viewer.build import build_viewer
from output_viewer.utils import rechmod
from output_viewer.index import OutputIndex, OutputPage, OutputFile, OutputRow, OutputGroup

path = '/Users/shaheen2/github/acme_diags/acme_diags/plotset5/'

ind = OutputIndex("Test Package")

page1 = OutputPage("Page 1")
page2 = OutputPage("Page 2")

'''
f1 = OutputFile(os.path.join(path, 'set5_ANN_PRECT_TRMM/test.png'))

r1 = OutputRow('Row', f1)

g1 = OutputGroup('Group')

page1.addGroup(g1)
page1.addRow(r1, 0)
'''

ind.addPage(page1)
ind.addPage(page2)
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
