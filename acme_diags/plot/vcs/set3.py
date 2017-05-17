import vcs
import cdms2

def plot(ref, test, diff):
    vcs_canvas = vcs.init()
    clean_template = vcs.createtemplate('clean_template')
    clean_template.blank(["mean", "max", "min", "zvalue", "dataname", "crtime", "ytic2", "xtic2", "xname", "yname", "legend"])
    
    graph1 = vcs_canvas.createxvsy('ref_test_plot')
    graph1.datawc_y1 = 1
    graph1.datawc_y2 = 10
    graph1.datawc_x1 = 1
    vcs_canvas.plot(ref, graph1, clean_template)
    vcs_canvas.png('set3vcs.png')

if __name__ == '__main__':
    pth = '/Users/shaheen2/github/acme_diags/tests/'
    f = cdms2.open(pth + 'precc.nc')
    ref = f('PRECC')
    f.close()

    f = cdms2.open(pth + 'precl.nc')
    test = f('PRECL')
    f.close()

    diff = ref - test
    plot(ref, test, diff)