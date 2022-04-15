#!/usr/bin/env python

import sys, string, copy, glob, os, json, argparse, textwrap, functools
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as pl

def is_none(a): return a == None
def is_empty(a): return len(a) == 0
def is_zero(a): return a == 0
def is_pod_number(n): return type(n) in (int, float)
def first(coll): return coll[0]
def second(coll): return coll[1]
def last(coll): return coll[-1]
def mapl(fun, *args): return list(map(fun, *args))

def sort(coll):
    c = copy.deepcopy(coll)
    c.sort()
    return c

def get(d, keys):
    try:
        for k in keys: d = d[k]
    except: d = d[keys]
    return d

def geton(d, keys):
    try: return get(d, keys)
    except: return None

def assoc_nested(d, keys, val):
    # Associate in a nested dict, creating new sub-dicts as needed.
    dn = d
    for key in keys[0:-1:]:
        if not key in dn: dn[key] = {}
        dn = dn[key]
    dn[keys[-1]] = val

def assoc_nested_append(d, keys, val):
    # Associate in a nested dict, creating new sub-dicts as needed. The value
    # is intended to be an item to go into a list. If the list exists, append
    # to it; if not, create it with the item as the only element.'
    try:
        v = get(d, keys)
    except:
        v = []
    assoc_nested(d, keys, v + [val])

def unit_test():
    d = {}
    keys = ('hi', 1, 'bye', 3.5)
    assoc_nested_append(d, keys, 4)
    assoc_nested_append(d, keys, 5)
    assert(geton(d, keys) == [4,5])

def dispfig(fn_prefix=None, format='pdf', tight=True):
    if tight: pl.tight_layout()
    if not fn_prefix or is_empty(fn_prefix):
        return pl.show()
    else:
        pl.savefig(fn_prefix + '.' + format, format=format, bbox_inches='tight')

def my_grid():
    pl.grid(True, lw=0.5, ls='-', color=(0.8, 0.8, 0.8), zorder=-1, which='both')
    return pl.gca().set_axisbelow(True)

def pad_lim(lim, pad=0.05, mult=False):
    if mult:
        v = first(lim) * (1 - pad), second(lim) * (1 + pad)
    else:
        d = second(lim) - first(lim)
        delta = pad * d
        v = first(lim) - delta, second(lim) + delta
    return v

def axis_tight_pad(pad=0.05, mult=False):
    pl.axis('tight')
    xl = pl.xlim()
    yl = pl.ylim()
    pl.xlim(pad_lim(xl, pad, mult))
    return pl.ylim(pad_lim(yl, pad, mult))

class pl_plot:
    def __init__(me, figsize, filename, format=None, tight=True):
        me.filename = filename
        me.format = 'pdf' if is_none(format) else format
        me.tight = tight
        pl.close()
        pl.figure(num=1, figsize=figsize)
    def cleanup(me):
        dispfig(me.filename, format=me.format, tight=me.tight)
    def __enter__(me): return me
    def __exit__(me, *args): pass
    def __del__(me): return me.cleanup()

class ClimoPlotter:
    def __init__(me):
        me.d = {}

    def parse(me, dirname):
        json_fns = glob.glob(dirname + '/*.json')
        for jfn in json_fns:
            _, name = os.path.split(jfn)
            stuff = name.split('-')
            region = first(last(stuff).split('.'))
            timespan = stuff[-2]
            name = functools.reduce(lambda a, e: a + '-' + e, stuff[1:-2:])
            try:
                with open(jfn, 'r') as fp:
                    d = json.load(fp)
                    assoc_nested(me.d, (name, region, timespan), d)
            except:
                print('Failed to parse ' + jfn)

    def plot(me, filename, exclude_names=None, save_format=None, ann_only=False,
             ylim_corr=None, ylim_rmse=None, ylim_mean=None, ylim_extr=None):
        if is_none(save_format): save_format = 'pdf'
        if is_none(exclude_names): exclude_names = []
        dx = 0.1
        nplot = 4
        names = sort(list(filter(lambda e: not e in exclude_names, me.d.keys())))
        with pl_plot((11, 9), filename, format=save_format):
            axs = []
            for ploti in range(nplot):
                axs.append(pl.subplot(nplot, 1, ploti+1))
            for ni, name in enumerate(names):
                try:
                    d1 = me.d[name]['ocean' if name in ('TAUXY',) else
                        'global']
                except:
                    d1 = None
                if is_none(d1): continue
                tss = d1.keys()
                nts = len(tss)
                x = ni + 0.5 + (dx * -2 if nts > 1 else 0)
                it = 'ANN', 'DJF', 'MAM', 'JJA', 'SON'
                it = it[0:1] if nts == 1 or ann_only else it
                tss_avail = it
                for ti, ts in enumerate(tss_avail):
                    for ploti in range(nplot):
                        pl.sca(axs[ploti])
                        d2 = geton(d1, ts)
                        if is_none(d2): continue
                        if ploti == 0: y = d2['misc']['corr']
                        elif ploti == 1: y = d2['misc']['rmse']
                        elif ploti == 2: y = d2['diff']['mean']
                        elif ploti == 3: y = (d2['test']['min'] - d2['ref']['min'],
                                              d2['test']['max'] - d2['ref']['max'])
                        else: y = None
                        if is_zero(ploti): den = 1
                        elif ploti in (1, 2): den = d2['ref_regrid']['std']
                        else: den = d2['ref']['max'] - d2['ref']['min']
                        color = 'kbgrm'[ti]
                        if is_zero(den): continue
                        if is_pod_number(y):
                            pl.plot(x, 100 * (y / den), color + ('o' if is_zero(ni % 2) else 's'))
                        else:
                            pl.plot(x, 100 * (second(y) / den), color + '^')
                            pl.plot(x, 100 * (first(y) / den), color + 'v', fillstyle='none')
                    x += dx
            for ploti in range(nplot):
                pl.sca(axs[ploti])
                axis_tight_pad()
                my_grid()
                xticks = mapl(lambda e: 0.5 + e, range(len(names)))
                if ploti == nplot-1: pl.xticks(xticks, names, rotation=90)
                else: pl.xticks(xticks, [])
                if   not is_none(ylim_corr) and ploti == 0: pl.ylim(ylim_corr, 101)
                elif not is_none(ylim_rmse) and ploti == 1: pl.ylim(-1, ylim_rmse)
                elif not is_none(ylim_mean) and ploti == 2: pl.ylim(-ylim_mean, ylim_mean)
                elif not is_none(ylim_extr) and ploti == 3: pl.ylim(-ylim_extr, ylim_extr)
                pl.xlim(0, len(names))
                if ploti == 0: t = 'Correlation (%)'
                elif ploti == 1: t = 'RMSE normalized w.r.t. std dev (%)'
                elif ploti == 2: t = 'Diff mean normalized w.r.t. std dev (%)'
                elif ploti == 3: t = 'Min, Max Diff normalized w.r.t. range (%)'
                else: t = None
                pl.title(t)

def main():
    help = textwrap.dedent("""
    Make a figure that summarizes e3sm_diags lat_lon json output.

    Read the lat_lon *.json files output by e3sm_diags and make a figure that
    summarizes these statistical data.

    Colors: black: ANN, blue: DJF, green: MAM, red: JJA, majenta: SON

    Squares and circles are to help with visual separation of points but
    otherwise have no meaning.

    Filled up-pointing triangle is max; open down-pointing triangle is min.

    RMSE and mean-difference values are normalized by ref_regrid std.

    Extrema values are normalized by ref (max - min).
    """)
    p = argparse.ArgumentParser(description=help, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('jsondir', type=str,
                   help='Path to e3sm_diags lat_lon *.json files')
    p.add_argument('figname', type=str,
                   help='Name of output image, excluding format suffix')
    p.add_argument('--format', type=str, default='png',
                   help='Image format (default png)')
    p.add_argument('--ann-only', action='store_true',
                   help='Plot ANN data only')
    for plot in ('corr', 'rmse', 'mean', 'extrema'):
        p.add_argument('--ylim-{}'.format(plot), type=float, default=0,
                       help='Y limit for {} plot; <=0 for none (default)'.format(plot))
    o = p.parse_args()
    c = ClimoPlotter()
    c.parse(o.jsondir)
    c.plot(o.figname,
           ann_only=o.ann_only,
           ylim_corr=o.ylim_corr if o.ylim_corr > 0 else None,
           ylim_rmse=o.ylim_rmse if o.ylim_rmse > 0 else None,
           ylim_mean=o.ylim_mean if o.ylim_mean > 0 else None,
           ylim_extr=o.ylim_extrema if o.ylim_extrema > 0 else None,
           save_format=o.format)

if __name__ == '__main__':
    unit_test()
    main()
