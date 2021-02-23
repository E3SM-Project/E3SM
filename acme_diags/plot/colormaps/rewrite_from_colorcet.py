try:
    import colorcet
except BaseException:
    print("Cannot convert from colorcet w/o colorcet")
    import sys

    sys.exit()

all_cms = colorcet.cm


def dump_cmap(name, mpl_cmap):
    nm = "cet_%s" % name
    with open("%s.rgb" % nm, "w") as f:
        f.write("# Converted from colorcet\n")
        f.write("#\n")
        f.write("# number of colors in table\n")
        f.write("#ncolors = %i\n" % mpl_cmap.N)
        f.write("#\n")
        f.write("#  r   g   b\n")
        for i in range(mpl_cmap.N):
            a = float(i) / float(mpl_cmap.N - 1)
            r, g, b, a = [int(x * 255) for x in mpl_cmap(a)]
            f.write(" %3s %3s %3s\n" % (r, g, b))
    print("Wrote %s" % nm)


for cmap in list(all_cms.keys()):
    dump_cmap(cmap, all_cms[cmap])
