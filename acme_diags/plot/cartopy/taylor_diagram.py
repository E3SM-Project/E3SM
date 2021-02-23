# Copyright: Modified from below document which has  been placed in the public domain.  #Jill Zhang, Dec 2017

"""
Taylor diagram (Taylor, 2001) test implementation.
http://www-pcmdi.llnl.gov/about/staff/Taylor/CV/Taylor_diagram_primer.htm
"""

__version__ = "Time-stamp: <2012-02-17 20:59:35 ycopin>"
__author__ = "Yannick Copin <yannick.copin@laposte.net>"

import matplotlib
import numpy as np

matplotlib.use("Agg")


class TaylorDiagram(object):
    """Taylor diagram: plot model standard deviation and correlation
    to reference (data) sample in a single-quadrant polar plot, with
    r=stddev and theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=111, label="_"):
        """Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using mpl_toolkits.axisartist.floating_axes. refstd is
        the reference standard deviation to be compared to.
        """

        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF
        from matplotlib.projections import PolarAxes

        self.refstd = refstd  # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.concatenate((np.arange(10) / 10.0, [0.95, 0.99]))
        tlocs = np.arccos(rlocs)  # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)  # Positions
        gl2_num = np.linspace(0, 1.5, 7)
        gl2 = GF.FixedLocator(gl2_num)
        gl2_ticks = [(x, str(x)) for x in gl2_num]
        gl2_ticks[-1] = [gl2_num[-1], ""]  # type: ignore
        gl2_ticks[0] = [gl2_num[0], "0"]  # type: ignore
        tf1 = GF.DictFormatter(dict(list(zip(tlocs, list(map(str, rlocs))))))
        tf2 = GF.DictFormatter(dict(gl2_ticks))

        # Standard deviation axis extent
        self.smin = min(gl2_num)
        # self.smax = 1.5*self.refstd
        self.smax = max(gl2_num)

        ghelper = FA.GridHelperCurveLinear(
            tr,
            extremes=(0, np.pi / 2, self.smin, self.smax),  # 1st quadrant
            grid_locator1=gl1,
            grid_locator2=gl2,
            tick_formatter1=tf1,
            tick_formatter2=tf2,
        )

        #        if fig is None:
        #            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Standard deviation (Normalized)")

        ax.axis["right"].set_axis_direction("top")  # "Y axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left")

        ax.axis["bottom"].set_visible(False)  # Useless

        # Contours along standard deviations
        ax.grid(False)

        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height * 0.8])

        # Put a legend to the right of the current axis
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

        self._ax = ax  # Graphical axes
        self.ax = ax.get_aux_axes(tr)  # Polar coordinates

        # Add reference point and stddev contour
        # print "Reference std:", self.refstd
        (l,) = self.ax.plot([0], self.refstd, "k*", ls="", ms=10, label=label)
        t = np.linspace(0, np.pi / 2)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, "k--", label="_")

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """Add sample (stddev,corrcoeff) to the Taylor diagram. args
        and kwargs are directly propagated to the Figure.plot
        command."""

        (l,) = self.ax.plot(
            np.arccos(corrcoef), stddev, *args, **kwargs
        )  # (theta,radius)
        self.samplePoints.append(l)

        return l

    def add_contours(self, levels=5, **kwargs):
        """Add constant centered RMS difference contours."""

        rs, ts = np.meshgrid(
            np.linspace(self.smin, self.smax), np.linspace(0, np.pi / 2)
        )
        # Compute centered RMS difference
        rms = np.sqrt(self.refstd ** 2 + rs ** 2 - 2 * self.refstd * rs * np.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours
