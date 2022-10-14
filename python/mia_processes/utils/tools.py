"""
Module that contains multiple functions used across mia_processes

:Contains:
    :Class:
        - PageNumCanvas
        - ReportLine
    :Functions:
        - dict4runtime_update
        - plot_qi2
        - slice_planes_plot
        -
"""

##########################################################################
# Populse_mia - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import datetime
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import os

from populse_mia.data_manager.project import COLLECTION_CURRENT

from matplotlib.cm import get_cmap
from mpl_toolkits.axes_grid1 import ImageGrid
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas
from reportlab.platypus import Flowable, Paragraph
from traits.api import Undefined

#from sys import exit
#from shutil import copyfile
# import readline as readlineComp
# import rlcompleter, getpass
# import logging # autocomp debug logfile

class PageNumCanvas(canvas.Canvas):
    """For add \"page number of total\" in each footer."""

    def __init__(self, *args, **kwargs):
        """Constructor."""
        canvas.Canvas.__init__(self, *args, **kwargs)
        self.pages = []

    def showPage(self):
        """On a page break, add information to the list."""
        self.pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        """Add the page number to each page (page x of y)."""
        page_count = len(self.pages)

        for page in self.pages:
            self.__dict__.update(page)
            self.draw_page_number(page_count)
            canvas.Canvas.showPage(self)

        canvas.Canvas.save(self)

    def draw_page_number(self, page_count):
        """Add the page number."""
        page = "Page %s of %s" % (self._pageNumber, page_count)
        self.setFont("Helvetica", 7)
        self.drawRightString(195*mm, 10*mm, page)


class ReportLine(Flowable):
    "Line flowable --- draws a line in a flowable"

    def __init__(self, width, height=0):
        Flowable.__init__(self)
        self.width = width
        self.height = height

    def __repr__(self):
        return "Line(w=%s)" % self.width

    def draw(self):
        "draw the line"
        self.canv.line(0, self.height, self.width, self.height)


def dict4runtime_update(dict4runtime, database, db_filename, *args):
    """blabla"""

    for tag in args:

        if tag in database.get_fields_names(COLLECTION_CURRENT):
            dict4runtime[tag] = database.get_value(COLLECTION_CURRENT,
                                                   db_filename,
                                                   tag)

        else:
            dict4runtime[tag] = "Undefined"

        if isinstance(dict4runtime[tag], datetime.date):
            dict4runtime[tag] = str(dict4runtime[tag])

def plot_qi2(x_grid, ref_pdf, fit_pdf, ref_data, cutoff_idx, out_file=None):
    """bla bla"""
    fig, ax = plt.subplots()

    ax.plot(
        x_grid,
        ref_pdf,
        linewidth=2,
        alpha=0.5,
        label="background",
        color="dodgerblue")

    refmax = np.percentile(ref_data, 99.95)
    x_max = x_grid[-1]

    ax.hist(
        ref_data,
        40 * max(int(refmax / x_max), 1),
        fc="dodgerblue",
        histtype="stepfilled",
        alpha=0.2,
        density=True)

    fit_pdf[fit_pdf > 1.0] = np.nan

    ax.plot(
        x_grid,
        fit_pdf,
        linewidth=2,
        alpha=0.5,
        label="chi2",
        color="darkorange")

    ylims = ax.get_ylim()

    ax.axvline(
        x_grid[-cutoff_idx],
        ymax=ref_pdf[-cutoff_idx] / ylims[1],
        color="dodgerblue")

    plt.xlabel('Intensity within "hat" mask')
    plt.ylabel("Frequency")
    ax.set_xlim([0, x_max])
    plt.legend()

    if out_file is None:
        out_file = os.path.abspath("qi2_plot.svg")

    fig.savefig(out_file, bbox_inches="tight", pad_inches=0, dpi=300)
    return out_file

def slice_planes_plot(data, fig_rows, fig_cols, inf_slice_start=None,
                      slices_gap=None, cmap="Greys_r", out_dir=None):
    "blablabla"

    brain_img = nib.as_closest_canonical(nib.load(data))
    brain_data = brain_img.get_fdata()
    brain_data = np.squeeze(brain_data)
    min_thres = np.percentile(brain_data, 5)
    mask_data = np.ones_like(brain_data)
    mask_data[brain_data <= min_thres] = 0
    ind_non_zero = np.argwhere(mask_data)
    (ystart, xstart, zstart), (ystop, xstop, zstop) = (ind_non_zero.min(0),
                                                      ind_non_zero.max(0) + 1)
    brain_data = brain_data[ystart:ystop, xstart:xstop, zstart:zstop]
    disp_slices = fig_rows * fig_cols

    if inf_slice_start in (None, Undefined) and slices_gap in (None, Undefined):
        slices_gap = brain_data.shape[2] // disp_slices
        memory = set()
        ind_slices = None

        while len(np.arange(start=0, stop=brain_data.shape[2],
                            step=slices_gap)) != disp_slices:

            if len(np.arange(start=0, stop=brain_data.shape[2],
                             step=slices_gap)) in memory:
                ind_slices = np.arange(start=0, stop=brain_data.shape[2],
                                       step=slices_gap)

                if len(ind_slices) > disp_slices:

                    while len(ind_slices) != disp_slices:
                        ind_slices = ind_slices[:-1]

                    ind_slices = ind_slices[:-1]
                    ind_slices = np.append(
                                    ind_slices,
                                    ind_slices[-1] + (
                                     brain_data.shape[2] - ind_slices[-1]) // 2)
                    break

                if len(ind_slices) < disp_slices:
                    slices_gap -= 1

            elif len(np.arange(start=0, stop=brain_data.shape[2],
                               step=slices_gap)) > disp_slices:
                memory.add(len(np.arange(start=0, stop=brain_data.shape[2],
                                         step=slices_gap)))
                slices_gap += 1

            elif len(np.arange(start=0, stop=brain_data.shape[2],
                               step=slices_gap)) < disp_slices:
                memory.add(len(np.arange(start=0, stop=brain_data.shape[2],
                                         step=slices_gap)))
                slices_gap -= 1

        if ind_slices is None:
            inf_slice_start = (brain_data.shape[2] -
                      np.arange(start=0,
                                stop=brain_data.shape[2],
                                step=slices_gap)[
                                      len(np.arange(start=0,
                                                    stop=brain_data.shape[2],
                                                    step=slices_gap)) - 1]) // 2
            ind_slices = np.arange(start=inf_slice_start,
                                   stop=brain_data.shape[2],
                                   step=slices_gap)

    elif len(np.array(list(range(0, brain_data.shape[2])))) > disp_slices:
        ind_slices = np.array(list(
                        range(0,
                        brain_data.shape[2])))[
                                   inf_slice_start:
                                   inf_slice_start + (slices_gap * disp_slices):
                                   slices_gap]

    else:
        ind_slices = np.array(list(range(0, brain_data.shape[2])))

    # Reminder: 19cm == 7.4803inch; 23cm == 9.0551
    fig = plt.figure(figsize=(7.4803, 9.0551))  # Width, height in inches.

    mask_data = np.logical_not(np.isnan(brain_data))
    vmin = np.percentile(brain_data[mask_data], 0.5)
    vmax = np.percentile(brain_data[mask_data], 99.5)

    zooms = brain_img.header.get_zooms()
    grid = ImageGrid(fig, 111, nrows_ncols=(fig_rows, fig_cols), axes_pad=0)
    cmap = get_cmap(cmap)

    if len(ind_slices) < 5:
        fontsize = 18

    elif len(ind_slices) < 10:
        fontsize = 16

    elif len(ind_slices) < 15:
        fontsize = 14

    elif len(ind_slices) < 20:
        fontsize = 12

    elif len(ind_slices) < 25:
        fontsize = 10

    elif len(ind_slices) < 40:
        fontsize = 8

    else:
        fontsize = 6

    for ax, ind_slice in zip(grid, ind_slices):
        phys_sp = np.array(zooms[:2]) * brain_data[:, :, ind_slice].shape

        ax.imshow(np.swapaxes(brain_data[:, :, ind_slice], 0, 1),
                  vmin=vmin,
                  vmax=vmax,
                  cmap=cmap,
                  interpolation="nearest",
                  origin="lower",
                  extent=[0, phys_sp[0], 0, phys_sp[1]])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(False)
        ax.axis("off")

        bgcolor = cmap(min(vmin, 0.0))
        fgcolor = cmap(vmax)

        ax.text(0.98,
                0.01,
                "%d" % ind_slice,
                color=fgcolor,
                transform=ax.transAxes,
                horizontalalignment="right",
                verticalalignment="bottom",
                fontsize=fontsize,
                bbox=dict(boxstyle="square,pad=0", edgecolor=bgcolor,
                          facecolor=bgcolor))

    fname, ext = os.path.splitext(os.path.basename(data))

    if out_dir is None:
        out_file = os.path.abspath(fname + "_slice_planes_plot.png")

    else:
        out_file = os.path.join(out_dir, fname + "_slice_planes_plot.png")

    fig.savefig(out_file, format="png", dpi=300, bbox_inches="tight")

    return out_file

