"""
Module that contains multiple functions used across mia_processes

:Contains:
    :Functions:
        - recupCover
"""

##########################################################################
# Populse_mia - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Paragraph
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
from reportlab.platypus import Flowable
import nibabel as nib
import numpy as np
from os.path import splitext, basename, join, abspath
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

#from os import listdir, system, makedirs, remove
#from os.path import isdir, isfile
#from datetime import datetime
#from sys import exit
#from shutil import copyfile
# import readline as readlineComp
# import rlcompleter, getpass
# import logging # autocomp debug logfile

class PageNumCanvas(canvas.Canvas):
    "For  add \"page number of total\" in each footer."

    def __init__(self, *args, **kwargs):
        "Constructor."
        canvas.Canvas.__init__(self, *args, **kwargs)
        self.pages = []

    def showPage(self):
        "On a page break, add information to the list."
        self.pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        "Add the page number to each page (page x of y)."
        page_count = len(self.pages)

        for page in self.pages:
            self.__dict__.update(page)
            self.draw_page_number(page_count)
            canvas.Canvas.showPage(self)

        canvas.Canvas.save(self)

    def draw_page_number(self, page_count):
        "Add the page number."
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


def recupCover(afile):
    """ Set up the coverage data for the report generated by MRIQC_report brick


    Read the \"afile\" file, split each line according to \":\"
    in a lx2 matrix.
    Then return the matrix ; column 0 = left from the colon ; column 1 = right
    from the colon.

    :param afile: A .txt file. Format == Parameter:Value.
    """

    styleSheet = getSampleStyleSheet()
    matrix = [[0, 0]]

    with open(afile, 'r') as txt:

        for line in txt:

            if bool(not line or line.isspace()):
                continue

            # Remove trailing characters
            line = line.rstrip()
            # Split the string according to the separator ":"
            header, doublePoint, info = line.partition(':')

            if matrix == [[0, 0]]:

                if header == " ":
                    matrix[0][0] = Paragraph("<para align=right>" + header +
                                             "</para>", styleSheet["Normal"])

                else:
                    # Class reportlab.platypus.paragraph.Paragraph
                    # with XML-like markup
                    matrix[0][0] = Paragraph("<para align=right><b>" + header +
                                             "&nbsp&nbsp&nbsp :</b></para>",
                                             styleSheet["BodyText"])
  
                matrix[0][1] = Paragraph("<para align=justify>" + info +
                                         "</para>", styleSheet["Normal"])

            else:

                if (header == " "):
                    temp = [Paragraph("<para align=right>" + header +
                                      "</para>", styleSheet["Normal"]),
                            Paragraph("<para align=justify>" + info +
                                      "</para>", styleSheet["Normal"])]

                else:
                    temp = [Paragraph("<para align=right><b>" + header +
                                      "&nbsp&nbsp&nbsp :</b></para>",
                                      styleSheet["BodyText"]),
                            Paragraph("<para align=justify>" + info +
                                      "</para>", styleSheet["Normal"])]

                matrix.append(temp)

    return matrix

def slice_planes_plot(data, fig_rows, fig_cols, inf_slice_start, slices_gap, cmap="Greys_r", out_dir=None):
    "blablabla"

    brain_img = nib.as_closest_canonical(nib.load(data))
    brain_data = brain_img.get_fdata()
    brain_data  = np.squeeze(brain_data)
    min_thres = np.percentile(brain_data, 5)
    mask_data = np.ones_like(brain_data)
    mask_data[brain_data <= min_thres] = 0
    ind_non_zero = np.argwhere(mask_data)
    (ystart, xstart, zstart), (ystop, xstop, zstop) = (ind_non_zero.min(0),
                                                      ind_non_zero.max(0) + 1)
    brain_data = brain_data[ystart:ystop, xstart:xstop, zstart:zstop]
    disp_slices = fig_rows * fig_cols

    if inf_slice_start == None and slices_gap == None:
        slices_gap = brain_data.shape[2] // disp_slices
        memory = set()
        ind_slices = None
        while len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)) != disp_slices:

            if len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)) in memory:
                ind_slices = np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)

                if len(ind_slices) > disp_slices:

                    while len(ind_slices) != disp_slices:
                        ind_slices = ind_slices[:-1]

                    ind_slices = ind_slices[:-1]
                    ind_slices = np.append(ind_slices, ind_slices[-1] + (brain_data.shape[2] - ind_slices[-1]) // 2)
                    break

                if len(ind_slices) < disp_slices:
                    slices_gap -= 1

            elif len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)) > disp_slices:
                memory.add(len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)))
                slices_gap += 1

            elif len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)) < disp_slices:
                memory.add(len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)))
                slices_gap -= 1

        if ind_slices is None:
            inf_slice_start = (brain_data.shape[2] - np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)[len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap)) - 1]) // 2
            ind_slices = np.arange(start=inf_slice_start, stop=brain_data.shape[2], step=slices_gap)

    elif len(np.array(list(range(0, brain_data.shape[2])))) > disp_slices:
        ind_slices = np.array(list(range(0, brain_data.shape[2])))[inf_slice_start:inf_slice_start + (slices_gap * disp_slices):slices_gap]

    else:
        ind_slices = np.array(list(range(0, brain_data.shape[2])))

    fig = plt.figure(figsize=(1.9 * fig_cols, 3 * fig_rows))

    mask_data = np.logical_not(np.isnan(brain_data))
    vmin = np.percentile(brain_data[mask_data], 0.5)
    vmax = np.percentile(brain_data[mask_data], 99.5)

    axis_numb = 1
    zooms = brain_img.header.get_zooms()

    for ind_slice in ind_slices:
        ax = fig.add_subplot(fig_rows, fig_cols, axis_numb)
        axis_numb += 1
        phys_sp = np.array(zooms[:2]) * brain_data[:, :, ind_slice].shape
        cmap = get_cmap(cmap)

        ax.imshow(np.swapaxes(brain_data[:, :, ind_slice], 0, 1),
                  vmin=vmin,
                  vmax=vmax,
                  cmap=cmap,
                  interpolation="nearest",
                  origin="lower",
                  extent=[0, phys_sp[0], 0, phys_sp[1]],)

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
                size=18,
                bbox=dict(boxstyle="square,pad=0", ec=bgcolor, fc=bgcolor))

    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
                        wspace=0.05, hspace=0.05)

    fname, ext = splitext(basename(data))
    if out_dir is None:
        out_file = abspath(fname + "_slice_planes_plot.png")
    else:
        out_file = join(out_dir, fname + "_slice_planes_plot.png")
    fig.savefig(out_file, format="png", dpi=300, bbox_inches="tight")

    return out_file

