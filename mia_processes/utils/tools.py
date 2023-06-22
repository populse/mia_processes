# -*- coding: utf-8 -*-
"""
Module that contains multiple functions used across mia_processes

:Contains:
    :Class:
        - PageNumCanvas
        - ReportLine
    :Functions:
        - checkFileExt
        - dict4runtime_update
        - get_dbFieldValue
        - mriqc_get_all_run
        - mriqc_group_iqms_tsv
        - plot_boxplot_points
        - plot_qi2
        - set_dbFieldValue
        - slice_planes_plot
"""

##########################################################################
# Populse_mia - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import datetime
import glob
import json
import os

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.cm import get_cmap
from mpl_toolkits.axes_grid1 import ImageGrid
from nilearn.plotting import plot_anat
from nitransforms.io.afni import _dicom_real_to_card
from populse_mia.data_manager.data_history_inspect import data_history_pipeline
from populse_mia.data_manager.project import (
    COLLECTION_CURRENT,
    COLLECTION_INITIAL,
)
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas
from reportlab.platypus import Flowable
from traits.api import Undefined


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
        self.drawRightString(195 * mm, 10 * mm, page)


class ReportLine(Flowable):
    """Line flowable --- draws a line in a flowable"""

    def __init__(self, width, height=0):
        Flowable.__init__(self)
        self.width = width
        self.height = height

    def __repr__(self):
        """Customise the string representation of the ReportLine object

        :returns: the string representation of the object
        """
        return "Line(w=%s)" % self.width

    def draw(self):
        "Draw the line"
        self.canv.line(0, self.height, self.width, self.height)


def checkFileExt(in_file, ext_dic):
    """Check file extension

    :param in_file: file name (a string)
    :param ext_dic: dictionary of the valid extensions for the file
                    (dictionary, ex:
                    EXT = {'NIFTI_GZ': 'nii.gz', 'NIFTI': 'nii'})
    :returns:
        - valid_bool: True if extension is valid (a boolean)
        - in_ext: file extension (a string)
        - file_name: file name without extension (a string)
    """

    # Get file extension
    valid_bool = False
    ifile = os.path.split(in_file)[-1]
    file_name, in_ext = ifile.rsplit(".", 1)
    if in_ext == "gz":
        (file_name_2, in_ext_2) = file_name.rsplit(".", 1)
        in_ext = in_ext_2 + "." + in_ext
        file_name = file_name_2

    valid_ext = list(ext_dic.values())

    if in_ext in valid_ext:
        valid_bool = True

    return valid_bool, in_ext, file_name


def dict4runtime_update(dict4runtime, database, db_filename, *args):
    """blabla"""

    for tag in args:
        if tag in database.get_fields_names(COLLECTION_CURRENT):
            dict4runtime[tag] = database.get_value(
                COLLECTION_CURRENT, db_filename, tag
            )

        else:
            dict4runtime[tag] = "Undefined"

        if isinstance(dict4runtime[tag], datetime.date):
            dict4runtime[tag] = str(dict4runtime[tag])

        if dict4runtime[tag] is None:
            dict4runtime[tag] = "Undefined"


def get_dbFieldValue(project, document, field):
    """Return, for a document, the field value from the database.

    :param project: the project.
    :param document: the absolute path of the document.
    :param : the field name.
    :returns: the value of the field for the document in the database.
    """
    file_position = (
        document.find(project.getName()) + len(project.getName()) + 1
    )
    database_filename = document[file_position:]

    return project.session.get_value(
        COLLECTION_CURRENT, database_filename, field
    )


def mriqc_get_all_run(modality, project, output_directory):
    """
    Get all raw files name with a mriqc run for one project and
    for one modality.

    :param modality:  modality (a string, 'bold' or 'anat')
    :param project:  the project (Project Object)
    :param output_directory: output directory where all json files are stored
                             (a string that representing a path)
    :returns:
        - files_names (a list)
    """
    # Get all json
    jsonfiles = glob.glob(
        os.path.join(output_directory, "*" + modality + "_qc.json")
    )
    if not jsonfiles:
        return None
    files_name = []
    for jsonfile in jsonfiles:
        file_position = (
            jsonfile.find(project.getName()) + len(project.getName()) + 1
        )
        json_database_filename = jsonfile[file_position:]
        try:
            file_name = data_history_pipeline(
                json_database_filename, project
            ).in_file
        except Exception:
            # if we removed all mriqc data execpt json and pdf
            file_name = json_database_filename
        files_name.append(file_name)
    return files_name


def mriqc_group_iqms_tsv(modality, output_directory):
    """
    Get all IQMs from mriqc json report for one project
    and for one modality and put them together in one tsv file.

    :param modality:  modality (a string, 'bold' or 'anat')
    :param output_directory: output directory where all json files are stored
                             (a sting that representing a path)
    :returns:
        - dataframe: all IQMs in a panda dataframe
        - out_tsv: out tsv file (a string that representing a file)
    """
    date = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S_%f")[:22]
    # Get all json
    jsonfiles = glob.glob(
        os.path.join(output_directory, "*" + modality + "_qc.json")
    )
    if not jsonfiles:
        return None, None

    tags_to_remove = [
        "histogram_qi2_x_grid",
        "histogram_qi2_ref_pdf",
        "histogram_qi2_fit_pdf",
        "histogram_qi2_ref_data",
        "histogram_qi2_cutoff_idx",
        "vec_outliers",
        "vec_fd",
        "vec_dvars",
    ]
    datalist = []
    for jsonfile in jsonfiles:
        with open(jsonfile) as f:
            data = json.load(f)
            for tag in tags_to_remove:
                if tag in list(data.keys()):
                    data.pop(tag)
            for tag in list(data.keys()):
                if "vec_spikes" in tag:
                    data.pop(tag)
            data["file_name"] = (
                jsonfile.split("/")[-1]
                .replace("anat_qc.json", "")
                .replace("bold_qc.json", "")
                .replace("mean_reg_", "")
            )
            datalist.append(data)

    dataframe = pd.DataFrame(datalist)
    cols = dataframe.columns.tolist()
    dataframe = dataframe.sort_values(by=["file_name"])

    out_tsv = os.path.join(
        output_directory, "mriqc_group_report_" + modality + date + ".tsv"
    )

    # Set filename at front
    cols.insert(0, cols.pop(cols.index("file_name")))
    dataframe[cols].to_csv(str(out_tsv), index=False, sep="\t")

    return dataframe, out_tsv


def plot_boxplot_points(dataframe, title, ylabel, out_file=None):
    """
    Plot boxplot with points data and save the figure in a png image

    :param dataframe:  tabular data. (a pandas dataframe)
    :param title: figure title (a string)
    :param ylabel: y axis label (a string)
    :param out_file: out figure path (a string)

    :returns:
        - out_file: out figure path (a string)
    """

    vals, names, xs = [], [], []
    for i, col in enumerate(dataframe.columns):
        vals.append(dataframe[col].values)
        names.append(col)
        xs.append(
            np.random.normal(i + 1.3, 0.04, dataframe[col].values.shape[0])
        )

    plt.figure()
    plt.boxplot(vals, sym="", labels=names)

    palette = ["r", "g", "b", "y", "c", "m"]
    while len(palette) < len(vals):
        palette += palette

    for x, val, c in zip(xs, vals, palette):
        plt.scatter(x, val, alpha=0.4, color=c)

    ax = plt.gca()
    ax.set(axisbelow=True, title=title, ylabel=ylabel)
    ax.tick_params(axis="x", labelrotation=45)
    if out_file is None:
        out_file = os.path.abspath("box_plot.png")
    plt.savefig(out_file, format="png", bbox_inches="tight")

    return out_file


def plot_qi2(x_grid, ref_pdf, fit_pdf, ref_data, cutoff_idx, out_file=None):
    """bla bla"""
    fig, ax = plt.subplots()

    ax.plot(
        x_grid,
        ref_pdf,
        linewidth=2,
        alpha=0.5,
        label="background",
        color="dodgerblue",
    )

    refmax = np.percentile(ref_data, 99.95)
    x_max = x_grid[-1]

    ax.hist(
        ref_data,
        40 * max(int(refmax / x_max), 1),
        fc="dodgerblue",
        histtype="stepfilled",
        alpha=0.2,
        density=True,
    )

    fit_pdf[fit_pdf > 1.0] = np.nan

    ax.plot(
        x_grid,
        fit_pdf,
        linewidth=2,
        alpha=0.5,
        label="chi2",
        color="darkorange",
    )

    ylims = ax.get_ylim()

    ax.axvline(
        x_grid[-cutoff_idx],
        ymax=ref_pdf[-cutoff_idx] / ylims[1],
        color="dodgerblue",
    )

    plt.xlabel('Intensity within "hat" mask')
    plt.ylabel("Frequency")
    ax.set_xlim([0, x_max])
    plt.legend()

    if out_file is None:
        out_file = os.path.abspath("qi2_plot.png")

    fig.savefig(out_file, bbox_inches="tight", pad_inches=0, dpi=300)
    return out_file


def plot_segmentation(anat_file, segmentation, name, out_dir=None, **kwargs):
    """
    Adapted from
    <https://github.com/nipreps/mriqc/blob/5a0f0408bd0c176dbc46088c6ffe279269180f3f/mriqc/viz/utils.py#L550>
    """

    vmax = kwargs.get("vmax")
    vmin = kwargs.get("vmin")

    anat_ras = nib.as_closest_canonical(nib.load(anat_file))
    anat_ras_plumb = anat_ras.__class__(
        anat_ras.dataobj, _dicom_real_to_card(anat_ras.affine), anat_ras.header
    )

    seg_ras = nib.as_closest_canonical(nib.load(segmentation))
    seg_ras_plumb = seg_ras.__class__(
        seg_ras.dataobj, _dicom_real_to_card(seg_ras.affine), seg_ras.header
    )

    if kwargs.get("saturate", False):
        vmax = np.percentile(anat_ras.get_fdata().reshape(-1), 70)

    if vmax is None and vmin is None:
        vmin = np.percentile(anat_ras.get_fdata().reshape(-1), 10)
        vmax = np.percentile(anat_ras.get_fdata().reshape(-1), 99)

    disp = plot_anat(
        anat_ras_plumb,
        display_mode=kwargs.get("display_mode", "ortho"),
        cut_coords=kwargs.get("cut_coords", 8),
        title=kwargs.get("title"),
        vmax=vmax,
        vmin=vmin,
    )
    disp.add_contours(
        seg_ras_plumb,
        levels=kwargs.get("levels", [1]),
        colors=kwargs.get("colors", "r"),
    )

    fname, ext = os.path.splitext(os.path.basename(anat_file))

    if out_dir is None:
        out_file = os.path.abspath(fname + "_slice_planes_plot.png")

    else:
        out_file = os.path.join(out_dir, fname + name + ".png")

    disp.savefig(out_file)
    disp.close()
    disp = None

    return out_file


def set_dbFieldValue(project, document, tag_to_add):
    """Creates for a document, the field and the corresponding value in the db.

    :param project: the project.
    :param document: the absolute path of the document.
    :param tag_to_add: a dictionary with 'default_value', 'description',
                      'field_type', 'name', 'origin', 'unit', 'value',
                      'visibility' keys.
    """
    file_position = (
        document.find(project.getName()) + len(project.getName()) + 1
    )
    database_filename = document[file_position:]

    if tag_to_add["name"] not in project.session.get_fields_names(
        COLLECTION_CURRENT
    ):
        project.session.add_field(
            COLLECTION_CURRENT,
            tag_to_add["name"],
            tag_to_add["field_type"],
            tag_to_add["description"],
            tag_to_add["visibility"],
            tag_to_add["origin"],
            tag_to_add["unit"],
            tag_to_add["default_value"],
        )

    if tag_to_add["name"] not in project.session.get_fields_names(
        COLLECTION_INITIAL
    ):
        project.session.add_field(
            COLLECTION_INITIAL,
            tag_to_add["name"],
            tag_to_add["field_type"],
            tag_to_add["description"],
            tag_to_add["visibility"],
            tag_to_add["origin"],
            tag_to_add["unit"],
            tag_to_add["default_value"],
        )

    if project.session.get_document(COLLECTION_CURRENT, database_filename):
        print("Path {0} already in database.".format(database_filename))

    else:
        project.session.add_document(COLLECTION_CURRENT, database_filename)
        project.session.add_document(COLLECTION_INITIAL, database_filename)

    project.session.set_values(
        COLLECTION_CURRENT,
        database_filename,
        {tag_to_add["name"]: tag_to_add["value"]},
    )
    project.session.set_values(
        COLLECTION_INITIAL,
        database_filename,
        {tag_to_add["name"]: tag_to_add["value"]},
    )


def slice_planes_plot(
    data,
    fig_rows,
    fig_cols,
    inf_slice_start=None,
    slices_gap=None,
    dyn=1,
    cmap="Greys_r",
    out_dir=None,
    only_noise=False,
    out_name=None,
):
    "blablabla"

    brain_img = nib.as_closest_canonical(nib.load(data))
    brain_data = brain_img.get_fdata()
    brain_data = np.squeeze(brain_data)

    if len(brain_data.shape) == 4:
        brain_data = brain_data[:, :, :, dyn]

    elif len(brain_data.shape) > 4:
        # TODO: what we do in this case?
        pass

    min_thres = np.percentile(brain_data, 5)
    mask_data = np.ones_like(brain_data)
    mask_data[brain_data <= min_thres] = 0
    ind_non_zero = np.argwhere(mask_data)
    (ystart, xstart, zstart), (ystop, xstop, zstop) = (
        ind_non_zero.min(0),
        ind_non_zero.max(0) + 1,
    )
    brain_data = brain_data[ystart:ystop, xstart:xstop, zstart:zstop]
    disp_slices = fig_rows * fig_cols

    if inf_slice_start in (None, Undefined) and slices_gap in (
        None,
        Undefined,
    ):
        slices_gap = brain_data.shape[2] // disp_slices
        memory = set()
        ind_slices = None

        while (
            len(np.arange(start=0, stop=brain_data.shape[2], step=slices_gap))
            != disp_slices
        ):
            if (
                len(
                    np.arange(
                        start=0, stop=brain_data.shape[2], step=slices_gap
                    )
                )
                in memory
            ):
                ind_slices = np.arange(
                    start=0, stop=brain_data.shape[2], step=slices_gap
                )

                if len(ind_slices) > disp_slices:
                    while len(ind_slices) != disp_slices:
                        ind_slices = ind_slices[:-1]

                    ind_slices = ind_slices[:-1]
                    ind_slices = np.append(
                        ind_slices,
                        ind_slices[-1]
                        + (brain_data.shape[2] - ind_slices[-1]) // 2,
                    )
                    break

                if len(ind_slices) < disp_slices:
                    slices_gap -= 1

            elif (
                len(
                    np.arange(
                        start=0, stop=brain_data.shape[2], step=slices_gap
                    )
                )
                > disp_slices
            ):
                memory.add(
                    len(
                        np.arange(
                            start=0, stop=brain_data.shape[2], step=slices_gap
                        )
                    )
                )
                slices_gap += 1

            elif (
                len(
                    np.arange(
                        start=0, stop=brain_data.shape[2], step=slices_gap
                    )
                )
                < disp_slices
            ):
                memory.add(
                    len(
                        np.arange(
                            start=0, stop=brain_data.shape[2], step=slices_gap
                        )
                    )
                )
                slices_gap -= 1

        if ind_slices is None:
            inf_slice_start = (
                brain_data.shape[2]
                - np.arange(
                    start=0, stop=brain_data.shape[2], step=slices_gap
                )[
                    len(
                        np.arange(
                            start=0, stop=brain_data.shape[2], step=slices_gap
                        )
                    )
                    - 1
                ]
            ) // 2
            ind_slices = np.arange(
                start=inf_slice_start,
                stop=brain_data.shape[2],
                step=slices_gap,
            )

    elif len(np.array(list(range(0, brain_data.shape[2])))) > disp_slices:
        # fmt: off
        ind_slices = np.array(list(range(0, brain_data.shape[2])))[
            inf_slice_start:
            inf_slice_start + (slices_gap * disp_slices):
            slices_gap
        ]
        # fmt: on

    else:
        ind_slices = np.array(list(range(0, brain_data.shape[2])))

    # Reminder: 19cm == 7.4803inch; 23cm == 9.0551
    fig = plt.figure(figsize=(7.4803, 9.0551))  # Width, height in inches.

    mask_data = np.logical_not(np.isnan(brain_data))
    if only_noise:
        vmin = np.percentile(brain_data[mask_data], 0)
        vmax = np.percentile(brain_data[mask_data], 61)
    else:
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

        ax.imshow(
            np.swapaxes(brain_data[:, :, ind_slice], 0, 1),
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            interpolation="nearest",
            origin="lower",
            extent=[0, phys_sp[0], 0, phys_sp[1]],
        )
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(False)
        ax.axis("off")

        bgcolor = cmap(min(vmin, 0.0))
        fgcolor = cmap(vmax)

        ax.text(
            0.98,
            0.01,
            "%d" % ind_slice,
            color=fgcolor,
            transform=ax.transAxes,
            horizontalalignment="right",
            verticalalignment="bottom",
            fontsize=fontsize,
            bbox=dict(
                boxstyle="square,pad=0", edgecolor=bgcolor, facecolor=bgcolor
            ),
        )

    fname, ext = os.path.splitext(os.path.basename(data))
    if out_name:
        fname += "_" + out_name

    if out_dir is None:
        out_file = os.path.abspath(fname + "_slice_planes_plot.png")

    else:
        out_file = os.path.join(out_dir, fname + "_slice_planes_plot.png")

    fig.savefig(out_file, format="png", dpi=300, bbox_inches="tight")

    return out_file
