# -*- coding: utf-8 -*-
"""
Module that contains multiple functions used across mia_processes

:Contains:
    :Class:
        - PageNumCanvas
        - ReportLine
    :Functions:
        - checkFileExt
        - del_dbFieldValue
        - dict4runtime_update
        - get_dbFieldValue
        - mriqc_get_all_run
        - mriqc_group_iqms_tsv
        - plot_boxplot_points
        - plot_qi2
        - plot_segmentation
        - plot_slice_planes
        - set_dbFieldValue
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
        """Draw the line"""
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


def del_dbFieldValue(project, document, tags2del):
    """Delete for a document, the value for a field in the db.

    :param project: the project.
    :param document: the absolute path of the document.
    :param tag_to_del: a list of fields where values will be deleted
    """

    file_position = (
        document.find(project.getName()) + len(project.getName()) + 1
    )
    database_filename = document[file_position:]

    for tag_to_del in tags2del:
        try:
            project.session.remove_value(
                COLLECTION_CURRENT, database_filename, tag_to_del
            )
        except ValueError:
            # The collection does not exist
            # or the field does not exist
            # or the document does not exist
            pass
        try:
            project.session.remove_value(
                COLLECTION_INITIAL, database_filename, tag_to_del
            )
        except ValueError:
            # The collection does not exist
            # or the field does not exist
            # or the document does not exist
            pass


def dict4runtime_update(dict4runtime, database, db_filename, *args):
    """Update the dict4runtime dictionary object with tags values

    :param dict4runtime: the dict used in mia_processes bricks to pass data
                         from the list_outputs method to the run_process_mia
                         method
    :param database: the database object (for example: self.project.session)
    :param db_filename: the name of the database file from which the tag value
                        will be retrieved
    :param args: the tags to be recovered

    """

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
    :param field: the field name.
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


def plot_slice_planes(
    data_1,
    data_2=None,
    fig_rows=4,
    fig_cols=4,
    slice_start=None,
    slice_step=None,
    dyn=1,
    cmap_1="Greys_r",
    vmin_1=None,
    vmax_1=None,
    cmap_2="rainbow",
    vmin_2=None,
    vmax_2=None,
    out_dir=None,
    only_noise=False,
    out_name=None,
):
    """Create a png file with a mosaic display of volume slices.

    fig_rows, fig_cols, slice_start and slice_step are optional. If they are
    not known, plot_slice_planes will try to adapt as best it can to obtain
    an acceptable display. The default values for fig_rows and
    fig_cols are `4` and for dyn it is '1'.

    :param data_1: an existing, image file (valid extensions: .nii).
    :param data_2: an existing, image file (valid extensions: .nii).
    :param fig_rows: the number of rows in the images panel (integer).
    :param fig_cols: the number of columns in the images panel (integer).
    :param slice_start: starting slice number to be displayed (integer).
    :param slice_step: the gap between each slice to be displayed (integer).
    :param dyn: the volume number in the space-time, for 4D (integer).
    :param cmap_1: the color map name (string) for data_1.
    :param vmin_1: the low value of the range used for data_1 display(float).
    :param vmax_1: the high value of the range used for data_1 display(float).
    :param cmap_2: the color map name (string) for data_2.
    :param vmin_2: the low value of the range used for data_2 display(float).
    :param vmax_2: the high value of the range used for data_2 display(float).
    :param out_dir: the output directory where the mosaic images will be saved.
    :param only_noise: if True, shows the noise (boolean).
    :param out_name: the suffix added to form the output name (string).

    Note: cmap for parametric display can be, gist_rainbow, RdYlBu, Spectral,
          rainbow_r, jet_r, seismic_r, bwr_r
    """

    brain_img_1 = nib.as_closest_canonical(nib.load(data_1))
    brain_data_1 = brain_img_1.get_fdata()
    brain_data_1 = np.squeeze(brain_data_1)

    if cmap_1 in (None, Undefined):
        cmap_1 = "Greys_r"

    cmap_1 = get_cmap(cmap_1)
    cmap_1.set_bad(color="black")  # Nan are black
    brain_data_2 = None

    if data_2 is not None:
        brain_img_2 = nib.as_closest_canonical(nib.load(data_2))
        brain_data_2 = brain_img_2.get_fdata()
        brain_data_2 = np.squeeze(brain_data_2)
        nan_indexes = np.isnan(brain_data_2)

        if cmap_2 in (None, Undefined):
            cmap_2 = "rainbow"

        cmap_2 = get_cmap(cmap_2)
        cmap_2.set_bad(alpha=0)
        # Remainder:
        # Amigo beta cmap vmin - vmax:0.01 - 0.25
        # Amigo SPMt cmap vmin - vmax:3 - 16

        if vmin_2 in (None, Undefined):
            vmin_2 = np.min(brain_data_2[~nan_indexes])

            if vmin_2 < 0:
                vmin_2 = 0.1

        if vmax_2 in (None, Undefined):
            vmax_2 = np.max(brain_data_2[~nan_indexes])

    if len(brain_data_1.shape) == 4:
        brain_data_1 = brain_data_1[:, :, :, dyn]

    elif len(brain_data_1.shape) > 4:
        # TODO: what we do in this case?
        pass

    min_thres = np.percentile(brain_data_1, 5)
    mask_data_1 = np.ones_like(brain_data_1)
    mask_data_1[brain_data_1 <= min_thres] = 0
    ind_non_zero = np.argwhere(mask_data_1)
    (ystart, xstart, zstart), (ystop, xstop, zstop) = (
        ind_non_zero.min(0),
        ind_non_zero.max(0) + 1,
    )
    brain_data_1 = brain_data_1[ystart:ystop, xstart:xstop, zstart:zstop]
    mask_data_1 = np.logical_not(np.isnan(brain_data_1))

    if only_noise:
        vmin_1 = np.percentile(brain_data_1[mask_data_1], 0)
        vmax_1 = np.percentile(brain_data_1[mask_data_1], 61)

    if vmin_1 in (None, Undefined):
        vmin_1 = np.percentile(brain_data_1[mask_data_1], 0.5)

    if vmax_1 in (None, Undefined):
        vmax_1 = np.percentile(brain_data_1[mask_data_1], 99.5)

    rowsxcols = disp_slices = fig_rows * fig_cols

    # The following part seems clearer and lighter than the part I'm just
    # commenting on just below. The commented part is too heavy and can
    # give unexpected results.
    # TODO: if this new part don't make issue we could delete the commented
    #      section completely.
    lst = np.arange(start=0, stop=brain_data_1.shape[2], step=1)

    if slice_start not in (None, Undefined):
        # fmt: off
        lst = lst[slice_start - 1 if slice_start != 0 else 0:]
        # fmt: on

    start = 0
    ind_slices = []

    if disp_slices > len(lst):
        disp_slices = len(lst)

    if slice_step in (None, Undefined):
        avg = len(lst) // disp_slices
        remainder = len(lst) % disp_slices

    else:
        avg = slice_step
        remainder = -1

    for i in range(disp_slices):
        end = start + avg + (1 if i < remainder else 0)

        if start >= len(lst):
            break

        ind_slices.append(lst[start:end])
        start = end

    if all(
        (isinstance(element, np.ndarray) and len(element) == 1)
        for element in ind_slices
    ):
        ind_slices = [element[0] for element in ind_slices]

    else:
        middle_idx = [(len(element) - 1) // 2 for element in ind_slices]
        ind_slices = [
            element[idx] for element, idx in zip(ind_slices, middle_idx)
        ]

    # # if slice_start and slice_step are unknown, we try to make a
    # # display of slices covering the entire volume, with equal spacing.
    # if slice_start in (None, Undefined) and slice_step in (
    #     None,
    #     Undefined,
    # ):
    #     slice_step = brain_data_1.shape[2] // disp_slices
    #
    #     if slice_step == 0:
    #         slice_step = 1
    #
    #     memory = set()
    #     ind_slices = None
    #
    #     while (
    #         len(np.arange(start=0,
    #                       stop=brain_data_1.shape[2], step=slice_step))
    #         != disp_slices
    #     ):
    #         if (
    #             len(
    #                 np.arange(
    #                     start=0, stop=brain_data_1.shape[2], step=slice_step
    #                 )
    #             )
    #             in memory
    #         ):
    #             ind_slices = np.arange(
    #                 start=0, stop=brain_data_1.shape[2], step=slice_step
    #             )
    #
    #             if len(ind_slices) > disp_slices:
    #                 while len(ind_slices) != disp_slices:
    #                     ind_slices = ind_slices[:-1]
    #
    #                 ind_slices = ind_slices[:-1]
    #                 ind_slices = np.append(
    #                     ind_slices,
    #                     ind_slices[-1]
    #                     + (brain_data_1.shape[2] - ind_slices[-1]) // 2,
    #                 )
    #                 break
    #
    #             if len(ind_slices) < disp_slices:
    #                 slice_step -= 1
    #
    #                 if len(np.arange(start=0,
    #                                  stop=brain_data_1.shape[2],
    #                                  step=1)) < disp_slices:
    #                     ind_slices = np.arange(start=0,
    #                                            stop=brain_data_1.shape[2],
    #                                            step=1)
    #                     break
    #
    #         elif (
    #             len(
    #                 np.arange(
    #                     start=0, stop=brain_data_1.shape[2], step=slice_step
    #                 )
    #             )
    #             > disp_slices
    #         ):
    #             memory.add(
    #                 len(
    #                     np.arange(
    #                         start=0,
    #                         stop=brain_data_1.shape[2], step=slice_step
    #                     )
    #                 )
    #             )
    #             slice_step += 1
    #
    #         elif (
    #             len(
    #                 np.arange(
    #                     start=0, stop=brain_data_1.shape[2], step=slice_step
    #                 )
    #             )
    #             < disp_slices
    #         ):
    #             memory.add(
    #                 len(
    #                     np.arange(
    #                         start=0,
    #                         stop=brain_data_1.shape[2], step=slice_step
    #                     )
    #                 )
    #             )
    #             slice_step -= 1
    #
    #             if slice_step == 0:
    #                 slice_step =1
    #
    #     if ind_slices is None:
    #         slice_start = (
    #             brain_data_1.shape[2]
    #             - np.arange(
    #                 start=0, stop=brain_data_1.shape[2], step=slice_step
    #             )[
    #                 len(
    #                     np.arange(
    #                         start=0,
    #                         stop=brain_data_1.shape[2], step=slice_step
    #                     )
    #                 )
    #                 - 1
    #             ]
    #         ) // 2
    #         ind_slices = np.arange(
    #             start=slice_start,
    #             stop=brain_data_1.shape[2],
    #             step=slice_step,
    #         )
    #
    # # if slice_start is known and slice_step is unknown, we try to make a
    # # display of slices covering the volume from slice_start to the end,
    # # with equal spacing.
    # elif slice_start not in (None, Undefined) and slice_step in (
    #     None,
    #     Undefined,
    # ):
    #     slice_step = (brain_data_1.shape[2] - slice_start) // disp_slices
    #
    #     if slice_step == 0:
    #         slice_step = 1
    #
    #     memory = set()
    #     ind_slices = None
    #
    #     while (
    #         len(
    #             np.arange(
    #                 start=slice_start,
    #                 stop=brain_data_1.shape[2],
    #                 step=slice_step,
    #             )
    #         )
    #         != disp_slices
    #     ):
    #         if (
    #             len(
    #                 np.arange(
    #                     start=slice_start,
    #                     stop=brain_data_1.shape[2],
    #                     step=slice_step,
    #                 )
    #             )
    #             in memory
    #         ):
    #             ind_slices = np.arange(
    #                 start=slice_start,
    #                 stop=brain_data_1.shape[2],
    #                 step=slice_step,
    #             )
    #
    #             if len(ind_slices) > disp_slices:
    #                 while len(ind_slices) != disp_slices:
    #                     ind_slices = ind_slices[:-1]
    #
    #                 ind_slices = ind_slices[:-1]
    #                 ind_slices = np.append(
    #                     ind_slices,
    #                     ind_slices[-1]
    #                     + (brain_data_1.shape[2] - ind_slices[-1]) // 2,
    #                 )
    #                 break
    #
    #             if len(ind_slices) < disp_slices:
    #                 slice_step -= 1
    #
    #         elif (
    #             len(
    #                 np.arange(
    #                     start=slice_start,
    #                     stop=brain_data_1.shape[2],
    #                     step=slice_step,
    #                 )
    #             )
    #             > disp_slices
    #         ):
    #             memory.add(
    #                 len(
    #                     np.arange(
    #                         start=slice_start,
    #                         stop=brain_data_1.shape[2],
    #                         step=slice_step,
    #                     )
    #                 )
    #             )
    #             slice_step += 1
    #
    #         elif (
    #             len(
    #                 np.arange(
    #                     start=slice_start,
    #                     stop=brain_data_1.shape[2],
    #                     step=slice_step,
    #                 )
    #             )
    #             < disp_slices
    #         ):
    #             memory.add(
    #                 len(
    #                     np.arange(
    #                         start=slice_start,
    #                         stop=brain_data_1.shape[2],
    #                         step=slice_step,
    #                     )
    #                 )
    #             )
    #             slice_step -= 1
    #
    #             if slice_step == 0:
    #                 slice_step = 1
    #                 break
    #
    #     if ind_slices is None:
    #         ind_slices = np.arange(
    #             start=slice_start,
    #             stop=brain_data_1.shape[2],
    #             step=slice_step,
    #         )
    #
    #     slice_step_bis = slice_step + 1
    #
    #     if (
    #         len(
    #             np.arange(
    #                 start=slice_start,
    #                 stop=brain_data_1.shape[2],
    #                 step=slice_step_bis,
    #             )
    #         )
    #         == disp_slices
    #     ):
    #         while (
    #             len(
    #                 np.arange(
    #                     start=slice_start,
    #                     stop=brain_data_1.shape[2],
    #                     step=slice_step_bis,
    #                 )
    #             )
    #             == disp_slices
    #         ):
    #             ind_slices_bis = np.arange(
    #                 start=slice_start,
    #                 stop=brain_data_1.shape[2],
    #                 step=slice_step_bis,
    #             )
    #
    #             if ind_slices_bis[-1] > ind_slices[-1]:
    #                 ind_slices = ind_slices_bis
    #
    #             slice_step_bis += 1
    #
    #     slice_step_bis = slice_step - 1
    #
    #     if slice_step_bis == 0:
    #         slice_step_bis = 1
    #
    #     if (
    #         len(
    #             np.arange(
    #                 start=slice_start,
    #                 stop=brain_data_1.shape[2],
    #                 step=slice_step_bis,
    #             )
    #         )
    #         == disp_slices
    #     ):
    #         while (
    #             len(
    #                 np.arange(
    #                     start=slice_start,
    #                     stop=brain_data_1.shape[2],
    #                     step=slice_step_bis,
    #                 )
    #             )
    #             == disp_slices
    #         ):
    #             ind_slices_bis = np.arange(
    #                 start=slice_start,
    #                 stop=brain_data_1.shape[2],
    #                 step=slice_step_bis,
    #             )
    #
    #             if ind_slices_bis[-1] > ind_slices[-1]:
    #                 ind_slices = ind_slices_bis
    #
    #             slice_step_bis -= 1

    # if slice_start and slice_step are known
    # elif (slice_start not in (None, Undefined) and
    #    slice_step not in (None, Undefined) and
    #    (len(np.array(list(range(0, brain_data_1.shape[2])))) > disp_slices)):
    #     # fmt: off
    #     ind_slices = np.array(list(range(0, brain_data_1.shape[2])))[
    #         slice_start:
    #         slice_start + (slice_step * disp_slices):
    #         slice_step
    #     ]
    #     # fmt: on
    #
    # else:
    #     ind_slices = np.array(list(range(0, brain_data_1.shape[2])))

    # Reminder: 19cm == 7.4803inch; 23cm == 9.0551
    fig = plt.figure(figsize=(7.4803, 9.0551))  # Width, height in inches.
    zooms = brain_img_1.header.get_zooms()

    if brain_data_2 is None:
        grid = ImageGrid(
            fig,
            111,
            nrows_ncols=(fig_rows, fig_cols),
            axes_pad=0,
        )
    else:
        grid = ImageGrid(
            fig,
            111,
            nrows_ncols=(fig_rows, fig_cols),
            axes_pad=0,
            cbar_mode="single",
            cbar_location="bottom",
            cbar_pad=0,
            cbar_size="2%",
        )

    if disp_slices < 5:
        fontsize = 18

    elif disp_slices < 10:
        fontsize = 16

    elif disp_slices < 15:
        fontsize = 14

    elif disp_slices < 20:
        fontsize = 12

    elif disp_slices < 25:
        fontsize = 10

    elif disp_slices < 40:
        fontsize = 8

    elif disp_slices < 100:
        fontsize = 6

    else:
        fontsize = 4

    rowsxcols_array = np.arange(0, rowsxcols, 1)

    for ax, slice_numb in zip(grid, rowsxcols_array):

        if slice_numb >= len(ind_slices):
            displ_1 = np.zeros((2, 2, 2))
            ind_slice = 1
            displ_2 = np.zeros((2, 2, 2))

        else:
            displ_1 = brain_data_1
            ind_slice = ind_slices[slice_numb]
            displ_2 = brain_data_2

        phys_sp = np.array(zooms[:2]) * brain_data_1[:, :, ind_slice].shape
        ax.imshow(
            np.swapaxes(displ_1[:, :, ind_slice], 0, 1),
            vmin=vmin_1,
            vmax=vmax_1,
            cmap=cmap_1,
            interpolation="nearest",
            origin="lower",
            extent=[0, phys_sp[0], 0, phys_sp[1]],
        )

        if displ_2 is not None:
            data = np.swapaxes(displ_2[:, :, ind_slice], 0, 1)
            im2 = ax.imshow(
                data,
                vmin=vmin_2,
                vmax=vmax_2,
                cmap=cmap_2,
                alpha=np.where(data < vmin_2, 0, 0.5),
                interpolation="nearest",
                origin="lower",
                extent=[0, phys_sp[0], 0, phys_sp[1]],
            )

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(False)
        ax.axis("off")
        bgcolor = cmap_1(min(vmin_1, 0.0))
        fgcolor = cmap_1(vmax_1)

        if slice_numb + 1 > len(ind_slices):
            ind_slice = ""

        ax.text(
            0.98,
            0.01,
            "{}".format(ind_slice),
            color=fgcolor,
            transform=ax.transAxes,
            horizontalalignment="right",
            verticalalignment="bottom",
            fontsize=fontsize,
            bbox=dict(
                boxstyle="square,pad=0", edgecolor=bgcolor, facecolor=bgcolor
            ),
        )

    if brain_data_2 is not None:
        cbar = grid.cbar_axes[0].colorbar(
            im2, orientation="horizontal", label="Intensity"
        )
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label("Intensity", size=8)

    if data_2 is None:
        fname, _ = os.path.splitext(os.path.basename(data_1))

    else:
        fname, _ = os.path.splitext(os.path.basename(data_2))

    if out_name:
        fname += "_" + out_name

    if out_dir is None:
        out_dir = os.path.dirname(data_2)

    out_file = os.path.join(out_dir, fname + "_slice_planes_plot.png")
    counter = 1

    while os.path.exists(out_file):
        b_name, ext = os.path.splitext(out_file)
        out_file = f"{b_name}_{counter:03d}{ext}"
        counter += 1

    fig.savefig(out_file, format="png", dpi=300, bbox_inches="tight")
    return out_file
