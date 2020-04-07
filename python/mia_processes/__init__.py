# -*- coding: utf-8 -*- #

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .info import __version__

_doc_path = None

def _init_doc_path():
    global _doc_path
    import os
    from .info import version_major
    from .info import version_minor
    import mia_processes

    opd = os.path.dirname
    p = os.path.join(
        opd(opd(opd(mia_processes.__file__))),
        'docs/html/process_docs/mia_processes')
    if os.path.exists(p):
        _doc_path = p
        return _doc_path
    _doc_path = 'https://populse.github.io/mia_processes/process_docs/mia_processes'
    return _doc_path

_init_doc_path()


