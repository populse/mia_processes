# -*- coding: utf-8 -*- #

"""The module dedicated to the main informations on mia_processes.

The info.py module is mainly used by the setup.py module.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os
import subprocess
import sys

# Current version
version_major = 2
version_minor = 0
version_micro = 1
version_extra = "" # leave empty for release

# Expected by setup.py: string of form "X.Y.Z"
if version_extra:
    __version__ = "{0}.{1}.{2}-{3}".format(version_major, version_minor, version_micro, version_extra)
    
else:
    __version__ = "{0}.{1}.{2}".format(version_major, version_minor, version_micro)
   
def get_gitversion():
    """Mia_processes version as reported by the last commit in git
    Returns the version or None if nothing was found
    """
    try:
        import mia_processes
        dir_mia_processes = os.path.realpath(
            os.path.join(os.path.dirname(mia_processes.__file__), os.path.pardir, os.path.pardir))

    except:
        dir_mia_processes = os.getcwd()

    dir_mia_processesgit = os.path.join(dir_mia_processes, ".git")

    if not os.path.exists(dir_mia_processesgit):
        return None

    ver = None

    try:
        gitversion, _ = subprocess.Popen(
            "git show -s --format=%h", shell=True, cwd=dir_mia_processes, stdout=subprocess.PIPE
        ).communicate()

    except Exception:
        pass

    else:
        ver = gitversion.decode().strip().split("-")[-1]

    return ver

if __version__.endswith("-dev"):
    gitversion = get_gitversion()

    if gitversion:
        __version__ = "{0}+{1}".format(__version__, gitversion)

# Expected by setup.py: the status of the project
CLASSIFIERS = ['Development Status :: 5 - Production/Stable',
               'Intended Audience :: Developers',
               'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre '
                 'License, version 2.1 (CeCILL-2.1)',
               'Topic :: Software Development :: Libraries :: Python Modules',
               'Operating System :: OS Independent',
               'Programming Language :: Python :: 3.5',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Python :: 3.7',
               'Programming Language :: Python :: 3.8',
               'Programming Language :: Python :: 3 :: Only',
               'Topic :: Scientific/Engineering',
               'Topic :: Utilities']

# project descriptions
DESCRIPTION = 'mia_processes'
LONG_DESCRIPTION = """
===============
mia_processes
===============

The pipelines library for the MIA [Multiparametric Image
Analysis, under the POPULSE umbrella] software
"""

# Other values used in setup.py
NAME = 'mia_processes'
ORGANISATION = 'populse'
MAINTAINER = 'Populse team'
MAINTAINER_EMAIL = 'populse-support@univ-grenoble-alpes.fr'
AUTHOR = 'Populse team'
AUTHOR_EMAIL = 'populse-support@univ-grenoble-alpes.fr'
URL = 'http://populse.github.io/mia_processes'
DOWNLOAD_URL = 'http://populse.github.io/mia_processes'
LICENSE = 'CeCILL'
VERSION = __version__
CLASSIFIERS = CLASSIFIERS
PLATFORMS = 'OS Independent'
REQUIRES = [
    'populse_mia==2.*',
    'populse_db>=2.0.0'
]
EXTRA_REQUIRES = {
    'doc': [
        'sphinx>=1.0',
    ],
}

brainvisa_build_model = 'pure_python'

# tests to run
# test_commands = ['%s -m mia_processes.test' % sys.executable]
