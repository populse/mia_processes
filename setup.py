# -*- coding: utf-8 -*- #

"""The module dedicated to install the mia_processes library."""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os
from setuptools import find_packages, setup

# Select modules to include in distribution
modules = find_packages('python')
print('\nmodules found: ', modules)

# Additional script files to include in distribution
scripts = []

# Data files to include in distribution
pkgdata = {
    # format is:
    # <python module>: <file extensions list>
}

# Read the info.py file in mia_processes module
release_info = {}
python_dir = os.path.join(os.path.dirname(__file__), 'python')

with open(os.path.join(python_dir, 'mia_processes', 'info.py')) as f:
    code = f.read()
    exec(code, release_info)


# Build the setup
setup(
    name=release_info["NAME"],
    description=release_info["DESCRIPTION"],
    long_description=release_info["LONG_DESCRIPTION"],
    license=release_info["LICENSE"],
    classifiers=release_info["CLASSIFIERS"],
    author=release_info["AUTHOR"],
    author_email=release_info["AUTHOR_EMAIL"],
    version=release_info["VERSION"],
    url=release_info["URL"],
    package_dir={'': 'python'},
    packages=modules,
    package_data=pkgdata,
    platforms=release_info["PLATFORMS"],
    install_requires=release_info["REQUIRES"],
    extras_require=release_info["EXTRA_REQUIRES"],
    scripts=scripts,
    zip_safe=False
)
