import sys

# Current version
version_major = 1
version_minor = 0
version_micro = 3
version_extra = ""

# Expected by setup.py: string of form "X.Y.Z"
__version__ = "{0}.{1}.{2}".format(version_major, version_minor, version_micro)


# Expected by setup.py: the status of the project
CLASSIFIERS = ['Development Status :: 5 - Production/Stable',
               'Intended Audience :: Developers',
               'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
               'Topic :: Software Development :: Libraries :: Python Modules',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
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
]
EXTRA_REQUIRES = {
    'doc': [
        'sphinx>=1.0',
    ],
}

brainvisa_build_model = 'pure_python'

# tests to run
test_commands = ['%s -m mia_processes.test' % sys.executable]
