# mia_processes

<!-- [![](https://codecov.io/github/populse/mia_processes/coverage.svg?branch=master)](https://codecov.io/github/populse/mia_processes) -->
[![](https://img.shields.io/badge/license-CeCILL-blue.svg)](https://github.com/populse/mia_processes/blob/master/LICENSE)
[![](https://img.shields.io/pypi/v/mia_processes.svg)](https://pypi.org/project/mia_processes/)
[![](https://img.shields.io/badge/python-3.5%2C%203.6%2C%203.7-yellow.svg)](#)
[![](https://img.shields.io/badge/platform-Linux%2C%20OSX%2C%20Windows-orange.svg)](#)

# Documentation

[The documentation is available in mia_processes's website](https://populse.github.io/mia_processes).

# Installation

* A compatible version of [Python](https://www.python.org/) (>= 3.5) and [pip](https://packaging.python.org/guides/tool-recommendations/) must be installed.

* Depending on the Python setup and OS, the “python3” command can be use as the default Python command.

    Try:
        `python -V`
	
    If it returns `Python 3.x.x`, replace all the `python3` commands below by `python`.
    
    If not installed, install it ...
 
 * Make sure to have pip installed.
 
    Try:
        `pip -V`

    If it returns a path including `/python3.x/`, replace all the `pip3` commands below by `pip`.
    
    If not installed, install it ...
    
* From PyPI

  * [Install the latest version of mia_processes and its dependencies from the Python Packaging Index](https://docs.python.org/3/installing/index.html):
  
        pip3 install mia_processes # depending of the setup, it could be necessary to add --user option

* From source, for Linux distributions.

  * Install a Version Control System, for example [git](https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control). Depending of your distribution, [package management system](https://en.wikipedia.org/wiki/Package_manager) can be different:
  
        sudo apt-get install git # Debian like
        sudo dnf install git # Fedora 22 and later
        # etc.
	
  * Clone the source codes.

    * Get source codes from Github. Replace [mia_processes_install_dir] with a directory of your choice:

          git clone https://github.com/populse/mia_processes.git [mia_processes_install_dir]

    * Or download the zip file (mia_processes-master.zip) of the project ([green button "Clone or download"](https://github.com/populse/mia_processes)), then extract the data in the directory of your choice ([mia_processes_install_dir]):

          unzip mia_processes-master.zip -d [mia_processes_install_dir]  # In this case [mia_processes_install_dir] becomes [mia_processes_install_dir]/mia_processes-master
	
  * Install the Python module distribution:

        cd [mia_processes_install_dir]  
        python3 setup.py install # depending of the setup, it could be necessary to add --user option

  * Remove the [mia_processes_install_dir] directory:

        cd ..  
        rm -r [mia_install_dir]  

# License

* The whole populse project is open source.
* Mia_processes is precisely released under the CeCILL software license.
* All license details can be found [here](http://cecill.info/licences/Licence_CeCILL_V2.1-en.html), or refer to the license file [here](https://github.com/populse/mia_processes/blob/master/LICENSE).

# Support and Communication

In case of a problem or to ask a question about how to do something in mia_processes, please [open an issue](https://github.com/populse/mia_processes/issues).

The developer team can even be contacted using populse-support@univ-grenoble-alpes.fr.
