# mia_processes

[![](https://codecov.io/github/populse/mia_processes/coverage.svg?branch=master)](https://codecov.io/github/populse/mia_processes)
[![](https://img.shields.io/badge/license-CeCILL-blue.svg)](https://github.com/populse/mia_processes/blob/master/LICENSE)
[![](https://img.shields.io/pypi/v/mia_processes.svg)](https://pypi.org/project/mia_processes/)
[![](https://img.shields.io/badge/python-3.5%2C%203.6%2C%203.7-yellow.svg)](#)
[![](https://img.shields.io/badge/platform-Linux%2C%20OSX%2C%20Windows-orange.svg)](#)

# Documentation

The documentation is available on mia_processes's website here: [https://populse.github.io/mia_processes](https://populse.github.io/mia_processes)

# Installation

 * Note: depending on your Python setup and OS, the “python3” command can be use as the default Python command.
     Try:
 
       python -V
       
     If it returns `Python 3.x.x`, replace all the `python3` commands below by `python`.

* From PyPI

  * Make sure to have pip3 installed. You can verify it by typing the following in a command line:
  
        pip3 --version
  
  * [Install the latest version of mia_processes and its dependencies from the Python Packaging Index:](https://docs.python.org/3/installing/index.html)
  
        python3 -m pip install mia_processes

* From source, for Linux distributions

  * A compatible version of Python must be installed
  
  * Install a Version Control System, for example [git](https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control). Depending of your distribution, [package management system](https://en.wikipedia.org/wiki/Package_manager) can be different
  
        sudo apt-get install git # Debian like
        sudo dnf install git # Fedora 22 and later
        # etc.
	
  * Clone the source codes

    * Get source codes from Github. Replace [mia_processes_install_dir] with a directory of your choice

          git clone https://github.com/populse/mia_processes.git [mia_processes_install_dir]

    * Or download the zip file (mia_processes-master.zip) of the project ([green button "Clone or download"](https://github.com/populse/mia_processes)), then extract the data in the directory of your choice [mia_processes_install_dir]

          unzip mia_processes-master.zip -d [mia_processes_install_dir]  # In this case [mia_processes_install_dir] becomes [mia_processes_install_dir]/mia_processes-master
	
  * Install the Python module distribution

        cd [mia_processes_install_dir]  
        python3 setup.py install --user # Ensure that you use python >= 3.5 (use python3.x to be sure)  

  * Remove the [mia_processes_install_dir] directory:

        cd ..  
        rm -r [mia_install_dir]  

# License

* The whole populse project is open source
* mia_processes is precisely released under the CeCILL software license
* You can find all the details on the license [here](http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html), or refer to the license file [here](https://github.com/populse/mia_processes/blob/master/LICENSE)

# Support and Communication

If you have a problem or would like to ask a question about how to do something in mia_processes, please [open an issue](https://github.com/populse/mia_processes/issues).

You can even contact the developer team by using populse-support@univ-grenoble-alpes.fr.
