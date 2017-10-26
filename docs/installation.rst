Installation
============

cMonkey\ :sub:`2`
-------------------
SYGNAL is designed to be used a follow up to a cMonkey\ :sub:`2` analysis for discovery of co-regulated biclusters (https://github.com/baliga-lab/cmonkey2). More information about cMonkey\ :sub:`2` can be found in the publication (https://www.ncbi.nlm.nih.gov/pubmed/25873626) and through the cMonkey\ :sub:`2` Wiki (https://github.com/baliga-lab/cmonkey2/wiki).

SYGNAL dependencies
-------------------

Python
~~~~~~

SYGNAL uses the multiprocessing capability of Python to speed up all the different computations required, and as such makes use of the ``multiprocessing`` module which is available in Python 2.7 and greater (https://www.python.org/downloads/). SYGNAL has been developed and tested to work with Python 2.7.x and Python 3.x.

R
~~~

In addition SYGNAL requires the installation of ``rpy2`` which can be found at (http://rpy2.bitbucket.org/). Each version of ``rpy2`` requires the installation of a specific version of the statistical package R. Please install the corresponding version of R (https://cran.r-project.org/).

Here are easy commands to install R and ``rpy`` on Ubuntu:

.. highlight:: none

::

   sudo apt-get install r-base r-base-dev
   sudo pip install rpy

R packages
~~~~~~~~~~

  * ``WGCNA``
  * ``impute``
  * ``getopt``
  * ``topGO``
  * ``org.Hs.eg.db`` or ``org.Mm.eg.db`` depending on species
  * ``GOSim``

Make sure to set your repositories to include the Bioconductor repositories with setRepositories() and install the apporiate packages through the commands below (you will also have to select the closest mirror as well):

.. highlight:: none

::

   setRepositories()
   install.packages(c('WGCNA','impute','getopt','topGO','org.Mm.eg.db','org.Hs.eg.db','GOSim'))

A host of other dependencies for installation of each package should be installed through this command. If there are failures in the installation of these packages please try to install each dependency separately, and contact the authors of these packages if problems persist.

In order to run the causality analysis portions of SYGNAL requires the use of the Network Edge Orienting (``NEO``) package   (https://labs.genetics.ucla.edu/horvath/aten/NEO/).  The ``neoDecember2015.txt`` is the source code needed for running ``NEO`` in SYGNAL.


Stand alone applications
~~~~~~~~~~~~~~~~~~~~~~~~

SYGNAL relies upon several sequence searching based algorithms ``MEME``, ``WEEDER`` and ``AME``. They will have to installed according to their documentation. Please find them at the following links:

  * MEME (http://meme-suite.org/doc/download.html?man_type=web)
  * WEEDER (https://github.com/baliga-lab/weeder_patched)
  * AME (http://bioinformatics.org.au/tools/ame/)

Installing SYGNAL
-----------------
As SYGNAL is a Python program it simply requires the downloading  the code from the github repository (https://github.com/baliga-lab/sygnal/archive/master.zip) into the directory you wish to run SYGNAL. As long as the dependencies are met then it should be possible to run SYGNAL. We have developed a best practices directory structure for using SYGNAL. Please refer to the section
:doc:`Preparing for SYGNAL run <preparing>`.
