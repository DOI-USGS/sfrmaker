============
Installation
============

``SFRmaker`` depends on a number of python packages, many of which have external C library dependencies. The easiest way to install most of these is with `conda`_. A few packages are not available via conda, and must be installed with `pip`_. Instructions for setting up an appropriate `conda environment`_, including considerations for users on the USGS internal network, are available |conda_instructions|.

Installing SFRmaker from the `python package index`_ using pip
--------------------------------------------------------------
At the command line:

.. code-block:: bash

    pip install sfrmaker

Installing SFRmaker from source
--------------------------------------------------------------
Clone or `download`_ the `sfrmaker repository`_. Then with a command window pointed at the root level (containing ``setup.py``):

.. code-block:: bash

    pip install -e .

This installs SFRmaker to the active python distribution by linking the source code in-place, making it easier to periodically pull updates using git.

Updating SFRMaker using Git
--------------------------------
To update SFRmaker (if it was cloned using Git), at the root level of the repository:

.. code-block:: bash

    git pull origin master

Alternatively, SFRmaker could be updated by downloading the repository again and installing via ``pip install -e .``.


.. _conda: https://docs.conda.io/en/latest/
.. _conda environment: https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html
.. |conda_instructions| raw:: html

   <a href="https://aleaf.github.io/modflow-setup/installation.html" target="_blank">here</a>

.. _download: https://github.com/aleaf/sfrmaker/archive/master.zip
.. _pip: https://packaging.python.org/tutorials/installing-packages/#use-pip-for-installing
.. _python package index: https://pypi.org/

.. _sfrmaker repository: https://github.com/aleaf/SFRmaker
