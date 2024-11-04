============
Installation
============

Installing python dependencies with Conda
-----------------------------------------
SFRmaker depends on a number of python packages, many of which have external C library dependencies. The easiest way to install most of these is with a package installer like `Conda`_. A few packages are not available via conda, and must be installed with `pip`_. If you are on the USGS internal network, see the `Considerations for USGS Users`_ section below first.

Download and install a python distribution and Conda-like package manager 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are many ways to do this:

    * The `Anaconda python distribution`_ comes with a large selection of popular data science and scientific packages pre-installed.

    * `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ is a minimal installer with a much smaller footprint, making it ideal for creating python environments dedicated to specific tasks (a recommended practice).

    * `Mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`_ is like Miniconda, but pre-configured to use the `Mamba`_ installer, and only the `conda-forge <https://conda-forge.org/docs/user/introduction.html>`_ channel for getting packages (more below). If the above two options don't work (for example, the Conda installer fails or gets stuck on the "solve" step), this may be your best option.

**Make sure to install the python distribution to your username** (not at the system level). More often than not, installing at the system level (for all users) seems to result in issues with library dependencies (for example, import of ``fiona`` or ``rasterio`` failing because gdal isn't found). It is also good practice to periodically do a `clean uninstall`_ of Anaconda, which at the system level requires admin. privileges.

Mambaforge-specific instructions:
  * If after downloading the Mambaforge installer you get a message that it "isn't commonly downloaded" and the option to "Trust" it is grayed out (because you don't have sufficient rights), try using a different browser like Google Chrome (for example, if you were using Microsoft Edge).
  * Once the installer is running, go against the recommendation and **select the option to add Mambaforge to the system path for your user.**

Anaconda-specific instructions:
  * In the installer, at the “Destination Select” step, select “Install for me only.” It should say something about how the software will be installed to your home folder.

  * If your installer skips the “Destination Select” step, when you get to "Installation Type", click “Change Install Location” and then “Install for me only.”

Download an environment file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  * `requirements.yml`_ for a `conda environment`_ with the minimum packages required to run SFRmaker, or

  * `gis.yml`_ for a more full set of packages in the python geospatial stack, including Jupyter Notebooks and packages needed for testing, documentation and packaging. Note that the environment described by ``requirements.yml`` is called `sfrmaker`, while the environment in ``gis.yml`` is called `gis`.

    .. note::
        To download the above YAML files, simply follow the links to get the raw text and then go to File > Save within your web browser, and save the text as a YAML file (with the `.yaml` or `.yml` extension).

  * Alternatively, clone (`using git`_) or `download`_ the SFRmaker repository, which includes the two environment files at the root level.

  * Note that both of these environment files contain a ``pip`` section of packages that will be installed with pip, after the Conda packages are installed.

Creating a `Conda environment`_ using `Mamba`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you are on the USGS internal network, see the `Considerations for USGS Users`_ section below first.

While Conda can of course be used to create a Conda environment, the Mamba package solver is generally faster and more robust, especially for larger, more complex environments like the included ``requirements.yml``. Mamba is a reimplementation of the conda package manager in C++.

Before using Mamba, you will need to `install it first <https://mamba.readthedocs.io/en/latest/installation.html>`_.

Python packages are available from conda via channels. Conda comes preconfigured to install packages from the default channel, which is maintained by Anaconda, Inc. In general, you may have better luck exclusively using the `conda-forge <https://conda-forge.org/docs/user/introduction.html>`_ channel instead, which is community-based and intended to provide a single location to get any package, with a minimum of hassle. In general, it is bad practice to mix package channels within a single environment. You can read more `here <https://conda-forge.org/docs/user/introduction.html>`__, but to set conda-forge as the default:

.. code-block:: bash

    conda config --add channels conda-forge

.. note::
    If you are having trouble installing Mamba (for example, the conda package solver fails when you try to install it, or takes an excessively long time), you may have better luck uninstalling `Anaconda completely <clean uninstall>`_ and installing `Mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`_ instead, as directed in the Mamba install instructions. Mambaforge solves both of the above problems by providing a minimal python distribution and conda-style package installer that is preconfigured to use both `conda-forge <https://conda-forge.org/docs/user/introduction.html>`_ and `Mamba`_.

Once you have a python distribution and mamba installed, to create the conda environment, open a new Anaconda Command Prompt on Windows or a new terminal window on OSX and point it to the location of ``requirements.yml`` or ``gis.yml`` and enter:

.. code-block:: bash

    mamba env create -f requirements.yml


    .. note::
        Creating the ``requirements.yml`` environment (or any environment with ``git+https: ...`` installs) requires Git to be installed and visible in the system path where ``env create`` is being run. If Git is installed and somehow not in the system path, it can be added to the system path on Windows 10 without admin. rights via the "environment variables" editor under User Accounts in the Control Panel (Google it).

Building the environment will probably take a while. If the build fails because of an SSL error, fix the problem (see `Considerations for USGS Users`_ below) and either:


    a) 	Update the environment

        .. code-block:: bash

            mamba env update -f requirements.yml

    b) 	or remove and reinstall it:

        .. code-block:: bash

            conda env remove -n sfrmaker
            mamba env create -f requirements.yml

Keeping the Conda environment up to date
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The python packages and other open source software libraries that SFRmaker depends on are continually changing. SFRmaker aims to mostly follow the `Numpy guidelines for package support <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_, which effectively means that the two latest minor versions of Python (e.g. 3.11 and 3.10) and their associated Numpy versions will be supported. However, occasionally backwards compatability with a particular package may be broken in a shorter timeframe, in which case the minimum required version of that package will be specified in the ``requirements.yml`` file. All of this to say that your Conda environment will eventually get out of date. The `Conda documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ has instructions for updating packages within a Conda environment, but at some point (perhaps a few times a year) it is good practice to simply delete the environment and rebuild it from the `.yml` file. Every so often, you may also want to reinstall Anaconda after a `clean uninstall`_.

Installing SFRmaker
------------------------
There are several ways to install SFRmaker. Regardless of the method, the installation must be performed in a python
environment with the required dependencies. In the case of the Conda environment created above, the environment must be activated, so that right version of python is called when ``python`` is entered at the command line:

.. code-block:: bash

    conda activate sfrmaker

Installing and updating SFRmaker from `PyPI <https://pypi.org/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once a suitable conda environment (that contains ALL of the dependencies) is made and activated, the simplest way to install SFRmaker is from the Python Package Index using pip.

.. code-block:: bash

    pip install sfrmaker

Subsequent releases of SFRmaker to PyPI can then be installed with

.. code-block:: bash

    pip install --upgrade sfrmaker

Note that in some situations you may have to ``pip uninstall sfrmaker`` and then ``pip install sfrmaker``. You can always check
what version of sfrmaker you have within a python session with

.. code-block:: python

    import sfrmaker
    sfrmaker.__version__

Or if you are using Conda, at the command line with

.. code-block:: bash

    conda list

Installing the latest develop version of SFRmaker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In some situations you may want the bleeding-edge version of SFRmaker that is being actively developed on GitHub. For example,
to incorporate a bug fix that was made after the latest release. Pip can also be used to fetch SFRmaker directly from GitHub:

.. code-block:: bash

    pip install git+https://github.com/doi-usgs/sfrmaker@develop

(for the develop branch). Subsequent updates can then be made with

.. code-block:: bash

    pip uninstall sfrmaker
    pip install git+https://github.com/doi-usgs/sfrmaker@develop

Installing the SFRmaker source code in-place
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, if you intend to contribute to SFRmaker (please do!) or update your install frequently, the best route is probably to clone the source code from git and install it in place.

.. code-block:: bash

    git clone https://github.com/doi-usgs/sfrmaker.git
    cd sfrmaker
    pip install -e .

.. note::
    Don't forget the ``.`` after ``pip install -e``!

Your local copy of the SFRmaker repository can then be subsequently updated with

.. code-block:: bash

    git pull origin master

.. note::
    If you are making local changes to SFRmaker that you want to contribute, the workflow is slightly different. See the :ref:`Contributing to SFRmaker` page for more details.


The advantage of installing the source code in-place is that any changes you make are automatically incorporated into your python environment, without any additional install commands. When debugging in an interactive development environment (IDE) such as Pycharm or VS Code, error tracebacks and inspection features go to the actual source code, not the version installed in the ``site-packages`` folder. Additionally, since this install is done through pip, ``pip uninstall``
will work to remove the package, and the current version of the package (including the latest commit information) will be visible with ``conda list``.

Installing the IPython kernel to use SFRmaker in Jupyter Notebooks
------------------------------------------------------------------------------------------------
This step may not be needed if you already have an existing Python environment with the packages required by SFRmaker *and* Jupyter Notebook installed. However, if you'd like to use SFRmaker in a Jupyter Notebook with the included ``sfrmaker`` environment (specified in ``requirements.yml``), you'll most likely need to install the IPython kernel in that environment. You can do this at the command line or terminal window (with ``sfrmaker`` activated):

.. code-block:: bash

    python -m ipykernel install --user --name sfrmaker --display-name "sfrmaker"


The first instance of ``sfrmaker`` in this command is the environment to install the kernel to, and the second instance (in quotes) is the name that will appear in the ``Kernel`` menu within Jupyter Notebook. To use the kernel, simply select it from the ``Kernel > Change kernel`` menu within  Jupyter Notebook.

Best practices
------------------------

* Install the \*conda distribution of your choice to your user account, NOT at the system level. Installing to your user means you have rights to delete and reinstall Anaconda as-needed, as well as to edit any configuration files for pip, Conda, etc. Installing at the system level also just seems to lead to more confusing problems with dependencies, at least in the USGS.
* Periodically (maybe a few times a year?) fully remove your \*conda distribution and reinstall it. If you just can't get things to work (packages won't import or produce DLL errors on import, adding or upgrading a package takes a very long time or results in excessive upgrades or downgrades of other packages, etc.), fully removing and reinstalling \*conda just may resolve your issues.
* Don't use your base environment; create and delete environments as needed. Conda is generally pretty good about managing packages between environments without wasting a lot of disk space.
* Use an environment file (as above) to create a conda environment, instead of installing packages ad-hoc.
* Use Mamba instead of Conda; it just works better for environments with a lot of packages.
* After setting up the above conda environment, scan the screen output to make sure that everything installed correctly, especially the packages installed through pip.
* Avoid mixing package channels within a Conda environment. Strictly sticking to conda-forge may yield the best results.
* Use `conda-pack`_, rather than an overly-detailed environment file, to guarantee reproducibility.


_`Considerations for USGS Users` (or anyone else with an organizational SSL certificate)
-------------------------------------------------------------------------------------------
Using conda or pip on the USGS network requires SSL verification, which can cause a number of issues.
If you are encountering persistant issues with creating the conda environment,
you may have better luck trying the install off of the USGS network (e.g. at home).
See `here <https://tst.usgs.gov/applications/application-and-script-signing/>`__ for more information
about SSL verification on the USGS network, and to download the DOI SSL certificate. It might be most helpful to consult with your local IT staff; some centers pre-install the appropriate certificate bundle on all new computers.

_`Installing the DOI SSL certificate for use with pip`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1) `Download the DOI SSL certificate (internal DOI network only) <https://tst.usgs.gov/applications/application-and-script-signing/>`_
2) Create a pip configuration file at your user level as described `here <https://pip.pypa.io/en/stable/topics/configuration/>`_
    * On Windows, this may be ``C:\Users\<your username>\AppData\Roaming\pip\pip.ini``.
    * On OSX, this may be ``/Users/<your username>/Library/Application Support/pip/pip.conf``.

::

    [global]
    cert = <path to DOI certificate bundle (e.g. cacert_usgs.pem)>


`More information on the pip configuration file and it's location <https://pip.pypa.io/en/stable/topics/configuration/>`_

You can verify the pip configuration that is being used at a particular location with

::  

    pip config list


.. note::
    When you are off the USGS network, you may have to comment out the ``cert=`` line with a leading ``#`` 
    in the above pip configuration file to get ``pip`` to work.

.. note::
    To test that pip is working, simply try ``pip install --upgrade`` on a package of choice; ideally one that is relatively downstream in the dependency stack (for example, ``gis-utils``), so as to minimize any additional changes to your environment that may occur.


_`Installing the DOI SSL certificate for use with git`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using git to reach out to the web also requires the appropriate certificate bundle. This is stored in the git configuration under the ``http.sslCAinfo`` variable. A valid certificate is needed here to clone and push/pull from git repositories, and also to install packages from git using pip (for example, the entries in ``requirements.yml`` starting with ``git+https://github.com/``). A missing certificate bundle will cause setup of the ``sfrmaker`` environment to fail on the pip installs section, even if a valid certificate bundle is specified in the pip configuration bundle. To add your certificate bundle to git:

.. code-block:: bash

    git config --global http.sslCAinfo <path to DOI certificate bundle (e.g. cacert_usgs.pem)>

You can then verify that the certificate was added by displaying your git configuration

.. code-block:: bash

    git config --list

.. note::
    To test that both pip and git are working, simply try ``pip install --upgrade`` on a package in the pip section of ``requirements.yml`` that starts with ``git+https://``; for example ``git+https://github.com/modflowpy/flopy@develop``.


Installing the DOI SSL certificate for use with conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See `these instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#ssl-verification-ssl-verify>`_.
This may or may not work. Basically, ``ssl_verify:`` needs to be set in your `condarc`_ file to point
to a valid SSL certificate, which may be different from the basic ``DOIRootCA2.cer`` file.

You can find the location of your `condarc`_ file with::

    conda info -a

which displays information about how Conda is configured. Note that you may have multiple `condarc`_
files at the system, user and possibly project levels.

Common issues:

* Conda Install fails on the USGS network without a certificate, or with an incorrectly formatted certificate.
  Possible solutions are to either get a correctly formatted SSL certificate from your IT person, or try installing off the network.
* Conda Install fails off the USGS network with a certificate (may or may not be correctly formatted). Solution:
  open your `condarc`_ file
  and comment out the SSL certificate file, if it is specified. E.g.::

    ssl_verify: #D:\certificates\DOIRootCA2.cer



Troubleshooting issues with the USGS network
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::
    If all else fails, you can simply try running the installs off-network.


SSL-related error messages when using conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(with ``SSL`` mentioned in the message and possibly ``bad handshake``)

Make sure that the ``Conda`` package installer is configured to use the USGS certificate
(see :ref:`Installing the DOI SSL certificate for use with conda` above).


SSL-related error messages when using building a conda environment or using pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If the error message looks like this:

.. code::

    Could not fetch URL https://pypi.org/simple/pip/: 
    There was a problem confirming the ssl certificate: HTTPSConnectionPool(host='pypi.org', port=443): 
    Max retries exceeded with url: /simple/pip/ (Caused by SSLError(SSLCertVerificationError(1, '[SSL: CERTIFICATE_VERIFY_FAILED] certificate 
    verify failed: unable to get local issuer certificate (_ssl.c:997)'))) - skipping 

Make sure that the ``pip`` package installer is configured to use the USGS certificate bundle
(see `Installing the DOI SSL certificate for use with pip`_ above).

If the error message looks like this:

.. code::

    Pip subprocess error:
    Running command git clone --filter=blob:none --quiet https://github.com/<some package>
    ...
    fatal: unable to access 'https://github.com/<some package>': SSL certificate problem: unable to get local issuer certificate
    error: subprocess-exited-with-error
    ...
    note: This error originates from a subprocess, and is likely not a problem with pip.

Make sure that git is configured to use the USGS certificate bundle
(see `Installing the DOI SSL certificate for use with git`_ above).


If you are on the USGS network, using Windows, and you get this error message:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
..

    CondaHTTPError: HTTP 500 INTERNAL ERROR for url <https://repo.anaconda.com/pkgs/msys2/win-64/m2w64-gettext-0.19.7-2.tar.bz2>
    Elapsed: 00:30.647993

    An HTTP error occurred when trying to retrieve this URL.
    HTTP errors are often intermittent, and a simple retry will get you on your way.

Adding the following line to ``environment.yml`` should work:

.. code-block:: yaml

    - msys2::m2w64-gettext


This tells conda to fetch ``m2w64-gettext`` from the ``msys2`` channel instead. Note that this is only a dependency on Windows,
so it needs to be commented out on other operating systems (normally it wouldn't need to be listed, but the above HTTP 500 error indicates that installation from the default source location failed.)


.. _Anaconda python distribution: https://www.anaconda.com/distribution/
.. _clean uninstall: https://docs.anaconda.com/anaconda/install/uninstall/
.. _Conda: https://docs.conda.io/en/latest/
.. _Mamba: https://mamba.readthedocs.io/en/latest/
.. _conda environment: https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html
.. _conda-pack: https://conda.github.io/conda-pack/
.. _condarc: https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html
.. _download: https://github.com/aleaf/sfrmaker/archive/master.zip
.. _gis.yml: https://raw.githubusercontent.com/aleaf/sfrmaker/master/gis.yml
.. _pip: https://packaging.python.org/tutorials/installing-packages/#use-pip-for-installing
.. _Readme file: https://github.com/aleaf/sfrmaker/blob/master/Readme.md
.. _requirements.yml: https://raw.githubusercontent.com/aleaf/sfrmaker/master/requirements.yml
.. _using git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git


.. _sfrmaker repository: https://github.com/aleaf/SFRmaker
