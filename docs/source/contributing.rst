Contributing to SFRmaker
=============================

(Note: much of this page was cribbed from the `geopandas <https://geopandas.org/>`_ project,
which has similar guidelines to `pandas <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_
and `xarray <http://xarray.pydata.org/en/stable/>`_.)

Getting started
----------------
All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome. If an issue that interests you isn't already listed in the `Issues tab`_, consider `filing an issue`_.

Bug reports and enhancement requests
------------------------------------------------
Bug reports are an important part of improving SFRmaker. Having a complete bug report
will allow others to reproduce the bug and provide insight into fixing. See
`this stackoverflow article <https://stackoverflow.com/help/mcve>`_ and
`this blogpost <https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>`_
for tips on writing a good bug report.

Before doing anything else, check the :ref:`Troubleshooting` page for solutions to common issues.

Trying the bug-producing code out on the *develop* branch is often a worthwhile exercise
to confirm the bug still exists. It is also worth searching existing bug reports and pull requests to see if the issue has already been reported and/or fixed.

To file a bug report or enhancement request, from the issues tab on the `SFRmaker GitHub page <https://github.com/usgs/sfrmaker>`_, select "New Issue".

Bug reports must:

#. Include a short, self-contained Python snippet reproducing the problem, along with the contents of your configuration file (if you are using one) and the full error traceback.
   You can format the code nicely by using `GitHub Flavored Markdown
   <https://github.github.com/github-flavored-markdown/>`_::

      ```python
      >>> import sfrmaker
      >>> lns = sfrmaker.Lines.from_nhdplus_v2(NHDPlus_paths='../tylerforks/NHDPlus/',
                                               filter='../tylerforks/grid.shp')
      ...
      ```

   e.g.::

      ```yaml
      <paste configuration file contents here>
      ```

      ```python
      <paste error traceback here>
      ```

#. Include the version of SFRmaker that you are running, which can be obtained with:

   .. code-block:: python

       import sfrmaker
       sfrmaker.__version__

   Depending on the issue, it may also be helpful to include information about the version
   of python and operating system.

#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then be visible on the `Issues tab`_ and open to comments/ideas from others.


Code contributions
------------------------------
Code contributions to SFRmaker to fix bugs, implement new features or improve existing code are encouraged. Regardless of the context, consider `filing an issue`_ first to make others aware of the problem and allow for discussion on potential approaches to addressing it.

In general, SFRmaker trys to follow the conventions of the pandas project where applicable. Contributions to SFRmaker are likely to
be accepted more quickly if they follow these guidelines.

In particular, when submitting a pull request:

- All existing tests should pass.  Please make sure that the test
  suite passes, both locally and on
  `GitHub Actions <https://github.com/usgs/sfrmaker/actions>`_.  Status on
  the GitHub Actions and code coverage checks will be visible on a pull request.

- New functionality should include tests.  Please write reasonable
  tests for your code and make sure that they pass on your pull request.

- Classes, methods, functions, etc. should have docstrings.  The first
  line of a docstring should be a standalone summary.  Parameters and
  return values should be documented explicitly. (Note: there are admittedly more than a few places in the existing code where docstrings are missing. Docstring contributions are especially welcome!

- Follow PEP 8 when possible. For more details see
  :ref:`below <contributing_style>`.

- Following the `FloPy Commit Message Guidelines <https://github.com/modflowpy/flopy/blob/develop/CONTRIBUTING.md>`_ (which are similar to the `Conventional Commits <https://www.conventionalcommits.org/en/v1.0.0/>`_ specification) is encouraged. Structured commit messages like these can result in more explicit commit messages that are more informative, and also facilitate automation of project maintenance tasks.

- Imports should be grouped with standard library imports first,
  3rd-party libraries next, and SFRmaker imports third.  Within each
  grouping, imports should be alphabetized.  Always use absolute
  imports when possible, and explicit relative imports for local
  imports when necessary in tests. Imports can be sorted automatically using the isort package with a pre-commit hook. For more details see :ref:`below <contributing_style>`.

- SFRmaker supports Python 3.9+ only.


Seven Steps for Contributing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are seven basic steps to contributing to *SFRmaker*:

1) Fork the *SFRmaker* git repository
2) Create a development environment
3) Install *SFRmaker* dependencies
4) Installing the SFRmaker source code
5) Make changes to code and add tests
6) Update the documentation
7) Submit a Pull Request

Each of these 7 steps is detailed below.


1) Forking the *SFRmaker* repository using Git
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To the new user, working with Git is one of the more daunting aspects of contributing to *SFRmaker*.
It can very quickly become overwhelming, but sticking to the guidelines below will help keep the process
straightforward and mostly trouble free.  As always, if you are having difficulties please
feel free to ask for help.

The code is hosted on `GitHub <https://github.com/usgs/sfrmaker>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* Software Carpentry's `Git Tutorial <http://swcarpentry.github.io/git-novice/>`_
* `Atlassian <https://www.atlassian.com/git/tutorials/what-is-version-control>`_
* the `GitHub help pages <http://help.github.com/>`_.
* Matthew Brett's `Pydagogue <http://matthew-brett.github.com/pydagogue/>`_.

Getting started with Git
*************************

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
*************************

You will need your own fork to work on the code. Go to the `SFRmaker project
page <https://github.com/sfrmaker/sfrmaker>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone git@github.com:your-user-name/sfrmaker.git sfrmaker-yourname
    cd sfrmaker-yourname
    git remote add upstream https://github.com/sfrmaker/sfrmaker.git

This creates the directory `sfrmaker-yourname` and connects your repository to
the upstream (main project) *SFRmaker* repository.

The testing suite should run automatically on GitHub Actions each time code is pushed to your fork,
and will also run on submittal of your pull request.

Creating a branch
*************************

You want your master branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *SFRmaker*. You can have many shiny-new-features
and switch in between them using the git checkout command.

To update this branch, you need to retrieve the changes from the develop branch::

    git fetch upstream
    git rebase upstream/develop

This will replay your commits on top of the latest SFRmaker git develop.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  **It's a good idea to move slowly while doing this and pay attention to the messages from git.** The wrong command at the wrong time can quickly spiral into a confusing mess.

If you have uncommitted changes, you will need to ``stash`` them prior
to updating.  This will effectively store your changes and they can be reapplied
after updating.

.. _contributing.dev_env:

2 & 3) Creating a development environment with the required dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A development environment is a virtual space where you can keep an independent installation of *SFRmaker*.
This makes it easy to keep both a stable version of python in one place you use for work, and a development
version (which you may break while playing with code) in another.

An easy way to create a *SFRmaker* development environment is as follows:

- Install either `Anaconda <http://docs.continuum.io/anaconda/>`_ or
  `miniconda <http://conda.pydata.org/miniconda.html>`_
- Make sure that you have :ref:`cloned the repository <contributing.forking>`
- ``cd`` to the *sfrmaker** source directory

Tell conda to create a new environment, named ``sfrmaker_dev``, that has all of the python packages needed to contribute to SFRmaker. Note that in the `geopandas instructions <https://geopandas.org/contributing.html>`_, this step is broken into two parts- 2) creating the environment, and 3) installing the dependencies. By using a yaml file that includes the environment name and package requirements, these two steps can be combined::

      conda env create -f requirements-dev.yml

This will create the new environment, and not touch any of your existing environments,
nor any existing python installation.

To work in this environment, you need to ``activate`` it. The instructions below
should work for both Windows, Mac and Linux::

      conda activate sfrmaker_dev

Once your environment is activated, you will see a confirmation message to
indicate you are in the new development environment.

To view your environments::

      conda info -e

To return to your home root environment::

      conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

At this point you can easily do a *development* install, as detailed in the next sections.


4) Installing the SFRmaker source code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once dependencies are in place, install the SFRmaker source code by navigating to the gitclone of the *SFRmaker* repository and (with the ``sfrmaker_dev`` environment activated) running::

    python install -e .


5) Making changes and writing tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*SFRmaker* is serious about testing and strongly encourages contributors to embrace
`test-driven development (TDD) <http://en.wikipedia.org/wiki/Test-driven_development>`_.
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

In general, tests are required for code pushed to *SFRmaker*.  Therefore,
it is worth getting in the habit of writing tests ahead of time so this is never an issue.

*SFRmaker* uses the `pytest testing system
<http://doc.pytest.org/en/latest/>`_ and the convenient
extensions in `numpy.testing
<http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_ and `pandas.testing <https://pandas.pydata.org/pandas-docs/stable/reference/general_utility_functions.html>`_.

Writing tests
***************

All tests should go into the ``tests`` directory. This folder contains many
current examples of tests, and we suggest looking to these for inspiration. In general,
the tests in this folder aim to be organized by module (e.g. ``test_routing.py`` for the functions in ``routing.py``).

Running the test suite
******************************

The tests can then be run directly inside your Git clone (without having to
install *SFRmaker*) by typing::

    pytest

6) Updating the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The *SFRmaker* documentation resides in the `docs` folder. Changes to the docs are
made by modifying the appropriate file in the `source` folder within `docs`.
The *SFRmaker* docs use reStructuredText syntax, `which is explained here <http://www.sphinx-doc.org/en/stable/rest.html#rst-primer>`_
and the docstrings follow the `Numpy Docstring standard <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

Once you have made your changes, you can try building the docs using sphinx. To do so, you can navigate to the `doc` folder and type::

    make -C docs html

The resulting html pages will be located in `docs/build/html`. It's a good practice to rebuild the docs often while writing to stay on top of any mistakes. The `reStructuredText extension in VS Code <https://marketplace.visualstudio.com/items?itemName=lextudio.restructuredtext>`_ is another way to continuously preview a rendered documentation page while writing.


7) Submitting a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you've made changes and pushed them to your forked repository, you then
submit a pull request to have them integrated into the *SFRmaker* code base.

You can find a pull request (or PR) tutorial in the `GitHub's Help Docs <https://help.github.com/articles/using-pull-requests/>`_.

.. note::
   Please make pull requests to the develop branch; the master branch is typically only updated when releases are made to PyPI.

.. note::
   If you haven't already, it's also good practice to rebase the feature branch that your are working on (e.g. ``shiny-new-feature`` above) to the develop branch on the main repository. This brings your branch into sync with the main develop branch, and places the commits you are adding (the new feature or bug fix, etc) on top of it. If you see a message in the pull request dialog on GitHub such as "This branch cannot be rebased due to conflicts", rebasing your feature branch to the develop branch on the main repository should fix it. Note that you may have to use the ``-f`` flag to force push to your fork after rebasing. For example:

   .. code-block:: 

     git fetch upstream
     git rebase upstream/develop
     git push origin shiny-new-feature

   (this assumes that your remote pointing to the main repository is called *upstream*, and your remote pointing to your fork is called *origin*) 

.. _contributing_style:

Style Guide & Linting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SFRmaker tries to follow the `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ standard. At this point, there's no enforcement of this, but I am considering implementing `Black <https://black.readthedocs.io/en/stable/>`_, which automates a code style that is PEP8-complient. Many editors perform automatic linting that makes following PEP8 easy.

SFRmaker does use the `isort <https://github.com/timothycrosley/isort>`_ package to automatically organize import statements. isort can installed via pip::

   $ pip install isort

And then run with::

   $ isort .

from the root level of the project.

Optionally (but recommended), you can setup `pre-commit hooks <https://pre-commit.com/>`_
to automatically run ``isort`` when you make a git commit. This
can be done by installing ``pre-commit``::

   $ python -m pip install pre-commit

From the root of the SFRmaker repository, you should then install the
``pre-commit`` included in *SFRmaker*::

   $ pre-commit install

Then ``isort`` will be run automatically each time you commit changes. You can skip these checks with ``git commit --no-verify``.

.. _filing an issue: https://docs.github.com/en/free-pro-team@latest/github/managing-your-work-on-github/creating-an-issue
.. _Issues tab: https://github.com/usgs/sfrmaker/issues