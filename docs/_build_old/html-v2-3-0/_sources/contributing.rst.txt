**************************
Contributing to E3SM Diags
**************************

.. highlight:: none

Code
====

We would welcome you feedback and suggestions on how to improve E3SM Diags.
Drop a line to Jill (zhang40 .at. llnl.gov) or Ryan (forsyth2 .at. llnl.gov).
 
Or would you like to directly modify E3SM Diags? We tried to build the code in such a way
to make it easy to modify and add new backends, diagnostics, plots, variables, etc.
We would welcome your contributions. Create a :ref:`development environment <dev-env>` 
and start coding away. Then submit your suggested modifications via a pull request.

Documentation
=============

Something incorrect? Something missing? Want to add some cool examples? 
We would welcome your help updating and maintaining this documentation.
Here is a quick guide on how to get started doing just that.

Create a conda environment
--------------------------

This documentation is created using 
`Sphinx <http://www.sphinx-doc.org/en/stable>`_. Sphinx is an open-source tool 
that makes it easy to create intelligent and beautiful documentation, written 
by Georg Brandl and licensed under the BSD license.

Create a new conda environment, install Sphinx as well as the 
Sphinx `readthedocs theme <https://github.com/rtfd/sphinx_rtd_theme>`_. ::

   $ conda create -n sphinx
   $ conda activate sphinx
   $ conda install -c anaconda sphinx
   $ pip install sphinx_rtd_theme

Optionally, if you need to convert Jupyter notebooks to import in Sphinx,
you'll also need to install `nbconvert <https://nbconvert.readthedocs.io/en/latest/#>`_ 
and `pandoc <https://pandoc.org/>`_: ::

   $ conda install -c conda-forge nbconvert pandoc

Checkout and edit documentation
-------------------------------

The documentation is maintained in the Github repository under a separate
branch named 'gh-pages'. This special branch is used by Github.com to directly
serve static web pages (see `GitHub Pages <https://pages.github.com/>`_).

Clone the repository and checkout the 'gh-pages' branch: ::

   $ cd <myDir>
   $ git clone git@github.com:E3SM-Project/e3sm_diags.git e3sm_diags
   $ cd e3sm_diags
   $ git fetch origin gh-pages
   $ git checkout -b <branch-name> origin/gh-pages

You should now see two sub-directories: `source` contains the documentation
source files, and `docs` the html web pages created by Sphinx.

To modify the documentation, simply edit the files under `source`.
Sphinx uses `reStructuredText <http://docutils.sourceforge.net/rst.html>`_ 
as its markup language. For more information on how to write documentation 
using Sphinx, you can refer to

* `First Steps with Sphinx <http://www.sphinx-doc.org/en/stable/tutorial.html>`_
* `reStructuredText Primer <http://www.sphinx-doc.org/en/stable/rest.html#external-links>`_

To see the changes you made to the documentation, rebuild the web pages ::

   $ cd <myDir>/e3sm_diags
   $ python source/quickguides/generate_quick_guides.py # Run if quick-guide files were updated
   $ make html # Make sure you have run `conda activate sphinx` first
 
and view them locally in a web browser at `file:///<myDir>/e3sm_diags/index.html`.

Sphinx may occasionally not build the new files properly. If that is the case,
try removing the `docs` sub-directory (be careful not to remove `source`)
and rebuild entirely: ::

   $ cd <myDir>/e3sm_diags
   $ rm -r docs
   $ python source/quickguides/generate_quick_guides.py # Run if quick-guide files were updated
   $ make html # Make sure you have run `conda activate sphinx` first
 

Once you are satisfied with your modifications, commit and push them back to 
the repository: ::

   $ cd <myDir>/e3sm_diags
   $ git add .
   $ git commit
   $ git push <your-fork-remote-name> <branch-name> # If not using a fork, use `origin`

Then, create a pull request from ``your-fork/e3sm_diags/branch-name`` to ``E3SM-Project/e3sm_diags/gh-pages``.
   
Once this pull request is merged, changes will immediately be available on the
`e3sm_diags documentation page <https://e3sm-project.github.io/e3sm_diags/>`_.

Converting Jupyter notebooks
----------------------------

If you have Jupyter notebooks that you'd like to import into the documentation,
they can easily be converted to rst format: ::

   $ jupyter nbconvert mygreatnotebook.ipynb --to rst

Initial setup (for reference only)
----------------------------------

The instructions below only apply for the initial configuration of the
Sphinx documentation on the Github repository. They are documented here
for reference only. Do not follow them unless you are setting up documentation
for a new repository. (Adapted from `Sphinx documentation on GitHub 
<http://datadesk.latimes.com/posts/2012/01/sphinx-on-github>`_.)

Create Sphinx conda environment (see above).

Create a new git branch (gh-pages): ::

  $ git branch gh-pages
  $ git checkout gh-pages

Clear out any­thing from the master branch and start fresh ::

  $ git symbolic-ref HEAD refs/heads/gh-pages
  $ rm .git/index
  $ git clean -fdx

Create documentation ::

  $ sphinx-quickstart

accept suggested default options, except ::

  Separate source and build directories (y/N) [n]: y

Edit Makefile and change BUILDIR ::

  BUILDDIR = docs

Remove old build directory ::

  $ rmdir build

Change the Sphinx theme to 'ReadTheDocs'. Edit 'source/conf.py and change ::

  html_theme = 'alabaster'

to ::

  import sphinx_rtd_theme
  html_theme = "sphinx_rtd_theme"
  html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

Try building documentation ::

  $ make html

Create an empty .nojekyll file to indicate to Github.com that this
is not a Jekyll static website: ::

  $ touch .nojekyll

Create a top-level re-direction file: ::

  $ vi index.html

with the following: ::

  <meta http-equiv="refresh" content="0; url=./docs/html/index.html" />

Commit and push back to Github: ::

  $ git add .
  $ git commit
  $ git push origin gh-pages

