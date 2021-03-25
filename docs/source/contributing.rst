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

Getting Started
--------------------------

This documentation is created using
`Sphinx <http://www.sphinx-doc.org/en/stable>`_. Sphinx is an open-source tool
that makes it easy to create intelligent and beautiful documentation, written
by Georg Brandl and licensed under the BSD license.

The documentation is maintained in the ``master`` branch of the GitHub repository.
You can include code and its corresponding documentation updates in a single pull request (PR).

After merging a PR, GitHub Actions automates the documentation building process.
It pushes the HTML build to the ``gh-pages`` branch, which is hosted on GitHub Pages.

Edit Documentation
-------------------------------

Sphinx uses `reStructuredText <http://docutils.sourceforge.net/rst.html>`_
as its markup language. For more information on how to write documentation
using Sphinx, you can refer to

* `First Steps with Sphinx <http://www.sphinx-doc.org/en/stable/tutorial.html>`_
* `reStructuredText Primer <http://www.sphinx-doc.org/en/stable/rest.html#external-links>`_

To begin editing:

1. Set up your local :ref:`development environment <dev-env>` if you haven't already.

2. To modify the documentation, simply edit the files under ``docs/source``.

3. To see the changes you made to the documentation, rebuild the web pages ::

    cd <myDir>/e3sm_diags
    conda activate e3sm_diags_env_dev
    cd docs
    python source/quickguides/generate_quick_guides.py # Run if quick-guide files were updated
    make html # no version selector dropdown will be generated

4. View them locally in a web browser at ``file:///<myDir>/e3sm_diags/docs/_build/html/index.html``.

5. Once you are satisfied with your modifications, commit and push changes to the repository: ::

    cd <myDir>/e3sm_diags
    # `docs/_build` is ignored by git since it does not need to be pushed
    git add .
    git commit "..."
    git push <fork-origin> <branch-name>

6. <`OPTIONAL`> If you want to generate and view versioned docs: ::

    # After commiting to your branch
    cd docs
    sphinx-multiversion source _build/html
    # Check the `_build/html` folder for all generated versioned docs
    # Open `_build/html/<your-branch>/index.html` to view in browser

   .. figure:: _static/docs-version-selector.png
      :alt: Docs version selector

      Docs version selector dropdown in the bottom left-hand corner

7. Create a pull request from ``<your-fork>/e3sm_diags/branch-name`` to ``E3SM-Project/e3sm_diags/master``.

Once this pull request is merged and GitHub Actions finishes building the docs, changes will be available on the
`e3sm_diags documentation page <https://e3sm-project.github.io/e3sm_diags/>`_.

How Documentation is Versioned
----------------------------
The `sphinx-multiversion <https://github.com/Holzhaus/sphinx-multiversion>`_ package manages documentation versioning.

``sphinx-multiversion`` is configured to generate versioned docs for available tags and
branches on local, ``origin`` (the likely name for your ``<fork-origin>``) and
``upstream`` (the likely name for your ``<upstream-origin>``).

Branches or tags that don’t contain both the sphinx ``source`` directory and the ``conf.py`` file will be skipped automatically.

    - Skipped versions include releases ``< v2.5.0`` since the documention source was not included in those tagged releases.
    - Run ``sphinx-multiversion source _build/html --dump-metadata`` to see which tags/branches matched.

Converting Jupyter Notebooks
----------------------------

If you have Jupyter notebooks that you'd like to import into the documentation,
they can easily be converted to rst format: ::

   $ jupyter nbconvert mygreatnotebook.ipynb --to rst

Initial setup (obsolete/for reference only)
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

Clear out anything from the master branch and start fresh ::

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

