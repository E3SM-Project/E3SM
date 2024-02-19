# Developing Docs for E3SM

The E3SM source code includes the documentation of the model in text files stored along with source code in the main repository. Please follow project guidelines to develop on feature branches.

Our documentation is written in Markdown language. See [Basic writing and formatting syntax - GitHub Docs](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to get started.

## Local development environment

Just like you compile and test your Fortran before committing it, you should "compile" and view your documentation before committing and pushing the branch. Unlike compiling and running the model, the documentation can be easily built and displayed on your laptop. Here are two way to do that.

### Virtual environment (venv)

You can create a local python virtual environment for the doc building packages (which is a one-time setup). On your local laptop in your home directory, create a python virtual environment called “mkdocsenv” for example `python3 -m venv mkdocsenv`. Now you can enter the environment with `source mkdocsenv/bin/activate`. To build and visualize the docs, please install the python packages needed:

```shell
pip install mkdocs-material mkdocs-monorepo-plugin pymdown-extensions mdutils mkdocs-bibtex
```

At any time, you can leave the environment by typing: `deactivate`.

The creation of mkdocsenv and the pip install commands only have to be done once on each platform you want to do documentation development.

Installing a python package in a virtual environment keeps it from being installed system wide which is good practice and also useful if you don’t have permissions to install system-level packages.

### Conda environment

The same python environment above can be achieve with conda. To do so, create a conda environment with a `conda create -n e3sm_docs mkdocs-material pymdown-extensions mkdocs-monorepo-plugin mdutils mkdocs-bibtex`, and then `conda activate e3sm_docs` to use the environment.

## Local building and testing

Activate your virtual environment that you made in the above step with either `source mkdocsenv/bin/activate` or `conda activate e3sm_docs`. Change directory to the top level E3SM directory and see if you can build the current docs by typing

```shell
mkdocs build
```

Note that EAMxx docs use an automated routine to populate the namelist parameters. In order to avoid warnings on this step, please ensure you have a copy of the CIME submodule inside your main repo. Then, you will be able to build with `mkdocs build --strict --verbose` (to convert any warning into an error, which we use in the PR checking).

To View the docs locally in your browser with the built-in server, you can type

```shell
mkdocs serve
```

Follow the instructions output by the above command to display in your browser. To exist the server, hit Control-C.

An alternative to mkdocs serve is to access to the docs directly without serving them through a URL. To achieve this, one could you `file:///E3SM_ROOT/site/index.html` in a modern browser like Chrome (in place of a URL). In this case, E3SM_ROOT is the full path to the clone E3SM repository on your local machine. One could directly go to other pages as needed using the same method (e.g., `file:///E3SM_ROOT/site/EAM/index.html`).

As you edit the documentation, rerun the `mkdocs build` command and refresh the website in your browser.  The output of `mkdocs build` will tell you if you have a syntax error in  your markdown files.  The website in your browser will let you know if it looks as you intended.

## Documentation organization

The documentation for a model should be stored within that model's source code directory.

The first file needed is mkdocs.yml which should be in `E3SM/components/<model>`.  See `E3SM/components/eam/mkdocs.yml` for an example:

```yaml
site_name: EAM

nav:
 - Introduction: 'index.md'
 - Users Guide: 'user-guide/index.md'
 - Developer Guide: 'dev-guide/index.md'
 - Technical Guide: 'tech-guide/index.md'
```

The "site_name" is the keyword used in the `E3SM/docs/index.md` file to reference all the model's documentation. (That file must be edited when adding a new model’s documentation.)

All other files should be stored under `E3SM/components/<model>/docs`

The basic structure under "docs" is

- `index.md`: top level introductory text and links to main sections
- `dev-guide/`: subdir for development guide
- `tech-guide/`: subidr for technical guide
- `user-guide/`: subdir for user guide

Each of the sub-dirs should also have an `index.md` file to organize their content.

## Files and sections

A subsection of your documentation should be contained entirely in a single file. Short sections can be grouped in to one file.

Equation numbering and references to equations are local to an individual markdown file.

## Figures

When you start adding figures to your documentation, keep in mind that binary files can lead to repo size bloat IF they change a lot so don’t commit a figure to the repo until you are confident its the near-final version.

Figures should be in png format and under 500 pixels on the longest side.

Keep figures in the same directory as the Markdown file referencing them.

## Markdown style guide

To help with maintainability, everyone is encouraged to follow best practices when it comes to writing markdown files. For example, popular IDEs support [`markdownlint`](https://github.com/DavidAnson/markdownlint).

## Automated preview during PRs

When you submit a PR, we have an automated process to build and test the docs (with `--strict --verbose`) and to preview them directly in the PR. An automated bot will leave a comment with a URL for you to review the docs.
