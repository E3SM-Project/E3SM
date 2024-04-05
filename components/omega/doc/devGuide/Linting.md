(omega-dev-linting)=

# Linting Code

## Enabling Lint Checks

First, please follow the procedure in {ref}`omega-dev-conda-env` to set up a
conda environment for linting the code and building the docs.

Once on each machine you are developing on, before committing code, please run:
```bash
pre-commit install
```

After this, your code will be linted with all the tools discussed below as part
of each call to `git commit`. Any modified files will also be linted when CI
(GitHub Actions) runs on each pull request and merge to the `develop` branch.

You can lint the full code using all the tools with:
```bash
pre-commit run --all-files
``````

## Bypassing Linting

If you wish to commit without first checking the code for lint (e.g. you plan
to fix the code in a later commit), run:
```bash
git commit --no-verify ...
```
with the same arguments you would normally use for `git commit`.


## Linting C++ Code

The tools used to lint C++ code are from
[cmake-pre-commit-hooks](https://github.com/Takishima/cmake-pre-commit-hooks)
and include [clang-format](https://clang.llvm.org/docs/ClangFormat.html),
[clang-tidy](https://clang.llvm.org/extra/clang-tidy/),
[cppcheck](https://cppcheck.sourceforge.io/), and
[include-what-you-use](https://github.com/include-what-you-use/include-what-you-use).
(Currently `clang-tidy`, `cppcheck` and `include-what-you-use` are disabled
but they will be enabled shortly.)

You can run these tools individually if you need to:
```bash
pre-commit run clang-format --all-files
pre-commit run clang-tidy --all-files
pre-commit run cppcheck --all-files
pre-commit run include-what-you-use --all-files
```
You can specify one more more files instead of `--all-files`.

## Linting Python Code

The tools used to lint python code include
[isort](https://pycqa.github.io/isort/) for sorting imports,
[flynt](https://github.com/ikamensh/flynt) for enforcing string formatting
with f-strings,
[flake8](https://flake8.pycqa.org/en/latest/) for enforcing the
[PEP8](https://peps.python.org/pep-0008/) style guide for python, and
[mypy](https://mypy-lang.org/) for performing variable type checking.

You can run these tools individually if you need to:
```bash
pre-commit run isort --all-files
pre-commit run flynt --all-files
pre-commit run flake8 --all-files
pre-commit run mypy --all-files
```
You can specify one more more files instead of `--all-files`.

## Updating the linting pacakge

To update the linting packages, you just need to recreate the development
conda environment.  See {ref}`omega-dev-update-conda-env`.
