.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

# To run these commands: make <COMMAND>
# ==================================================

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

# Clean local repository
# ----------------------
clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr conda-build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr tests_coverage_reports/
	rm -f .coverage
	rm -fr htmlcov/
	rm -f coverage.xml
	rm -fr .pytest_cache
	rm -rf .mypy_cache

clean-test-int-res: ## remove integration test results and image check failures
	rm -rf tests/integration/all_sets_results_test
	rm -rf tests/integration/image_check_failures

clean-test-int-data:  # remove integration test data and images (expected) -- useful when they are updated
	rm -rf tests/integration/integration_test_data
	rm -rf tests/integration/integration_test_images

# Quality Assurance
# ----------------------
pre-commit:  # run pre-commit quality assurance checks
	pre-commit run --all-files

lint: ## check style with flake8
	flake8 e3sm_diags tests

test: ## run tests quickly with the default Python and produces code coverage report
	pytest
	$(BROWSER) tests_coverage_reports/htmlcov/index.html

# Documentation
# ----------------------
docs: ## generate Sphinx HTML documentation, including API docs
	rm -rf docs/generated
	cd docs && make html
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

# Build
# ----------------------
install: clean ## install the package to the active Python's site-packages
	python -m pip install .
