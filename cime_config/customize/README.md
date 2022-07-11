# CIME customization
This directory contains the code for customizing CIME's runtime and provenance collection.

## flags.py
Contains the flags that alter CIME's runtime behavior.

## provenance.py
Contains the code for capturing build and pre/post run provenance.

Implements three provenance hooks.
- `save_build_provenance`
- `save_prerun_provenance`
- `save_postrun_provenance`

## Testing
The `tests/` directory contains unit tests for provenance code.

Requirements:
- `pytest`

### Example
```bash
pip install pytest
# or
# conda install -c conda-forge pytest

pytest -vvv --machine docker tests/
```
