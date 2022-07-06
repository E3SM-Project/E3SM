# CIME customization
This module contains the code for customizing CIME's runtime.

## flags.py
This file contains the flags that alter CIME's behavior.

## provenance.py
Contains the code for capturing build and pre/post run provenance.

## Testing
Tests for the provenance code are contained in the `tests/` directory.

To run the tests install the `pytest` package and run the following command.

```bash
pytest -vvv --machine <machine_name> tests
```
