# this allows landgen to be run as a standalone script, e.g. `python -m landgen <path/config.json>`

from landgen import landgen
import sys
from pathlib import Path

if len(sys.argv) != 2:
    print(f"Usage: python {Path(__file__).name} <path/config.json>")
    sys.exit(1)
print('Executing as standalone script')
landgen.main(sys.argv[1])