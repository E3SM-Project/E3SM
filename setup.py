import sys
from setuptools import find_packages, setup

data_files = [
    (sys.prefix + '/share/acme_diags/set5',
     ['acme_diags/plotset5/set5_diags_AMWG_default.json',
      'acme_diags/plotting/set5/plot_set_5_new.json',
      'acme_diags/plotting/set5/plot_set_5.json'
     ]
    ),
]

setup(
    name="acme_diags",
    version="0.1",
    author="Zeshawn Shaheen, Jill Zhang",
    author_email="aims@llnl.gov",
    description="ACME Diagnostics.",
    scripts=["acme_diags/plotset5/set5_driver.py"],
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    data_files=data_files
)
