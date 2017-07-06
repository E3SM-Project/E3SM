import sys
from setuptools import find_packages, setup

data_files = [

    (sys.prefix + '/share/acme_diags/set4',
     ['acme_diags/driver/set4_diags_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_4_new.json',
      'acme_diags/plot/vcs/plot_set_4.json'
     ]),
    (sys.prefix + '/share/acme_diags/set3',
     ['acme_diags/driver/set3_diags_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_3.json'
     ]),
    (sys.prefix + '/share/acme_diags/set5',
     ['acme_diags/driver/set5_diags_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_5_new.json',
      'acme_diags/plot/vcs/plot_set_5.json'
     ]),
    (sys.prefix + '/share/acme_diags/set7',
     ['acme_diags/driver/set7_diags_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_7_new.json',
      'acme_diags/plot/vcs/plot_set_7.json'
     ]),
    (sys.prefix + '/share/acme_diags',
     ['acme_diags/driver/acme_ne30_ocean_land_mask.nc'
     ])
]

setup(
    name="acme_diags",
    version="0.1",
    author="Zeshawn Shaheen, Chengzhu (Jill) Zhang",
    author_email="aims@llnl.gov",
    description="ACME Diagnostics.",
    scripts=["acme_diags/acme_diags_driver.py"],
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    data_files=data_files
)
