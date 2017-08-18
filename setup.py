import sys
from setuptools import find_packages, setup

data_files = [

    (sys.prefix + '/share/acme_diags/zonal_mean_xy',
     ['acme_diags/driver/zonal_mean_xy_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_3.json'
     ]),
    (sys.prefix + '/share/acme_diags/zonal_mean_2d',
     ['acme_diags/driver/zonal_mean_2d_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_4.json',
      'acme_diags/plot/vcs/plot_set_4_new.json'
     ]),
    (sys.prefix + '/share/acme_diags/lat_lon',
     ['acme_diags/driver/lat_lon_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_5_new.json',
      'acme_diags/plot/vcs/plot_set_5.json'
     ]),
    (sys.prefix + '/share/acme_diags/polar',
     ['acme_diags/driver/polar_AMWG_default.json',
      'acme_diags/plot/vcs/plot_set_7_new.json',
      'acme_diags/plot/vcs/plot_set_7.json'
     ]),
    (sys.prefix + '/share/acme_diags/cosp_histogram',
     ['acme_diags/driver/cosp_histogram_AMWG_default.json'
     ]),
    (sys.prefix + '/share/acme_diags',
     ['acme_diags/driver/acme_ne30_ocean_land_mask.nc',
      'docs/ACME_Logo.png'
     ])
]

setup(
    name="acme_diags",
    version="0.1.0",
    author="Chengzhu (Jill) Zhang, Zeshawn Shaheen",
    author_email="aims@llnl.gov",
    description="ACME Diagnostics.",
    scripts=["acme_diags/acme_diags_driver.py"],
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    data_files=data_files
)
