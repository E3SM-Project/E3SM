import sys
import os
import glob
from setuptools import find_packages, setup


def get_all_files_in_dir(directory, pattern):
    return glob.glob(os.path.join(directory, pattern))


rgb_files = get_all_files_in_dir('acme_diags/plot/colormaps', '*.rgb')
control_runs_files = get_all_files_in_dir('acme_diags/driver/control_runs', '*.csv')

data_files = [
    (sys.prefix + '/share/acme_diags/zonal_mean_xy',
     ['acme_diags/driver/default_diags/zonal_mean_xy_AMWG.json',
      'acme_diags/driver/default_diags/zonal_mean_xy_ACME.json',
      'acme_diags/plot/vcs/plot_set_3.json'
      ]),
    (sys.prefix + '/share/acme_diags/zonal_mean_2d',
     ['acme_diags/driver/default_diags/zonal_mean_2d_AMWG.json',
      'acme_diags/driver/default_diags/zonal_mean_2d_ACME.json',
      'acme_diags/plot/vcs/plot_set_4.json',
      'acme_diags/plot/vcs/plot_set_4_new.json'
      ]),
    (sys.prefix + '/share/acme_diags/lat_lon',
     ['acme_diags/driver/default_diags/lat_lon_AMWG.json',
      'acme_diags/driver/default_diags/lat_lon_ACME.json',
      'acme_diags/plot/vcs/plot_set_5_new.json',
      'acme_diags/plot/vcs/plot_set_5.json'
      ]),
    (sys.prefix + '/share/acme_diags/polar',
     ['acme_diags/driver/default_diags/polar_AMWG.json',
      'acme_diags/driver/default_diags/polar_ACME.json',
      'acme_diags/plot/vcs/plot_set_7_new.json',
      'acme_diags/plot/vcs/plot_set_7.json'
      ]),
    (sys.prefix + '/share/acme_diags/cosp_histogram',
     ['acme_diags/driver/default_diags/cosp_histogram_AMWG.json',
      'acme_diags/driver/default_diags/cosp_histogram_ACME.json'
      ]),
    (sys.prefix + '/share/acme_diags',
     ['acme_diags/driver/acme_ne30_ocean_land_mask.nc',
      'misc/ACME_Logo.png'
      ]),
    (sys.prefix + '/share/acme_diags/colormaps',
     rgb_files
     ),
    (sys.prefix + '/share/acme_diags/control_runs', 
     control_runs_files
     )
]

setup(
    name="acme_diags",
    version="1.2.1",
    author="Chengzhu (Jill) Zhang, Zeshawn Shaheen",
    author_email="zhang40@llnl.gov, shaheen2@llnl.gov",
    description="ACME Diagnostics.",
    scripts=["acme_diags/acme_diags_driver.py"],
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    data_files=data_files
)
