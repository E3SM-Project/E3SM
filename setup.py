import sys
import os
import glob
from setuptools import find_packages, setup


def get_all_files_in_dir(directory, pattern):
    return glob.glob(os.path.join(directory, pattern))


zonal_mean_xy_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'zonal_mean_xy*')
zonal_mean_xy_files.append('acme_diags/plot/vcs/plot_set_3.json')

zonal_mean_2d_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'zonal_mean_2d*')
zonal_mean_2d_files.append('acme_diags/plot/vcs/plot_set_4.json')
zonal_mean_2d_files.append('acme_diags/plot/vcs/plot_set_4_new.json')

lat_lon_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'lat_lon*')
lat_lon_files.append('acme_diags/plot/vcs/plot_set_5.json')
lat_lon_files.append('acme_diags/plot/vcs/plot_set_5_new.json')

polar_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'polar*')
polar_files.append('acme_diags/plot/vcs/plot_set_7.json')
polar_files.append('acme_diags/plot/vcs/plot_set_7_new.json')

cosp_histogram_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'cosp_histogram*')

rgb_files = get_all_files_in_dir('acme_diags/plot/colormaps', '*.rgb')
control_runs_files = get_all_files_in_dir('acme_diags/driver/control_runs', '*.csv')


data_files = [
    (sys.prefix + '/share/acme_diags/zonal_mean_xy',
     zonal_mean_xy_files
     ),
    (sys.prefix + '/share/acme_diags/zonal_mean_2d',
     zonal_mean_2d_files
     ),
    (sys.prefix + '/share/acme_diags/lat_lon',
     lat_lon_files
     ),
    (sys.prefix + '/share/acme_diags/polar',
     polar_files
     ),
    (sys.prefix + '/share/acme_diags/cosp_histogram',
     cosp_histogram_files
     ),
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
