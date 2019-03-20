import sys
import os
import glob
from setuptools import find_packages, setup


def get_all_files_in_dir(directory, pattern):
    return glob.glob(os.path.join(directory, pattern))


zonal_mean_xy_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'zonal_mean_xy*')
zonal_mean_xy_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'zonal_mean_xy*')
zonal_mean_xy_files.append('acme_diags/plot/vcs/plot_set_3.json')

zonal_mean_2d_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'zonal_mean_2d*')
zonal_mean_2d_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'zonal_mean_2d*')
zonal_mean_2d_files.append('acme_diags/plot/vcs/plot_set_4.json')
zonal_mean_2d_files.append('acme_diags/plot/vcs/plot_set_4_new.json')

lat_lon_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'lat_lon*')
lat_lon_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'lat_lon*')
lat_lon_files.append('acme_diags/plot/vcs/plot_set_5.json')
lat_lon_files.append('acme_diags/plot/vcs/plot_set_5_new.json')

polar_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'polar*')
polar_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'polar*')
polar_files.append('acme_diags/plot/vcs/plot_set_7.json')
polar_files.append('acme_diags/plot/vcs/plot_set_7_new.json')

cosp_histogram_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'cosp_histogram*')
cosp_histogram_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'cosp_histogram*')

rgb_files = get_all_files_in_dir('acme_diags/plot/colormaps', '*.rgb')
control_runs_files = get_all_files_in_dir('acme_diags/driver/control_runs', '*.csv')

INSTALL_PATH = 'share/e3sm_diags/'

data_files = [
    (os.path.join(INSTALL_PATH, 'zonal_mean_xy'),
     zonal_mean_xy_files
     ),
    (os.path.join(INSTALL_PATH, 'zonal_mean_2d'),
     zonal_mean_2d_files
     ),
    (os.path.join(INSTALL_PATH, 'lat_lon'),
     lat_lon_files
     ),
    (os.path.join(INSTALL_PATH, 'polar'),
     polar_files
     ),
    (os.path.join(INSTALL_PATH, 'cosp_histogram'),
     cosp_histogram_files
     ),
    (INSTALL_PATH,
     ['acme_diags/driver/acme_ne30_ocean_land_mask.nc',
      'misc/e3sm_logo.png'
      ]),
    (os.path.join(INSTALL_PATH, 'colormaps'),
     rgb_files
     ),
    (os.path.join(INSTALL_PATH, 'control_runs'),
     control_runs_files
     )
]

setup(
    name="e3sm_diags",
    version="1.6.1",
    author="Chengzhu (Jill) Zhang, Zeshawn Shaheen",
    author_email="zhang40@llnl.gov, shaheen2@llnl.gov",
    description="E3SM Diagnostics",
    scripts=["acme_diags/acme_diags_driver.py"],
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    data_files=data_files,
    entry_points={
        'console_scripts': [
            'e3sm_diags=acme_diags.acme_diags_driver:main',
            'acme_diags=acme_diags.acme_diags_driver:main',
            'e3sm_diags_vars=acme_diags.acme_diags_vars:main'
    ]}
)
