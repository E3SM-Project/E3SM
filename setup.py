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

meridional_mean_2d_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'meridional_mean_2d*')

lat_lon_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'lat_lon*')
lat_lon_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'lat_lon*')
lat_lon_files.append('acme_diags/plot/vcs/plot_set_5.json')
lat_lon_files.append('acme_diags/plot/vcs/plot_set_5_new.json')

lat_lon_vector_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'lat_lon_vector*')

polar_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'polar*')
polar_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'polar*')
polar_files.append('acme_diags/plot/vcs/plot_set_7.json')
polar_files.append('acme_diags/plot/vcs/plot_set_7_new.json')

cosp_histogram_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'cosp_histogram*')
cosp_histogram_files += get_all_files_in_dir('acme_diags/driver/default_diags/legacy_diags', 'cosp_histogram*')

area_mean_time_series = get_all_files_in_dir('acme_diags/driver/default_diags', 'area_mean_time_series*')
qbo = get_all_files_in_dir('acme_diags/driver/default_diags', 'qbo*')

enso_diags_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'enso_*')

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
    (os.path.join(INSTALL_PATH, 'meridional_mean_2d'),
     meridional_mean_2d_files
     ),
    (os.path.join(INSTALL_PATH, 'lat_lon'),
     lat_lon_files
     ),
    (os.path.join(INSTALL_PATH, 'polar'),
     polar_files
     ),
    (os.path.join(INSTALL_PATH, 'lat_lon_vector'),
     lat_lon_vector_files
     ),
    (os.path.join(INSTALL_PATH, 'cosp_histogram'),
     cosp_histogram_files
     ),
    (os.path.join(INSTALL_PATH, 'area_mean_time_series'),
     area_mean_time_series
     ),
    (os.path.join(INSTALL_PATH, 'enso_diags'),
     enso_diags_files
     ),
    (os.path.join(INSTALL_PATH, 'qbo'),
     qbo
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
     ),
    (os.path.join(INSTALL_PATH, 'viewer'),
     ['acme_diags/viewer/index_template.html'
      ])
]

setup(
    name="e3sm_diags",
    version="2.0.0",
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
