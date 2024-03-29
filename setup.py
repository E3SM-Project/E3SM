import glob
import os

from setuptools import find_packages, setup


def get_all_files_in_dir(directory, pattern):
    return glob.glob(os.path.join(directory, pattern))


zonal_mean_xy_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "zonal_mean_xy*"
)
zonal_mean_xy_files += get_all_files_in_dir(
    "e3sm_diags/driver/default_diags/legacy_diags", "zonal_mean_xy*"
)

zonal_mean_2d_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "zonal_mean_2d_model*"
)
zonal_mean_2d_files += get_all_files_in_dir(
    "e3sm_diags/driver/default_diags/legacy_diags", "zonal_mean_2d*"
)
zonal_mean_2d_stratosphere_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "zonal_mean_2d_stratosphere*"
)

meridional_mean_2d_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "meridional_mean_2d*"
)

lat_lon_files = get_all_files_in_dir("e3sm_diags/driver/default_diags", "lat_lon*")
lat_lon_files += get_all_files_in_dir(
    "e3sm_diags/driver/default_diags/legacy_diags", "lat_lon*"
)

lat_lon_vector_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "lat_lon_vector*"
)

lat_lon_land_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "lat_lon_land*"
)
lat_lon_river_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "lat_lon_river*"
)

polar_files = get_all_files_in_dir("e3sm_diags/driver/default_diags", "polar*")
polar_files += get_all_files_in_dir(
    "e3sm_diags/driver/default_diags/legacy_diags", "polar*"
)

cosp_histogram_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "cosp_histogram*"
)
cosp_histogram_files += get_all_files_in_dir(
    "e3sm_diags/driver/default_diags/legacy_diags", "cosp_histogram*"
)

area_mean_time_series = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "area_mean_time_series*"
)
qbo = get_all_files_in_dir("e3sm_diags/driver/default_diags", "qbo*")
streamflow = get_all_files_in_dir("e3sm_diags/driver/default_diags", "streamflow*")
enso_diags_files = get_all_files_in_dir("e3sm_diags/driver/default_diags", "enso_*")
diurnal_cycle_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "diurnal_cycle_*"
)
arm_diags_files = get_all_files_in_dir("e3sm_diags/driver/default_diags", "arm_diags_*")
tc_analysis_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "tc_analysis_*"
)
annual_cycle_zonal_mean_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "annual_cycle_zonal_mean_*"
)
aerosol_aeronet_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "aerosol_aeronet*cfg"
)
aerosol_aeronet_files += get_all_files_in_dir(
    "e3sm_diags/driver/default_diags/aerosol_aeronet_data", "*.txt"
)
aerosol_budget_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "aerosol_budget*cfg"
)
mp_partition_files = get_all_files_in_dir(
    "e3sm_diags/driver/default_diags", "mp_partition*cfg"
)
rgb_files = get_all_files_in_dir("e3sm_diags/plot/colormaps", "*.rgb")
control_runs_files = get_all_files_in_dir("e3sm_diags/driver/control_runs", "*")

INSTALL_PATH = "share/e3sm_diags/"

data_files = [
    (os.path.join(INSTALL_PATH, "zonal_mean_xy"), zonal_mean_xy_files),
    (os.path.join(INSTALL_PATH, "zonal_mean_2d"), zonal_mean_2d_files),
    (
        os.path.join(INSTALL_PATH, "zonal_mean_2d_stratosphere"),
        zonal_mean_2d_stratosphere_files,
    ),
    (
        os.path.join(INSTALL_PATH, "meridional_mean_2d"),
        meridional_mean_2d_files,
    ),
    (os.path.join(INSTALL_PATH, "lat_lon"), lat_lon_files),
    (os.path.join(INSTALL_PATH, "polar"), polar_files),
    (os.path.join(INSTALL_PATH, "lat_lon_vector"), lat_lon_vector_files),
    (os.path.join(INSTALL_PATH, "lat_lon_land"), lat_lon_land_files),
    (os.path.join(INSTALL_PATH, "lat_lon_river"), lat_lon_river_files),
    (os.path.join(INSTALL_PATH, "cosp_histogram"), cosp_histogram_files),
    (
        os.path.join(INSTALL_PATH, "area_mean_time_series"),
        area_mean_time_series,
    ),
    (os.path.join(INSTALL_PATH, "enso_diags"), enso_diags_files),
    (os.path.join(INSTALL_PATH, "qbo"), qbo),
    (os.path.join(INSTALL_PATH, "streamflow"), streamflow),
    (os.path.join(INSTALL_PATH, "diurnal_cycle"), diurnal_cycle_files),
    (os.path.join(INSTALL_PATH, "arm_diags"), arm_diags_files),
    (os.path.join(INSTALL_PATH, "tc_analysis"), tc_analysis_files),
    (
        os.path.join(INSTALL_PATH, "annual_cycle_zonal_mean"),
        annual_cycle_zonal_mean_files,
    ),
    (
        INSTALL_PATH,
        [
            "e3sm_diags/driver/acme_ne30_ocean_land_mask.nc",
            "misc/e3sm_logo.png",
        ],
    ),
    (os.path.join(INSTALL_PATH, "aerosol_aeronet"), aerosol_aeronet_files),
    (os.path.join(INSTALL_PATH, "aerosol_budget"), aerosol_budget_files),
    (os.path.join(INSTALL_PATH, "mp_partition"), mp_partition_files),
    (os.path.join(INSTALL_PATH, "colormaps"), rgb_files),
    (os.path.join(INSTALL_PATH, "control_runs"), control_runs_files),
    (
        os.path.join(INSTALL_PATH, "viewer"),
        ["e3sm_diags/viewer/index_template.html"],
    ),
]

setup(
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD 3-Clause License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    name="e3sm_diags",
    version="2.11.0",
    author="Chengzhu (Jill) Zhang, Tom Vo, Ryan Forsyth, Chris Golaz and Zeshawn Shaheen",
    author_email="zhang40@llnl.gov",
    description="E3SM Diagnostics",
    license="BSD 3-Clause",
    scripts=["e3sm_diags/e3sm_diags_driver.py"],
    packages=find_packages(include=["e3sm_diags", "e3sm_diags.*"]),
    data_files=data_files,
    test_suite="tests",
    entry_points={
        "console_scripts": [
            "e3sm_diags=e3sm_diags.e3sm_diags_driver:main",
            "e3sm_diags_vars=e3sm_diags.e3sm_diags_vars:main",
        ]
    },
)
