from setuptools import find_packages, setup

setup(
    name="acme_diags",
    version="0.1",
    author="Zeshawn Shaheen, Jill Zhang",
    author_email="aims@llnl.gov",
    description="ACME Diagnostics.",
    scripts=["acme_diags/plotset5/set5_driver.py"],
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    package_data={"": ["*.json"]},
)
