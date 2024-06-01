from setuptools import setup, Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True

setup(
    name='pyeamxx',
    version='0.0.1.dev5',
    author='E3SM SCREAM',
    description='EAMxx wrapper',
    packages=['', 'pyeamxx'],
    package_data={
        '': ['libpyeamxx/*.so'],
    },
    distclass=BinaryDistribution,
    zip_safe=False
)
