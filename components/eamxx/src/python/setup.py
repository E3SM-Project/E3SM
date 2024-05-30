from setuptools import setup, Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True

setup(
    name='pyeamxx',
    version='0.0.1',
    author='E3SM SCREAM',
    description='EAMxx wrapper',
    packages=['','pyeamxx'],
    package_data={
        '': ['*.*so'],
    },
    distclass=BinaryDistribution,
    zip_safe=False
)

# TODOs:
# 1. We need to figure out how to handle the complex env needed by scream
# 2. We need to think forward about packaging this neatly

# for both issues, I think the most straightforward solution for now is the following:
# We simply build manually-ish for support machines and just upload the wheels
# to pypi. That way, we can ensure the build works on the machine,
# we can set its environment as we wish, etc..
# The downside is that it is manual...