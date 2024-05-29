from setuptools import setup, Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True

# For all of this to work, we either need to:
# 1. build the extension as we build the python package; or
# 2. simply package the extension as "data".
# The latter is actually quite annoying, but we may have to do it anyway.
# So, bring the extension to the dir of this file, before triggering

setup(
    name='screaminpy',
    version='1.0.0',
    author="screamers",
    description='screaminpy wrapper',
    packages=['', 'screaminpy'],
    package_data={
        'screaminpy': ['*.py'],
        '': ['pyscream.*so'],
    },
    distclass=BinaryDistribution
)

# TODOs:
# 1. We need to figure out how to handle the complex env needed by scream
# 2. We need to think forward about packaging this neatly

# for both issues, I think the most straightforward solution for now is the following:
# We simply build manually-ish for support machines and just upload the wheels
# to pypi. That way, we can ensure the build works on the machine,
# we can set its environment as we wish, etc..
# The downside is that it is manual...