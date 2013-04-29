from distutils.core import setup

# Read the version number
with open("coo_utils/_version.py") as f:
    exec(f.read())

setup(
    name='coo_utils',
    version=__version__, # use the same version that's in _version.py
    author='David N. Mashburn',
    author_email='david.n.mashburn@gmail.com',
    packages=['coo_utils'],
    scripts=[],
    url='http://pypi.python.org/pypi/coo_utils/',
    license='LICENSE.txt',
    description='utilities for managing nested lists of lists of scipy.sparse matrices',
    long_description=open('README.rst').read(),
    install_requires=[
                      'numpy>=1.0',
                      'scipy>=0.8',
                     ],
)
