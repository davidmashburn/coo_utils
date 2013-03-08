from distutils.core import setup

setup(
    name='coo_utils',
    version='0.1.2',
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
