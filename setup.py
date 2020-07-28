# -*- coding: utf-8 -*-
"""
AnisoCADO: A python package to simulate ELT SCAO PSFs

How to compile and put these on pip::

    $ python setup.py sdist bdist_wheel
    $ twine upload dist/*

"""
from setuptools import setup, find_packages

with open('README.md') as f:
    __readme__ = f.read()

with open('LICENSE') as f:
    __license__ = f.read()

with open('anisocado/version.py') as f:
    __version__ = f.readline().split("'")[1]


print(__version__)


def setup_package():
    setup(name='anisocado',
          version=__version__,
          description='Generate off-axis SCAO PSFs for the ELT',
          long_description=__readme__,
          long_description_content_type="text/markdown",
          author='Eric Gendron, Kieran Leschinski',
          author_email='kieran.leschinski@univie.ac.at',
          url='https://github.com/astronomyk/anisocado',
          license="GNU general public license",
          include_package_data=True,
          packages=find_packages(exclude=('tests', 'docs')),
          package_dir={'anisocado': 'anisocado'},
          install_requires=['numpy', 'astropy', 'matplotlib'],
          classifiers=["Programming Language :: Python :: 3",
                       "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                       "Operating System :: OS Independent",
                       "Intended Audience :: Science/Research",
                       "Topic :: Scientific/Engineering :: Astronomy", ]
          )


if __name__ == '__main__':
    setup_package()
