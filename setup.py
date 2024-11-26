# -*- coding: utf-8 -*-

# Copyright (C) 2021 Michael Hogg

# This file is part of pygeodesic - See LICENSE.txt for information on usage and redistribution

from setuptools import setup, find_packages
from setuptools.extension import Extension
from codecs import open
from os import path
import sys

# Get current path
here = path.abspath(path.dirname(__file__))


# Function to open the readme file
def readme():
    with open(path.join(here, "README.md")) as f:
        return f.read()


long_description = readme()

# Find the version
exec(open(path.join("pygeodesic", "version.py")).read())

try:
    from Cython.Distutils import build_ext
except ImportError:
    from setuptools.command.build_ext import build_ext

    use_cython = False
else:
    use_cython = True


# Custom class to add required dependencies for building the project
class build_ext_dependencies(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)
        import numpy

        self.include_dirs.append(numpy.get_include())
        self.include_dirs.append("pygeodesic\geodesic_kirsanov")


cmdclass = {"build_ext": build_ext_dependencies}
ext_modules = []
if use_cython:
    ext_modules += [
        Extension(
            "pygeodesic.geodesic",
            sources=[
                "pygeodesic/geodesic.pyx",
            ],
            language="c++",
        )
    ]
else:
    ext_modules += [
        Extension(
            "pygeodesic.geodesic",
            sources=[
                "pygeodesic/geodesic.cpp",
            ],
            language="c++",
        )
    ]

setup(
    name="pygeodesic",
    version=__version__,
    description="Python library for calculating geodesic distance on triangular based surface meshes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT license",
    keywords=["geodesic", "distance", "path", "triangle", "mesh", "python", "cython"],
    author="Michael Hogg",
    author_email="michael.christopher.hogg@gmail.com",
    url="https://github.com/mhogg/pygeodesic",
    download_url="https://github.com/mhogg/pygeodesic/releases",
    packages=find_packages(),
    include_package_data=True,
    package_data={"": ["README.md", "LICENSE.txt"], "pygeodesic": ["Examples\*"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Cython",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.9",
    ext_modules=ext_modules,
    cmdclass=cmdclass,
    setup_requires=["numpy"],
    install_requires=["numpy>=2,<3"],
)
