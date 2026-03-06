"""
Setup file for parp14_tools package.

Installation:
    cd parp14_tools
    pip install -e .

This installs the package in "editable" mode, so changes to the
source files take effect immediately without reinstalling.
"""

from setuptools import setup, find_packages

setup(
    name='parp14_tools',
    version='0.1.0',
    description='Analysis and visualization tools for CALVADOS simulations',
    author='Your Name',
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        'numpy',
        'matplotlib',
        'mdtraj',
        'MDAnalysis',
        'pandas',
        'scipy',
    ],
)
