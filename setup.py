# ----------------------------------------------------------------------------
# Copyright (c) 2021-, Andersen Lab development team.
#
# Distributed under the terms of the XXX
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()


description = ("Freyja recovers relative lineage abundances from mixed "
               "SARS-CoV-2 samples")


setup(
    name="freyja",
    version="2021.10",
    packages=find_packages(),
    author="Joshua Levy",
    license='XXX',
    author_email="jolevy@scripps.edu",
    url="https://github.com/joshuailevy/freyja",
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points='''
        [console_scripts]
        freyja=freyja._cli:cli
        ''',
    package_data={
        'freyja': ['data/*', ],
        },
    install_requires=["click", "numpy", "pandas", "cvxpy"],
)
