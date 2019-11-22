# -*- coding: utf-8 -*-
#
# This file is part of REANA.
# Copyright (C) 2017, 2018, 2019 CERN.
#
# REANA is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""REANA client"""

from __future__ import absolute_import, print_function

import os
import re

from setuptools import find_packages, setup

readme = open('README.rst').read()

setup_requires = [
]

install_requires = [
    'matplotlib>=3.1.2',
    'mpmath>=1.1.0'
]

packages = find_packages()


# Get the version string. Cannot be done with import!
with open(os.path.join('reana_client', 'version.py'), 'rt') as f:
    version = re.search(
        '__version__\s*=\s*"(?P<version>.*)"\n',
        f.read()
    ).group('version')

setup(
    name='mieze-simulation',
    version=version,
    description=__doc__,
    long_description=readme,
    author='FRM2',
    author_email='daniel.prelipcean@frm2.tum.de',
    url='https://github.com/dprelipcean/mieze-simulation',
    zip_safe=False,
    include_package_data=True,
    install_requires=install_requires,
    setup_requires=setup_requires,
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python',
    ],
)