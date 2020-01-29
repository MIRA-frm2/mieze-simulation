# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""MIEZE simulation."""

from __future__ import absolute_import, print_function

from setuptools import find_packages, setup

readme = open('README.rst').read()

setup_requires = [
]

install_requires = [
    'matplotlib>=3.1.2',
    'mpmath>=1.1.0',
    'pandas>=0.25.3'
]

packages = find_packages()


setup(
    name='mieze-simulation',
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
