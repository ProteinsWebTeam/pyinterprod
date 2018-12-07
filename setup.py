#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import Extension, find_packages, setup

setup(
    name="pyinterprod",
    version="1.0",
    description="",
    long_description="",
    packages=find_packages(),
    ext_modules=[
        Extension("pyinterprod.proteinupdate.sprot",
                  ["pyinterprod/proteinupdate/sprotmodule.c"],
                  extra_link_args=["-lsqlite3"])
    ],
    install_requires=["cx-Oracle>=6.0"],
    zip_safe=False,
)
