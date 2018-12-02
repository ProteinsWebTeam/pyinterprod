#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, Extension

setup(
    name="pyinterprod",
    version="1.0",
    description="",
    long_description="",
    ext_modules=[
        Extension("pyinterprod.proteinupdate.sprot",
                  ["pyinterprod/proteinupdate/sprotmodule.c"],
                  extra_link_args=["-lsqlite3"])
    ],
    zip_safe=False
)
