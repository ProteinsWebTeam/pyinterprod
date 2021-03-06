#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import Extension, find_packages, setup

from pyinterprod import __version__


setup(
    name="pyinterprod",
    version=__version__,
    description="",
    long_description="",
    packages=find_packages(),
    ext_modules=[
        Extension(
            "pyinterprod.proteinupdate.sprot",
            ["pyinterprod/proteinupdate/sprotmodule.c"],
            extra_link_args=["-lsqlite3"],
        )
    ],
    install_requires=["cx-Oracle>=7.1", "mundone>=0.2.0"],
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "ipr-clans = pyinterprod.clan:main",
            "ipr-pronto = pyinterprod.pronto:main",
            "ipr-uniprot = pyinterprod.proteinupdate:main",
            "ipr-memberdb = pyinterprod.memberdbupdate:main",
        ]
    },
)
