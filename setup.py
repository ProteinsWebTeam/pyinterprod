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
        Extension(name="pyinterprod.uniprot.sprot",
                  sources=["pyinterprod/uniprot/sprotmodule.c"],
                  extra_link_args=["-lsqlite3"])
    ],
    install_requires=["cx-Oracle>=7.3", "mundone>=0.4.0", "psycopg2>=2.8.4"],
    entry_points={
        "console_scripts": [
            # "ipr-clans = pyinterprod.clan:main",
            "ipr-pronto = pyinterprod.cli:run_pronto_update",
            "ipr-uniprot = pyinterprod.cli:run_uniprot_update"
        ]
    }
)
