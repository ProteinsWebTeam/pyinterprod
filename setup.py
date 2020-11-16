#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import Extension, find_packages, setup

from pyinterprod import __version__


def get_requirements():
    filepath = os.path.join(os.path.dirname(__file__), "requirements.txt")

    with open(filepath) as fh:
        requirements = fh.read().splitlines()

    return requirements


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
    install_requires=get_requirements(),
    entry_points={
        "console_scripts": [
            # "ipr-clans = pyinterprod.clan:main",
            "ipr-pre-memdb = pyinterprod.cli:update_database",
            "ipr-memdb = pyinterprod.cli:run_member_db_update",
            "ipr-pronto = pyinterprod.cli:run_pronto_update",
            "ipr-uniprot = pyinterprod.cli:run_uniprot_update"
        ]
    }
)
