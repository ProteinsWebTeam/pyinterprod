# -*- coding: utf-8 -*-

from dataclasses import dataclass


@dataclass
class Method:
    accession: str
    sig_type: str
    name: str = None
    description: str = None
    abstract: str = None
