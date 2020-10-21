# -*- coding: utf-8 -*-

from dataclasses import dataclass
from datetime import datetime


@dataclass
class Method:
    accession: str
    sig_type: str
    name: str = None
    description: str = None
    abstract: str = None
    date: datetime = None
