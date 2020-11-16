# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
from datetime import datetime


@dataclass
class Method:
    accession: str
    sig_type: str
    name: str = None
    description: str = None
    abstract: str = None
    date: datetime = None


@dataclass
class Clan:
    accession: str
    name: str = None
    description: str = None
    members: list = field(default_factory=list)
