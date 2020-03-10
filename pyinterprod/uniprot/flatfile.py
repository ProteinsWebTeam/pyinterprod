# -*- coding: utf-8 -*-

from typing import Tuple

from . import sprot


def load(swissp: str, trembl: str, database: str) -> Tuple[int, int]:
    swissp_cnt = sprot.load(swissp, database, "protein")
    trembl_cnt = sprot.load(trembl, database, "protein")

    return swissp_cnt, trembl_cnt
