# -*- coding: utf-8 -*-

import os
from tempfile import mkstemp
from urllib.request import urlopen


def download(url: str) -> str:
    fd, dst = mkstemp()
    os.close(fd)

    with urlopen(url) as res, open(dst, "wb") as fh:
        fh.write(res.read())

    return dst
