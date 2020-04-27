# -*- coding: utf-8 -*-

import re
from typing import Any, Iterator

from psycopg2.errors import UndefinedObject


class CsvIO:
    def __init__(self, records: Iterator, sep: str="\t", null="\\N"):
        self.records = records
        self.buffer = ''
        self.sep = sep
        self.null = null

    def readline(self):
        try:
            record = next(self.records)
        except StopIteration:
            return ''
        else:
            return self.sep.join(map(self.format, record)) + '\n'

    def read(self, size: int=-1):
        output = ""
        if size is None or size < 0:
            while True:
                line = self.readline()
                if line:
                    output += line
                else:
                    break
        else:
            output = self.buffer
            while len(output) < size:
                line = self.readline()
                if line:
                    output += line
                else:
                    break

            self.buffer = output[size:]
            output = output[:size]
        return output

    def format(self, value: Any):
        if value is None:
            return self.null
        else:
            return str(value).replace('\n', '\\n')


def drop_index(con, index: str):
    with con.cursor() as cur:
        try:
            cur.execute(f"DROP INDEX {index}")
        except UndefinedObject as exc:
            con.rollback()
        else:
            con.commit()


def url2dict(url: str) -> dict:
    m = re.match(r'([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', url)

    if m is None:
        raise RuntimeError(f"invalid connection string: {url}")

    return dict(
        user=m.group(1),
        password=m.group(2),
        host=m.group(3),
        port=int(m.group(4)),
        dbname=m.group(5)
    )
