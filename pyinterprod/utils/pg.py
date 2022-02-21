import re
import time
from typing import Any, Iterator

from psycopg2.errors import DiskFull


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


def cluster(con, table: str, index: str, **kwargs):
    max_attempts = kwargs.get("max_attempts", 1)
    seconds = kwargs.get("seconds", 600)

    with con.cursor() as cur:
        num_attempts = 0

        while True:
            num_attempts += 1
            try:
                cur.execute(f"CLUSTER {table} USING {index}")
            except DiskFull:
                con.rollback()

                if num_attempts == max_attempts:
                    con.close()
                    raise

                time.sleep(seconds)
            else:
                con.commit()
                break
