from typing import Union


class Table:
    def __init__(self, con, query: str, autocommit: bool = False,
                 buffer_size: int = 100000, depends_on=None):
        self.con = con
        self.cur = con.cursor()
        self.query = query
        self.autocommit = autocommit
        self.buffer_size = buffer_size
        self.depends_on = depends_on
        self.rows = []
        self.count = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def _execute(self, record: Union[dict, list, tuple]):
        self.rows.append(record)
        self.count += 1

        if len(self.rows) == self.buffer_size:
            self.flush()

    def insert(self, record: Union[dict, list, tuple]):
        self._execute(record)

    def update(self, record: Union[dict, list, tuple]):
        self._execute(record)

    def delete(self, record: Union[dict, list, tuple]):
        self._execute(record)

    def flush(self):
        if not self.rows:
            return
        elif self.depends_on:
            self.depends_on.flush()

        self.cur.executemany(self.query, self.rows)
        self.rows = []

        if self.autocommit:
            self.con.commit()

    def close(self):
        if self.con is not None:
            self.flush()
            self.cur.close()
            self.con = None
