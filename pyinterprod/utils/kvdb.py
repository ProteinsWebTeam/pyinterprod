# -*- coding: utf-8 -*-

import pickle
import sqlite3


class KVdb:
    def __init__(self, filepath: str, writeback: bool = False):
        self.filepath = filepath
        self.writeback = writeback
        self.con = sqlite3.connect(self.filepath)
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS data (
                id TEXT PRIMARY KEY NOT NULL,
                value TEXT NOT NULL
            )
            """
        )
        self.cache = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __len__(self) -> int:
        return self.con.execute("SELECT COUNT(*) FROM data").fetchone()[0]

    def __delitem__(self, key):
        try:
            del self.cache[key]
        except KeyError:
            pass
        finally:
            self.con.execute("DELETE FROM data WHERE id = ?", (key,))
            self.con.commit()

    def __getitem__(self, key):
        try:
            return self.cache[key]
        except KeyError:
            pass

        sql = "SELECT value FROM data WHERE id = ?"
        row = self.con.execute(sql, (key,)).fetchone()

        if row is None:
            raise KeyError(key)

        value = pickle.loads(row[0])
        if self.writeback:
            self.cache[key] = value

        return value

    def __setitem__(self, key, value):
        if self.writeback:
            self.cache[key] = value
        else:
            sql = "INSERT OR REPLACE INTO data (id, value) VALUES (?, ?)"
            self.con.execute(sql, (key, pickle.dumps(value)))
            self.con.commit()

    def __iter__(self):
        return self.keys()

    def keys(self):
        for key, in self.con.execute("SELECT id FROM data ORDER BY id"):
            yield key

    def values(self):
        for obj, in self.con.execute("SELECT value FROM data ORDER BY id"):
            yield pickle.loads(obj)

    def items(self):
        for row in self.con.execute("SELECT id, value FROM data ORDER BY id"):
            yield row[0], pickle.loads(row[1])

    def close(self):
        if self.con is None:
            return

        self.sync()
        self.con.close()
        self.con = None

    def sync(self):
        if not self.cache:
            return

        sql = "INSERT OR REPLACE INTO data (id, value) VALUES (?, ?)"
        self.con.executemany(sql, ((key, pickle.dumps(value))
                                   for key, value in self.cache.items()))
        self.con.commit()
        self.cache = {}
