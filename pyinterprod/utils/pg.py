import re
import time

from psycopg.errors import DiskFull


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


def get_primary_key(cur, table: str):
    cur.execute(
        """
        SELECT constraint_name
        FROM information_schema.table_constraints
        WHERE lower(table_name) = %s AND constraint_type = 'PRIMARY KEY'        
        """,
        [table.lower()]
    )
    row = cur.fetchone()
    return row[0] if row else None
