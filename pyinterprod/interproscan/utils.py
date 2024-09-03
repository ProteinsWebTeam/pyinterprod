from concurrent.futures import ThreadPoolExecutor, as_completed

import oracledb

from pyinterprod import logger
from pyinterprod.utils import oracle
from .analyses import get_analyses


def clean_tables(uri: str, analysis_ids: list[int] | None = None):
    con = oracledb.connect(uri)
    cur = con.cursor()

    analyses = get_analyses(cur)

    table2analyses = {}
    for analysis_id, analysis in analyses.items():
        table = analysis["tables"]["matches"]
        try:
            table2analyses[table].append(analysis_id)
        except KeyError:
            table2analyses[table] = [analysis_id]

        table = analysis["tables"]["sites"]
        if table:
            try:
                table2analyses[table].add(analysis_id)
            except KeyError:
                table2analyses[table] = {analysis_id}

    actions = []
    for table in sorted(table2analyses):
        table = table.upper()
        for p in oracle.get_partitions(cur, "IPRSCAN", table):
            if p["value"] == "DEFAULT":
                continue

            analysis_id = int(p["value"])

            if analysis_ids and analysis_id not in analysis_ids:
                continue

            if analysis_id not in analyses:
                # Obsolete analysis: remove data
                actions.append((
                    f"  - {p['name']:<30}: drop partition",
                    [(
                        f"ALTER TABLE {table} DROP PARTITION {p['name']}", []
                    )]
                ))
                continue

            analysis = analyses[analysis_id]
            max_upi = analysis["max_upi"]

            if analysis_id not in table2analyses[table]:
                # Obsolete analysis: remove data
                actions.append((
                    f"  - {p['name']:<30}: drop partition",
                    [(
                        f"ALTER TABLE {table} DROP PARTITION {p['name']}", []
                    )]
                ))
            elif max_upi:
                cur.execute(
                    """
                    SELECT COUNT(*) 
                    FROM IPRSCAN.ANALYSIS_JOBS 
                    WHERE ANALYSIS_ID = :1
                      AND UPI_FROM > :2
                    """,
                    [analysis_id, max_upi]
                )
                cnt, = cur.fetchone()

                if cnt > 0:
                    # Delete jobs after the max UPI
                    actions.append((
                        f"  - {p['name']:<30}: delete jobs > {max_upi}",
                        [(
                            """
                            DELETE FROM IPRSCAN.ANALYSIS_JOBS
                            WHERE ANALYSIS_ID = :1
                            AND UPI_FROM > :2
                            """,
                            [analysis_id, max_upi]
                        ), (
                            f"""
                            DELETE FROM {table} PARTITION ({p['name']})
                            WHERE UPI_FROM > :1
                            """,
                            [max_upi]
                        )]
                    ))
            else:
                # No max UPI: remove data
                actions.append((
                    f"  - {p['name']:<30}: delete jobs",
                    [(
                        f"DELETE FROM IPRSCAN.ANALYSIS_JOBS "
                        f"WHERE ANALYSIS_ID = :1",
                        [analysis_id]
                    ), (
                        f"ALTER TABLE {table} TRUNCATE PARTITION {p['name']}",
                        []
                    )]
                ))

    if actions:
        print("The following actions will be performed:")
        for descr, queries in actions:
            print(descr)

        if input("Proceed? [y/N] ").lower().strip() == "y":
            for descr, queries in actions:
                for sql, params in queries:
                    cur.execute(sql, params)

            con.commit()
        else:
            print("Canceled")

    cur.close()
    con.close()


def rebuild_indexes(uri: str, analysis_ids: list[int] | None = None):
    con = oracledb.connect(uri)
    cur = con.cursor()

    analyses = get_analyses(cur)

    tables = set()
    for analysis_id, analysis in analyses.items():
        if analysis_ids and analysis_id not in analysis_ids:
            continue

        tables.add(analysis["tables"]["matches"])

        if analysis["tables"]["sites"]:
            tables.add(analysis["tables"]["sites"])

    to_rebuild = set()
    for table in tables:
        for index in oracle.get_indexes(cur, "IPRSCAN", table):
            if index["is_unusable"]:
                to_rebuild.add(index["name"])

    cur.close()
    con.close()

    errors = 0
    with ThreadPoolExecutor(max_workers=8) as executor:
        fs = {}

        for index in to_rebuild:
            f = executor.submit(_rebuild_index, uri, index)
            fs[f] = index

        for f in as_completed(fs):
            index = fs[f]

            try:
                f.result()
            except Exception as exc:
                logger.error(f"{index} rebuild failed: {exc}")
                errors += 1
            else:
                logger.info(f"{index} rebuilt")

    if errors > 0:
        raise RuntimeError(f"{errors} errors occurred")


def _rebuild_index(uri: str, name: str):
    logger.info(f"rebuilding {name}")
    con = oracledb.connect(uri)
    cur = con.cursor()
    oracle.rebuild_index(cur, name)
    cur.close()
    con.close()
