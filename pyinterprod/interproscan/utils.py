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
                    f"  - {p['name']:<30} ({analysis_id:>5}): delete data",
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
                    f"  - {p['name']:<30} ({analysis_id:>5}): delete data",
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
                        f"  - {p['name']:<30} ({analysis_id:>5}): delete jobs/data > {max_upi}",
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
                    f"  - {p['name']:<30} ({analysis_id:>5}): delete jobs/data",
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

    for table in tables:
        for index in oracle.get_indexes(cur, "IPRSCAN", table):
            if index["is_unusable"]:
                logger.info(f"rebuilding {index['name']}")

                try:
                    oracle.rebuild_index(cur, index["name"])
                except Exception as exc:
                    cur.close()
                    con.close()
                    raise

    cur.close()
    con.close()
