from .. import orautils
import cx_Oracle
import os


def execute_query(cur, query):
    results = []
    cur.execute(query)
    query_res = cur.fetchall()
    res_count = len(query_res)
    if res_count > 0:
        for row in query_res:
            if len(row) > 1:
                results.append("\t".join(row))
            else:
                results.append(str(row[0]))
        if len(results) > 1:
            results.append("".join(["\n", str(res_count), " rows selected."]))
    else:
        results.append("no rows selected.")
    return results


def generate_old_report(cur, dbcode, outdir):
    filename = os.path.join(outdir, f"oldSigStatsReport.{dbcode}.txt")
    lines = []

    sqlqueries = {
        "Total number of signatures:": "SELECT count(1) FROM interpro.method WHERE dbcode = '{}'",
        "Check how many of them skip_flagged and checked as Yes:": """SELECT m.method_ac, e.entry_ac
       FROM interpro.method m
       ,interpro.entry e
       ,interpro.entry2method e2m
       WHERE e2m.method_ac = m.method_ac
       AND m.skip_flag = 'Y'
       AND e2m.entry_ac = e.entry_ac
       AND e.checked = 'Y'
       AND m.dbcode = '{}'""",
        "Number of integrated signatures:": """SELECT count(1)
       FROM interpro.method m
       ,interpro.entry2method e2m
       WHERE m.method_ac = e2m.method_ac
       AND m.dbcode = '{}'""",
        "Unintegrated signatures:": """SELECT m.method_ac
       FROM interpro.method m
       WHERE m.dbcode = '{0}'
       MINUS
       SELECT e2m.method_ac
       FROM interpro.entry2method e2m
       ,interpro.method m
       WHERE m.method_ac = e2m.method_ac
       AND m.dbcode = '{0}'""",
        "Number of existing matches:": "SELECT count(1) FROM interpro.match WHERE dbcode = '{}'",
    }

    for key in sqlqueries:
        header = "\n" + key
        lines.append(header)
        query = sqlqueries[key].format(dbcode)
        lines += execute_query(cur, query)

    with open(filename, "w") as output:
        output.write("\n".join(lines))
        output.write("\n")


def get_signatures_to_delete(cur, dbcode):
    query = """SELECT method_ac
       FROM interpro.method
       WHERE dbcode = '{0}' 
       MINUS
       SELECT method_ac
       FROM interpro.method_stg
       WHERE dbcode = :1"""
    cur.execute(query)
    query_res = cur.fetchall()
    res_count = len(query_res)

    toprint = "Signatures To Be Deleted\n"
    if res_count > 0:
        toprint += f"{str(res_count)} rows selected."
    else:
        toprint += "no rows selected."

    return toprint


def generate_new_report(cur, dbcode, outdir):
    filename = os.path.join(outdir, f"newSigStatsReport.{dbcode}.txt")
    lines = []
    lines.append("Statistics for new data\n")
    sqlqueries = {
        "Signatures To Be Deleted": """SELECT method_ac
       FROM interpro.method
       WHERE dbcode = '{0}' 
       MINUS
       SELECT method_ac
       FROM interpro.method_stg
       WHERE dbcode = '{0}'""",
        "New Signatures": """SELECT method_ac
       FROM interpro.method_stg
       WHERE dbcode = '{0}'
       MINUS
       SELECT method_ac
       FROM interpro.method
       WHERE dbcode = '{0}'""",
        "Name Change\nMethod ac	Old name	New name": """SELECT m.method_ac, m.name, s.name
       FROM interpro.method_stg s
       ,interpro.method m
       WHERE m.dbcode = '{0}'
       AND m.dbcode = s.dbcode
       AND m.method_ac = s.method_ac
       AND m.name != s.name""",
        "Description Change\nMethod ac	Old name	New name": """SELECT m.method_ac, m.description, s.description
       FROM interpro.method_stg s
       ,interpro.method m
       WHERE m.dbcode = '{0}'
       AND m.dbcode = s.dbcode
       AND m.method_ac = s.method_ac
       AND m.description != s.description""",
        "Signature Type Change\nMethod ac        Old_type        New_type": """SELECT m.method_ac, m.sig_type, s.sig_type
       FROM interpro.method_stg s
       ,interpro.method m
       WHERE m.dbcode = '{0}'
       AND m.dbcode = s.dbcode
       AND m.method_ac = s.method_ac
       AND m.sig_type != s.sig_type""",
        "List of SignatureLess Entries as a Result of Deletion": """SELECT entry_ac,method_ac FROM (
       SELECT e2m.entry_ac
       ,e2m.method_ac
       ,CASE
       WHEN count(e2m.method_ac) over (partition by e2m.entry_ac) = a.n_db THEN 'DB Only'
       ELSE 'Other DB as well'
       END ding
       FROM interpro.method m
       ,interpro.entry2method e2m
       ,(SELECT e2m.entry_ac
       ,count(e2m.method_ac) n_db
       FROM interpro.entry e
       ,interpro.entry2method e2m
       WHERE e2m.method_ac IN (SELECT method_ac
       FROM interpro.method
       WHERE dbcode = '{0}'
       MINUS
       SELECT method_ac
       FROM interpro.method_stg
       WHERE dbcode = '{0}')
       AND e.entry_ac = e2m.entry_ac
       GROUP BY e2m.entry_ac ) a
       WHERE m.method_ac = e2m.method_ac
       AND e2m.entry_ac = a.entry_ac
       ) a WHERE a.ding = 'DB Only'""",
        "Number Matches Affected As a Result Of Deletion": """SELECT count(1)
       FROM interpro.match
       WHERE method_ac IN (SELECT method_ac
       FROM interpro.method
       WHERE dbcode = '{0}'
       MINUS
       SELECT method_ac
       FROM interpro.method_stg
       WHERE dbcode = '{0}')""",
    }
    todelete = 0

    for key in sqlqueries:
        header = "\n" + key
        lines.append(header)
        query = sqlqueries[key].format(dbcode)
        result = execute_query(cur, query)
        lines += result
        print(key, result)
        if key == "Signatures To Be Deleted":
            if result[0] != "no rows selected.":
                todelete = result[0].split(" ")[0]

    with open(filename, "w") as output:
        output.write("\n".join(lines))
        output.write("\n")

    return todelete


def generate_report(user: str, dsn: str, outdir: str, memberdb: list):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    for member in memberdb:
        if member["dbcode"]:
            dbcode = member["dbcode"]
            print(f"\nprocessing dbcode: {dbcode}")
            generate_old_report(cur, dbcode, outdir)
            deleted = generate_new_report(cur, dbcode, outdir)
            if deleted > 0:
                print(f"{deleted} dead signatures found")
                # delete_dead_signatures(dbcode)
        else:
            print(f"Error dbcode not found for {member['name']}")

    cur.close()
    con.close()
