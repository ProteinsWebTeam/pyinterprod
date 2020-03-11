from .. import orautils
import cx_Oracle


def chunks(l, n):
    """Yield chunks of size n from iterable.

    Args:
        l (Iterable): The iterable to chunk
        n (int): Maximum number of items in chunk

    Yields:
        Iterable: Tuples of length n or less (final bucket will have only the remaining items)

    """
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i : i + n]


def populate_method_stg(user: str, dsn: str, dat_file: str, memberdb: list):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    data_list = list()

    try:
        with open(dat_file, "r") as fd:
            for line in fd:
                line = line.strip(":\n").split("::")
                line[5] = line[5].replace("|", "")
                if line[5] == "":
                    line[5] = None
                print(line)
                data_list.append(line)
    except FileNotFoundError as e:
        print(e)

    query = "DELETE FROM INTERPRO.METHOD_STG"
    cur.execute(query)
    con.commit()

    data_list_chunk = list(chunks(data_list, 1000))
    query = "INSERT INTO INTERPRO.METHOD_STG VALUES(:1,:2,:3,:4,:5,:6)"
    for chunk in data_list_chunk:
        cur.executemany(query, chunk)
    con.commit()

    cur.close()
    con.close()

    for member in memberdb:
        if member["incremental"] == True:
            incremental_update(user, dsn, member["dbcode"])


def update_method(user: str, dsn: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    query = """ MERGE
                INTO INTERPRO.METHOD M
                USING INTERPRO.METHOD_STG S
                ON (M.METHOD_AC = S.METHOD_AC)
                WHEN MATCHED THEN UPDATE SET M.ABSTRACT = S.ABSTRACT, M.DESCRIPTION=S.DESCRIPTION, M.SIG_TYPE=S.SIG_TYPE, M.NAME=S.NAME
                WHEN NOT MATCHED
                    THEN INSERT (M.METHOD_AC,M.NAME,M.DBCODE,M.DESCRIPTION,M.SIG_TYPE,M.METHOD_DATE,M.CANDIDATE,M.ABSTRACT)
                    VALUES (S.METHOD_AC,S.NAME,S.DBCODE,S.DESCRIPTION,
                            S.SIG_TYPE,SYSDATE,'Y',S.ABSTRACT)
      """

    cur.execute(query)
    con.commit()

    cur.close()
    con.close()


def incremental_update(user: str, dsn: str, dbcode: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    query = """
        INSERT INTO INTERPRO.METHOD_STG(METHOD_AC,DBCODE,NAME,DESCRIPTION,SIG_TYPE,ABSTRACT)
        SELECT METHOD_AC,DBCODE,NAME,DESCRIPTION,SIG_TYPE,ABSTRACT
        FROM INTERPRO.METHOD
        WHERE DBCODE=:1 AND METHOD_AC NOT IN (SELECT METHOD_AC FROM INTERPRO.METHOD_STG);
    """

    cur.execute(query, (dbcode,))
    con.commit()

    cur.close()
    con.close()


def get_dbcodes_memberdb(user: str, dsn: str, memberdb: list):
    memberwithdbcode = list()
    for member in memberdb:
        dbcode = get_dbcode(user, dsn, member["name"])
        if dbcode:
            member["dbcode"] = dbcode
        memberwithdbcode.append(member)
    return memberwithdbcode


def get_dbcode(user: str, dsn: str, memberdb: str):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    query_dbcode = """ SELECT DBCODE FROM INTERPRO.CV_DATABASE WHERE DBNAME=:1"""

    cur.execute(query_dbcode, (memberdb,))
    result = cur.fetchone()

    cur.close()
    con.close()

    if result:
        return result[0]
    else:
        print(f"Member database not found {memberdb}")
        return None


def update_db_version(user: str, dsn: str, memberdb: list):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    query_count = """ SELECT COUNT(*) FROM INTERPRO.METHOD WHERE DBCODE=:1"""

    for member in memberdb:
        if member["dbcode"]:
            dbcode = member["dbcode"]
            cur.execute(query_count, (dbcode,))
            countmethod_result = cur.fetchone()
            countmethod = 0
            if countmethod_result:
                countmethod = countmethod_result[0]
            try:
                query = "INSERT INTO INTERPRO.DB_VERSION(DBCODE,VERSION,ENTRY_COUNT,FILE_DATE) VALUES(:1,:2,:3,:4)"
                cur.execute(
                    query,
                    (dbcode, member["version"], countmethod, member["release-date"]),
                )
                print(f"{dbcode} inserted")
            except Exception as e:
                print("Couldn't insert, trying to update", e)
                query = "UPDATE INTERPRO.DB_VERSION SET ENTRY_COUNT=:1, VERSION=:2 , FILE_DATE=:3 ,LOAD_DATE=SYSDATE WHERE DBCODE=:4"
                cur.execute(
                    query,
                    (countmethod, member["version"], member["release-date"], dbcode),
                )

    con.commit()

    cur.close()
    con.close()


def update_iprscan2dbcode(user: str, dsn: str, memberdb: list):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    for member in memberdb:
        try:
            query = "UPDATE INTERPRO.IPRSCAN2DBCODE SET IPRSCAN_SIG_LIB_REL_ID = (SELECT MAX(ID) FROM IPRSCAN.MV_SIGNATURE_LIBRARY_RELEASE WHERE LIBRARY=:1) WHERE DBCODE=:2"
            cur.execute(query, (member["name"], member["dbcode"]))
            if cur.rowcount == 0:
                raise "Row not found"
            print(f"\nIPRSCAN2DBCODE for {member['name']} updated")
        except:
            print("\nCouldn't update, trying to insert")
            query = "INSERT INTO INTERPRO.IPRSCAN2DBCODE VALUES((SELECT MAX(ID) FROM IPRSCAN.MV_SIGNATURE_LIBRARY_RELEASE WHERE LIBRARY=:name),:dbcode,:evidence,(SELECT MAX(ID) FROM IPRSCAN.MV_SIGNATURE_LIBRARY_RELEASE WHERE LIBRARY=:name))"
            cur.execute(
                query,
                name=member["name"],
                dbcode=member["dbcode"],
                evidence=member["evidence"],
            )
            print(f"\n{member['name']} inserted")

    con.commit()

    cur.close()
    con.close()
