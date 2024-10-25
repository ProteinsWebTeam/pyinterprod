import tarfile
import uuid

import oracledb

from pyinterprod import logger
from pyinterprod.utils import Table
from pyinterprod.utils.oracle import drop_table, get_partitions


def load_matches(uri: str, databases: dict[str, str]):
    con = oracledb.connect(uri)
    cur = con.cursor()

    partitions = {}
    for p in get_partitions(cur, "INTERPRO", "TOAD_MATCH"):
        name = p["name"]
        dbcode = p["value"][1:-1]  # 'X' -> X
        partitions[dbcode] = name

    for dbshort, filepath in databases.items():
        sql = "SELECT DBCODE FROM INTERPRO.CV_DATABASE WHERE DBSHORT = :1"
        cur.execute(sql, [dbshort.upper()])
        row = cur.fetchone()
        if row is None:
            cur.close()
            con.close()
            raise ValueError(f"No database found for {dbshort}")

        dbcode, = row

        try:
            partition = partitions[dbcode]
        except KeyError:
            cur.close()
            con.close()
            err = f"No partition in TOAD_MATCH for database {dbshort}"
            raise KeyError(err)

        load_database_matches(cur, dbcode, partition, filepath)

    cur.close()
    con.close()


def load_database_matches(cur: oracledb.Cursor, dbcode: str, partition: str,
                          filepath: str):
    # Insert "raw" matches (as provided by DeepMind)
    drop_table(cur, "INTERPRO.TOAD_MATCH_NEW", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.TOAD_MATCH_NEW NOLOGGING
        AS SELECT * FROM INTERPRO.TOAD_MATCH WHERE 1 = 0
        """
    )

    query = """
        INSERT /*+ APPEND */ 
        INTO INTERPRO.TOAD_MATCH_NEW 
        VALUES (:1, :2, :3, :4, :5, :6, :7)
    """
    with Table(con=cur.connection, query=query, autocommit=True) as table:
        for uniprot_acc, method_acc, frags, score in iter_matches(filepath):
            if len(frags) > 1:
                group_uuid = str(uuid.uuid4())
            else:
                group_uuid = None

            for pos_from, pos_to in frags:
                table.insert((
                    uniprot_acc,
                    method_acc,
                    dbcode,
                    pos_from,
                    pos_to,
                    group_uuid,
                    score
                ))

    # Filter out obsolete matches (involving deleted proteins or signatures)
    drop_table(cur, "INTERPRO.TOAD_MATCH_TMP", purge=True)
    cur.execute(
        """
        CREATE TABLE INTERPRO.TOAD_MATCH_TMP NOLOGGING
        AS SELECT * FROM INTERPRO.TOAD_MATCH WHERE 1 = 0
        """
    )
    cur.execute(
        """
        INSERT /*+ APPEND */ INTO INTERPRO.TOAD_MATCH_TMP
        SELECT *
        FROM (
            SELECT T.PROTEIN_AC,
                   T.METHOD_AC,
                   T.DBCODE
                   T.POS_FROM,
                   CASE WHEN M.POS_TO <= P.LEN
                        THEN M.POS_TO 
                        ELSE P.LEN 
                        END POS_TO,
                   T.UUID,
                   T.SCORE 
            FROM INTERPRO.TOAD_MATCH_NEW T
            INNER JOIN INTERPRO.PROTEIN P ON T.PROTEIN_AC = P.PROTEIN_AC
            INNER JOINT INTERPRO.METHOD M ON T.METHOD_AC = M.METHOD_AC
        )
        WHERE POS_FROM <= POS_TO
        """,
    )
    cur.connection.commit()

    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT CK_TOAD_MATCH_TMP$FROM 
        CHECK (POS_FROM >= 1)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT CK_TOAD_MATCH_TMP$TO 
        CHECK (POS_FROM <= POS_TO)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT PK_TOAD_MATCH_TMP
        PRIMARY KEY (PROTEIN_AC, METHOD_AC, DBCODE, POS_FROM, POS_TO)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT FK_TOAD_MATCH_TMP$PROTEIN 
        FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT FK_TOAD_MATCH_TMP$PROTEIN 
        FOREIGN KEY (METHOD_AC) REFERENCES INTERPRO.METHOD (METHOD_AC)
        """
    )
    cur.execute(
        """
        ALTER TABLE INTERPRO.TOAD_MATCH_TMP
        ADD CONSTRAINT FK_TOAD_MATCH_TMP$DBCODE 
        FOREIGN KEY (DBCODE) REFERENCES INTERPRO.CV_DATABASE (DBCODE)
        """
    )

    cur.execute(
        f"""
        ALTER TABLE INTERPRO.TOAD_MATCH
        EXCHANGE PARTITION ({partition})
        WITH TABLE INTERPRO.TOAD_MATCH_TMP
        """
    )


def iter_matches(filepath: str):
    """Iterate TOAD inferences
    """
    with tarfile.open(filepath, mode="r") as tar:
        for member in tar:
            if member.name.endswith(".tsv"):
                br = tar.extractfile(member)
                lines = br.read().decode("utf-8").splitlines(keepends=False)

                # First line is a header
                for line in lines[1:]:
                    values = line.split("\t")
                    if len(values) == 5:
                        # No discontinuous domains
                        uniprot_acc, signature_acc, start, end, score = values
                        yield (uniprot_acc, signature_acc,
                               [(int(start), int(end))], float(score))
                    elif len(values) == 23:
                        # Discontinuous domains (up to ten fragments)
                        uniprot_acc, signature_acc = values[:2]
                        fragments = []
                        for i in range(10):
                            start = values[2+i*2]
                            end = values[3+i*2]
                            if start == "NULL":
                                break
                            else:
                                fragments.append((int(start), int(end)))

                        score = float(values[-1])
                        fragments.sort()
                        yield uniprot_acc, signature_acc, fragments, score
                    else:
                        err = f"Unexpected number of columns: {values}"
                        raise ValueError(err)
