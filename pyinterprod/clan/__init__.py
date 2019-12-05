#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from typing import Sequence

import cx_Oracle

from pyinterprod import logger, orautils
from pyinterprod.clan import cdd, pfam, panther, pirsf


def create_tables(cur: cx_Oracle.Cursor):
    orautils.drop_table(cur, "INTERPRO", "CLAN2MEMBER", purge=True)
    orautils.drop_table(cur, "INTERPRO", "CLAN", purge=True)

    cur.execute(
        """
        CREATE TABLE INTERPRO.CLAN
        (
            CLAN_AC VARCHAR2(25) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            NAME VARCHAR2(100) DEFAULT NULL,
            DESCRIPTION VARCHAR2(4000) DEFAULT NULL,
            CONSTRAINT PK_CLAN PRIMARY KEY (CLAN_AC),
            CONSTRAINT FK_CLAN$DBCODE
              FOREIGN KEY (DBCODE)
              REFERENCES INTERPRO.CV_DATABASE (DBCODE)
              ON DELETE CASCADE 
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE INTERPRO.CLAN2MEMBER
        (
            CLAN_AC VARCHAR2(25) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            SCORE FLOAT NOT NULL,
            CONSTRAINT PK_CLAN2MEMBER PRIMARY KEY (CLAN_AC, METHOD_AC),
            CONSTRAINT UQ_CLAN2MEMBER$METHOD_AC UNIQUE (METHOD_AC),
            CONSTRAINT FK_CLAN2MEMBER$CLAN_AC
              FOREIGN KEY (CLAN_AC)
              REFERENCES INTERPRO.CLAN (CLAN_AC),
            CONSTRAINT FK_CLAN2MEMBER$METHOD_AC
              FOREIGN KEY (METHOD_AC)
              REFERENCES INTERPRO.METHOD (METHOD_AC)
              ON DELETE CASCADE
        )
        """
    )


def _insert_clans(con: cx_Oracle.Connection, clans: Sequence[dict], dbcode: str):
    cur = con.cursor()
    cur.execute("SELECT METHOD_AC FROM INTERPRO.METHOD WHERE DBCODE=:1", (dbcode,))
    members = {row[0] for row in cur}
    cur.close()

    t1 = orautils.TablePopulator(con, "INSERT INTO INTERPRO.CLAN "
                                      "VALUES (:1, :2, :3, :4)")
    t2 = orautils.TablePopulator(con, "INSERT INTO INTERPRO.CLAN2MEMBER "
                                      "VALUES (:1, :2, :3)", depends_on=t1)

    for c in clans:
        clan_members = []
        for member in c["members"]:
            if member["accession"] in members:
                clan_members.append(member)
            else:
                logger.warning(f"{member['accession']} is not "
                               f"a valid signature accession")

        if not clan_members:
            continue

        t1.insert((c["accession"], dbcode, c["name"], c["description"]))
        for member in clan_members:
            t2.insert(
                (c["accession"], member["accession"], member["score"]))

    t1.close()
    t2.close()


def _main():
    dsn = ""
    user = ""
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    create_tables(cur)
    cur.close()

    pfam_url = ""
    logger.info("pfam")
    _insert_clans(con, cdd.get_clans(), 'J')
    logger.info("cdd")
    _insert_clans(con, pfam.get_clans(pfam_url), 'H')
    logger.info("panther")
    _insert_clans(con, panther.get_clans(con), 'V')
    logger.info("pirsf")
    _insert_clans(con, pirsf.get_clans(), 'U')
    con.commit()
    con.close()


if __name__ == '__main__':
    _main()



# def main():
#     parser = argparse.ArgumentParser(description="Pronto schema update")
#     parser.add_argument("config", metavar="CONFIG.JSON",
#                         help="config JSON file")
#     parser.add_argument("-t", "--tasks", nargs="*",
#                         help="tasks to run (default: all)")
#     parser.add_argument("-o", "--output", default=_REPORT,
#                         help=f"output report for curators (default: {_REPORT})")
#     args = parser.parse_args()
#
#     success = run(args.config, report=args.output, tasks=args.tasks)
#     sys.exit(0 if success else 1)
