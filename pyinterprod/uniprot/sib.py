import oracledb

from pyinterprod import logger
from pyinterprod.utils import email


def export(url: str, emails: dict):
    logger.info("exporting data")
    con = oracledb.connect(url)
    cur = con.cursor()

    """
    The SIB_EXPORT package belongs to INTERPRODP and exports the tables:
    * INTERPRO.DB_VERSION
    * INTERPRO.ENTRY
    * INTERPRO.ENTRY2ENTRY
    * INTERPRO.ENTRY2METHOD
    * INTERPRO.MATCH
    * INTERPRO.METHOD
    * INTERPRO.PROTEIN
    * INTERPRO.XREF_CONDENSED
    """
    cur.callproc("SIB_EXPORT.EXP_SIB_DATA")

    logger.info("data successfully exported")

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'u'")
    version, = cur.fetchone()
    cur.close()
    con.close()

    email.send(
        info=emails,
        to=["interpro"],
        subject=f"Protein update {version}: data for SIB ready",
        content=f"""\
The data required by SIB was successfully exported.

Please, archive the dump on the FTP, and inform SIB that the archive is ready.

Recipients
----------
To: {emails['sib']}

Subject
-------
InterPro data for UniProt private release available

Body 
----
Dear Swiss-Prot team,

The interpro.tar.gz archive for UniProt private release {version} \
is available at ftp://ftp-private.ebi.ac.uk/interpro/

Kind regards,
The InterPro Production Team
"""
    )
