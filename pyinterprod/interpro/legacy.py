# -*- coding: utf-8 -*-

import time
from datetime import datetime, timedelta

import cx_Oracle

from pyinterprod.utils import email


def refresh_mviews(url: str, send_email: bool=True, notice: int=3600):
    if send_email:
        date = datetime.now() + timedelta(seconds=notice)
        email.send(
            to=[email.INTERPRO],
            subject="MV update scheduled",
            content=f"""\
Dear curators,

Legacy materialized views will start being refreshed \
on {date:%a %d %b %Y} at {date:%H:%M}.

You may continue to use Pronto but please do not perform \
any of the following actions:
    - integrate a member database signature;
    - unintegrate a member database signature;
    - delete an Interpro entry (unless it does not have any signature).
    
InterPro entries may be created (but do NOT integrate any signature), \
modified (name, short name, abstracts, GO terms, etc.), or checked/unchecked.

An email notification will be sent at the end of the update.

The InterPro Production Team
"""
        )

        time.sleep(notice)

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    # Relies on INTERPRO.ENTRY2METHOD and INTERPRO.MATCH
    cur.callproc("INTERPRO.REFRESH_MATCH_COUNTS.REFRESH")
    cur.close()
    con.close()

    if send_email:
        email.send(
            to=[email.INTERPRO],
            subject="MV update complete",
            content="""\
Dear curators,

Legacy materialized views are now updated. \
You may resume integrating or unintegrating signatures. Have fun!

The InterPro Production Team
"""
        )
