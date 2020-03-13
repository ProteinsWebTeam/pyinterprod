# -*- coding: utf-8 -*-

import mimetypes
import os
from email.message import EmailMessage
from getpass import getuser
from smtplib import SMTP
from typing import Sequence


# Email addresses constant
# AA = "automated_annotation@ebi.ac.uk"
# AA_DEV = "aa_dev@ebi.ac.uk"
# INTERPRO = "interpro-team@ebi.ac.uk"
# UNIPROT_DB = "uniprot-database@ebi.ac.uk"
# UNIPROT_PROD = "uniprot-prod@ebi.ac.uk"
# UNIRULE = "unirule@ebi.ac.uk"
# SIB = "uniprot-prod@isb-sib.ch"
AA = "mblum@ebi.ac.uk"
AA_DEV = "mblum@ebi.ac.uk"
INTERPRO = "mblum@ebi.ac.uk"
UNIPROT_DB = "mblum@ebi.ac.uk"
UNIPROT_PROD = "mblum@ebi.ac.uk"
UNIRULE = "mblum@ebi.ac.uk"
SIB = "mblum@ebi.ac.uk"

# Email server info
_EBI_EMAIL_SERVER = "outgoing.ebi.ac.uk"
_EBI_EMAIL_PORT = 587
_EBI_EMAIL = "@ebi.ac.uk"


def send(to: Sequence[str], subject: str, content: str, **kwargs):
    cc_addrs = kwargs.get("cc")
    bcc_addrs = kwargs.get("bcc")
    attachments = kwargs.get("attachments")
    sender = getuser() + _EBI_EMAIL

    msg = EmailMessage()
    msg.set_content(content)
    msg["Sender"] = sender
    msg["To"] = set(to)

    if cc_addrs:
        msg["Cc"] = ','.join(set(cc_addrs))

    if bcc_addrs:
        msg["Bcc"] = ','.join(set(bcc_addrs))

    msg["Subject"] = subject

    if attachments:
        for path in attachments:
            ctype, encoding = mimetypes.guess_type(path)

            if ctype is None or encoding is not None:
                # Use a generic type.
                ctype = "application/octet-stream"

            maintype, subtype = ctype.split('/', 1)

            with open(path, "rb") as fh:
                msg.add_attachment(fh.read(),
                                   maintype=maintype,
                                   subtype=subtype,
                                   filename=os.path.basename(path))

    with SMTP(_EBI_EMAIL_SERVER, port=_EBI_EMAIL_PORT) as s:
        s.send_message(msg)
