# -*- coding: utf-8 -*-

import mimetypes
from email.message import EmailMessage
from getpass import getuser
from smtplib import SMTP
from typing import Optional


_EBI_EMAIL_SERVER = "outgoing.ebi.ac.uk"
_EBI_EMAIL_PORT = 587
_EBI_EMAIL = "@ebi.ac.uk"


def send_mail(to_addrs: list, subject: str, content: str,
              cc_addrs: Optional[list]=None,
              bcc_addrs: Optional[list]=None,
              attachments: Optional[list]=None):
    # TODO: remove after debug
    to_addrs = ["mblum@ebi.ac.uk"]
    cc_addrs = None
    bcc_addrs = None
    
    msg = EmailMessage()
    msg.set_content(content)

    msg["Sender"] = getuser() + _EBI_EMAIL
    msg["To"] = set(to_addrs)
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
                msg.add_attachment(fh.read(), maintype=maintype,
                                   subtype=subtype, filename=path)

    with SMTP(_EBI_EMAIL_SERVER, port=_EBI_EMAIL_PORT) as s:
        s.send_message(msg)
