# -*- coding: utf-8 -*-

import mimetypes
import os
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
    sender = getuser() + _EBI_EMAIL

    msg = EmailMessage()
    msg.set_content(content)
    msg["Sender"] = sender
    msg["To"] = set(to_addrs)

    if cc_addrs:
        msg["Cc"] = ','.join(set(cc_addrs))

    if bcc_addrs:
        bcc_addrs = set(bcc_addrs)
        bcc_addrs.add(sender)
    else:
        bcc_addrs = {sender}

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
                                   subtype=subtype,
                                   filename=os.path.basename(path))

    with SMTP(_EBI_EMAIL_SERVER, port=_EBI_EMAIL_PORT) as s:
        s.send_message(msg)
