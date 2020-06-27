# -*- coding: utf-8 -*-

import mimetypes
import os
from email.message import EmailMessage
from smtplib import SMTP
from typing import Optional, Sequence

_osos = Optional[Sequence[str]]  # optional sequence of strings


def send(info: dict, subject: str, content: str, attachments: _osos = None):
    if not info["To"]:
        return

    msg = EmailMessage()
    msg.set_content(content)
    msg["Sender"] = info["Sender"]
    msg["To"] = info["To"]

    if info["Cc"]:
        msg["Cc"] = ','.join(info["Cc"])

    if info["Bcc"]:
        msg["Bcc"] = ','.join(info["Bcc"])

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

    with SMTP(info["Server"], port=info["Port"]) as s:
        s.send_message(msg)
