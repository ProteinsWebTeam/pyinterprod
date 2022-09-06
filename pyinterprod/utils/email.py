import mimetypes
import os
from email.message import EmailMessage
from smtplib import SMTP
from typing import Sequence


def send(info: dict, to: Sequence[str], subject: str, content: str, **kwargs):
    cc = kwargs.get("cc", [])
    bcc = kwargs.get("bcc", [])
    attachments = kwargs.get("attachments", [])

    if not info:
        return

    try:
        host, port = info["server"].split(':')
    except ValueError:
        host = info["server"]
        port = 0
    else:
        port = int(port)

    msg = EmailMessage()
    msg.set_content(content)
    msg["From"] = info["sender"]
    msg["To"] = ", ".join({info[key] for key in to})

    if cc:
        msg["Cc"] = ", ".join({info[key] for key in cc})

    if bcc:
        msg["Bcc"] = ", ".join({info[key] for key in bcc})

    msg["Subject"] = subject
    msg.set_content(content)

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

    with SMTP(host=host, port=port) as s:
        s.send_message(msg)
