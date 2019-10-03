import logging


__version__ = "1.0.3"


logger = logging.getLogger(__name__)

if not logger.hasHandlers():
    logger.setLevel(logging.INFO)
    _ch = logging.StreamHandler()
    _ch.setFormatter(
        logging.Formatter(
            fmt='%(asctime)s: %(levelname)s: %(message)s',
            datefmt='%y-%m-%d %H:%M:%S'
        )
    )
    logger.addHandler(_ch)
