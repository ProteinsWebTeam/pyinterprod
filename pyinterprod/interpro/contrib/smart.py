from typing import List

from .common import Method, parse_xml

_TYPE = 'D'


def parse_signatures(db_sources: dict) -> List[Method]:
    """
    Parse the SIGNATURE_LIBRARY.xml provided by SMART.
    """
    return parse_xml(db_sources["sig_source"], _TYPE)
