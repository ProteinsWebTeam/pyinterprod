from typing import List

from .common import Method, parse_xml

_TYPE = 'D'


def parse_signatures(filepath: str) -> List[Method]:
    """
    Parse the SIGNATURE_LIBRARY.xml provided by SMART.
    """
    return parse_xml(filepath, _TYPE)
