from typing import List

from .common import Method, parse_xml

_TYPE_PATTERNS = 'C'  # conserved sites
_TYPE_PROFILES = 'D'


def parse_patterns(filepath: str) -> List[Method]:
    """
    Parse the pattern_model.cgi file provided by PROSITE (XML format)
    http://prosite.expasy.org/cgi-bin/prosite/pattern_model.cgi

    :param filepath:
    :return:
    """
    return parse_xml(filepath, _TYPE_PATTERNS)


def parse_profiles(filepath: str) -> List[Method]:
    """
    Parse the profile_model.cgi file provided by PROSITE (XML format)
    http://prosite.expasy.org/cgi-bin/prosite/profile_model.cgi

    :param filepath:
    :return:
    """
    return parse_xml(filepath, _TYPE_PROFILES)
