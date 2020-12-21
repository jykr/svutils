"""Utilities to extract information from ALT string or pysam.VariantRecord"""

import re
from pysam import VariantRecord

def strands(alt_string):
    """Extract strand orientation from ALT column of vcf file.

    Args:
        alt_string: ALT column of vcf file

    Returns:
        str: length 2 string with strand information. Returns None if alt_string doesn't match VCF spec of SVs.
    """
    pp_pattern = re.compile(r'\D\].+\]')
    pm_pattern = re.compile(r'\D\[.+\[')
    mp_pattern = re.compile(r'\].+\]\D')
    mm_pattern = re.compile(r'\[.+\[\D')

    if pp_pattern.match(alt_string): return "++"
    elif pm_pattern.match(alt_string): return "+-"
    elif mp_pattern.match(alt_string): return "-+"
    elif mm_pattern.match(alt_string): return "--"

    return None


def alt_pos_string(rec):
    """Extract mate breakend position string from ALT column.

    Args:
        rec: VariantRecord object of pysam

    Returns:
        str: string with mate breakend position (ex."1:2343253")
    """
    return re.search(r"[\[\]](.*)[\[\]]", rec.alts[0]).group(1)

