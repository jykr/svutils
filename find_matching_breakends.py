"""Compare how many of svtk standardized MANTA vcf output has match in SVTYPE-assigned svaba output vcf.

Search if target vcf file has breakend that corresponds to the each breakend entry in query vcf. Note that a breakend here is a single row of vcf file. A pair of breakends will correspond to a SV for most cases.
Input vcfs for this script should both have SVTYPE assignment. This script ignores which sequence is inserted for INS variants. 

Query breakend has its corresponding target breakend if for a query breakpoint,
    1. Target breakends' position is within dist (default 100bp) and
    2. Breakend orientation ("strands") match.
    3. ALT position matches (see VCF file specification > SV annotation for more info)

Example:
    $ python find_matching_breakends.py query.vcf target.vcf output.txt
"""

from pysam import VariantFile
import re
import sys
import alt_string_utils

vcf_query = VariantFile(sys.argv[1])
vcf_target = VariantFile(sys.argv[2])


dist = 100


def record_matches(query, target, ignore_strands = False, ignore_alt_pos = False):
    """Returns if a query VariantRecord object matches with the target VariantRecord object.
    It is assumed that query and target breakend position is closer than 100bp. For query and target breakend to be matched, their strand should match and mate breakend position should be less than dist (default 100bp).

    Args:
        query (pysam.VariantRecord): query breakend
        target (pysam.VariantRecord): target breakend
        ignore_strands (bool): If True, query and target matches even when their strand orientation isn't the same.
        ignore_alt_pos (bool): If True, mate breakend position don't need to be matched for query and target to be matched.

    Returns:
        bool: True if query and target matches, False otherwise.
    """
    ## Strand pattern matches
    assert len(query.alts) == 1 & len(target.alts) == 1
    query_strand = strands(query.alts[0])
    target_strand = strands(query.alts[0])
    if not ignore_strands and query_strand != target_strand: return False

    ## ALT position matches

    if not ignore_alt_pos:
        try:
            query_alt_chr, query_alt_pos = alt_pos_string(query).split(":")
            target_alt_chr, target_alt_pos = alt_pos_string(target).split(":")
        except ValueError as e:
            if target.info["SVTYPE"] == "INS" and query.info["SVTYPE"] != "INS" : return False
            elif target.info["SVTYPE"] == "INS" and query.info["SVTYPE"] == "INS": return True
            elif target.info["SVTYPE"] != "INS": 
                print("Non-insertion have wrong ALT grammar.", e)
                exit(1)

        if query_alt_chr != target_alt_chr: return False
        if abs( int(query_alt_pos) - int(target_alt_pos)) > dist: return False
    return True


def has_match(query, targets, ignore_strands = False, ignore_alt_pos = False):
    """Returns if the query has any match within target.
    Calls record_matches for each target VariantRecord object that is within dist of the query breakend position.

    Args: 
        query (pysam.VariantRecord): query breakend
        targets (pysam.VariantFile): vcf file with target breakends
        ignore_strands, ignore_alt_pos (bool): argument to pass to record_matches

    Returns:
        bool: If there is any match, returns True. If there is no match, returns False.
    """
    within_distance_targets = targets.fetch(query.chrom, query.start - dist, query.start + dist)
    for candidate_hit in within_distance_targets:
        if record_matches(query, candidate_hit, ignore_strands = ignore_strands, ignore_alt_pos = ignore_alt_pos): return True

    return False




## Read vcf records from input file
out_file = open(sys.argv[3], "w")
for rec in vcf_query.fetch():
    if rec.info["SVTYPE"] == "BND" or rec.info["SVTYPE"] == "DEL" or rec.info["SVTYPE"] == "DUP":
        rec_has_match = has_match(rec, vcf_target)
    elif rec.info["SVTYPE"] == "INS":
        rec_has_match = has_match(rec, vcf_target, ignore_alt_pos = True, ignore_strands = True)
    out_file.write(rec.id + "\t" + rec.info["SVTYPE"] + "\t" + str(rec_has_match) + "\n")

out_file.close()


