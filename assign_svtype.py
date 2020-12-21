"""Assigns SVTYPE to DUP,DEL,INVs in Svaba output vcf. 

All breakends are assigned as BND in raw Svaba output. This script can assign corresponding SVTYPE to the adjacencies that has strand information of DUP, DEL, INV, but not INS.
Some output BNDs will still contain INS. Note that for complex rearrangements, the SVTYPE assignment will be inaccurate. For more discussion on SVTYPE assignment, see https://github.com/walaj/svaba/issues/4.
"""

from pysam import VariantFile
import re
import sys
import alt_string_utils

vcf_in_path = sys.argv[1]
vcf_out_path = sys.argv[2]

vcf_in = VariantFile(vcf_in_path)
vcf_out = VariantFile(vcf_out_path, 'w', header = vcf_in.header)
vcf_out.header.add_meta(key = "INFO", items = [("ID", "SVLEN"), ("Number", "."), ("Type", "Integer"), ("Description", "Length of structural variant")])
vcf_out.header.add_meta(key="INFO", items=[("ID", "END"), ("Number", "."), ("Type", "Integer"), ("Description", "End coordinate")])

mates_dict = {}

def determine_svtype(list_mate_recs):
    """Determines SV type of a mate pair based on strands.

    Args:
        list_mate_recs (list): list of vcf records that are mates

    Returns:
        VariantRecord: Breakend record with new SVTYPE
    """
    #get alt, chrom1, chrom2, position (pos), id, old SVTYPE (should be BND if virgin svaba vcf) from line
    assert len(list_mate_recs) == 2, "Not a mate pair" + list_mate_recs[0]

    if list_mate_recs[0].id.split(":")[1] == "1":
        rec = list_mate_recs[0]
        mate = list_mate_recs[1]
    else:
        rec = list_mate_recs[1]
        mate = list_mate_recs[0] 

    oldType = rec.info["SVTYPE"]

    # get new type
    if oldType == 'BND' and rec.chrom == mate.chrom:
        assert rec.pos <= mate.pos, "Mate 1 has larger position than Mate 2 at: " + str(rec)
        
        rec_alt = rec.alts[0]
        mate_alt = mate.alts[0]

        # INV
        INV_PATTERN_THIS1 = re.compile(r'\D\].+\]')
        INV_PATTERN_MATE1 = re.compile(r'\D\].+\]')
        INV_PATTERN_THIS2 = re.compile(r'\[.+\[\D')
        INV_PATTERN_MATE2 = re.compile(r'\[.+\[\D')
        if INV_PATTERN_THIS1.match(rec_alt) and INV_PATTERN_MATE1.match(mate_alt) or INV_PATTERN_THIS2.match(rec_alt) and INV_PATTERN_MATE2.match(mate_alt):
                new_record = rec.copy()
                new_record.info["SVTYPE"] = "INV"
                new_record.info["SVLEN"] = str(mate.pos - rec.pos + 1)
                new_record.stop = mate.pos
                return [new_record]

        # DEL
        DEL_PATTERN_THIS = re.compile(r'\D\[.+\[')
        DEL_PATTERN_MATE = re.compile(r'\].+\]\D')
        if DEL_PATTERN_THIS.match(rec_alt) and DEL_PATTERN_MATE.match(mate_alt):
                new_record = rec.copy()
                new_record.info["SVTYPE"] = "DEL"
                new_record.info["SVLEN"] = str(mate.pos - rec.pos + 1)
                new_record.stop = mate.pos
                return [new_record]

        # DUP
        DUP_PATTERN_THIS = re.compile(r'\].+\]\D')
        DUP_PATTERN_MATE = re.compile(r'\D\[.+\[')
        if DUP_PATTERN_THIS.match(rec_alt)  and DUP_PATTERN_MATE.match(mate_alt):
                new_record = rec.copy()
                new_record.info["SVTYPE"] = "DUP"
                new_record.info["SVLEN"] = str(mate.pos - rec.pos + 1)
                new_record.stop = mate.pos
                return [new_record]
                
    return [rec, mate]

n_record = 0
for rec in vcf_in.fetch():
    n_record += 1
    record_id = rec.id
    record_mate_id = rec.info["MATEID"]
    mate_common_id = record_id.split(":")[0]

    assert mate_common_id == record_mate_id.split(":")[0]

    if mate_common_id in mates_dict.keys():
        mates_dict[mate_common_id].append(rec)
    else:
        mates_dict[mate_common_id] = [rec]

n_bnds = 0
n_dels = 0
n_dups = 0
n_invs = 0

for mate_common_id in mates_dict.keys():
    type_converted_BNDs = determine_svtype(mates_dict[mate_common_id])
    for converted_rec in type_converted_BNDs:
        if converted_rec.info["SVTYPE"] == "BND": n_bnds += 1
        elif converted_rec.info["SVTYPE"] == "DEL": n_dels += 1
        elif converted_rec.info["SVTYPE"] == "DUP": n_dups += 1
        elif converted_rec.info["SVTYPE"] == "INV": n_invs += 1
        vcf_out.write(converted_rec)

print("Converted %d BNDs into:"% n_record)
print("%d BNDs, %d DELs, %d DUPs, %d INVs"% (n_bnds, n_dels, n_dups, n_invs))

assert n_record == n_bnds + n_dels *2 + n_dups *2 + n_invs*2
