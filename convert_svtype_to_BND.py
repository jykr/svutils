""" Converts all SVTYPE into BNDs.

This script can be used to compare SV caller output without assigning SVTYPE to callers with no SVTYPE annotation (ex. Svaba).
"""

from pysam import VariantFile
import re
import sys
import alt_string_utils

input_vcf_path = sys.argv[1]
output_vcf_path = sys.argv[2]
svlen_string = "SVLEN" # "SPAN" for Svaba

def convert_indel(rec):
    """Assign INFO/SVTYPE and INFO/SVLEN to INS/DEL that are of indel encoding (VCF spec). """
    if "SVTYPE" in rec.info.keys(): return rec
    alt = rec.alts[0]
    ref = rec.ref
    if len(alt) > len(ref):
        rec.info["SVTYPE"] = "INS"
        rec.info["SVLEN"] = rec.info["SPAN"]
        return(rec)
    else:
        rec.info["SVTYPE"] = "DEL"
        rec.info["SVLEN"] = -rec.info["SPAN"]
        return(rec)


def get_length(rec, source):
    """Get SVLEN from MANTA or Svaba vcf breakend record"""
    try:
        return abs(rec.info[svlen_string][0])
    except KeyError as e:
        if rec.info[svlen_string] == "INS":
            return len(rec.info["CONTIG"][0])
        print(e)


def convert_to_bnd(rec, source = "manta"):
    """Convert all breakends into BNDs."""

    if rec.info["SVTYPE"] == "INS":
        mate1 = rec.copy()
        mate2 = rec.copy()
        inserted_contig = rec.alts[0][len(rec.ref):]

        inserted_length = get_length(rec, source)
        mate2.pos = mate1.pos + 1
        
        # In case of insertions, we don't have information of contig to compare (for now)."
        mate1.alts = ("[ctg]",)
        mate2.alts = ("[ctg]",)

        return [mate1, mate2]

    elif rec.info["SVTYPE"] == "DEL":
        mate1 = rec.copy()
        mate2 = rec.copy()
        deleted_contig = rec.ref[len(rec.alts[0]):]
        
        deleted_length = get_length(rec, source)
        try: 
            mate2.pos = mate1.pos + deleted_length + 1
        except TypeError:
            print(mate2.pos, mate1.pos, deleted_length)
        
        mate1.ref = mate1.alts[0]
        mate1.alts = (mate1.ref + "[" + mate2.chrom + ":" + str(mate2.pos) + "[",)

        mate2.ref = "N"
        mate2.alts = ("]" + mate1.chrom + ":" + str(mate1.pos) + "]" + mate2.ref,)

        return [mate1, mate2]        

    elif rec.info["SVTYPE"] == "DUP":
        mate1 = rec.copy()
        mate2 = rec.copy()
        mate2.pos = mate1.pos + get_length(rec, source) + 1
        
        mate1.alts = ("]" + mate2.chrom + ":" + str(mate2.pos) + "]N", )
        mate2.alts = ("N[" + mate1.chrom + ":" + str(mate1.pos) + "[",)

        return [mate1, mate2]

    elif rec.info["SVTYPE"] == "INV":
        assert source == "manta", "No rule set for non-MANTA INV conversion"
        mate1 = rec.copy()
        mate2 = rec.copy()
        mate2.pos = mate1.stop

        if "INV5" in rec.info.keys():
            mate1.alts = ("[" + mate2.chrom + ":" + str(mate2.pos) + "[N",)
            mate2.alts = ("[" + mate1.chrom + ":" + str(mate1.pos) + "[N",)

        elif "INV3" in rec.info.keys():
            mate1.alts = ("N]" + mate2.chrom + ":" + str(mate2.pos) + "]",)
            mate2.alts = ("N]" + mate1.chrom + ":" + str(mate1.pos) + "]",)

        else:
            raise IOError("VCF file format does not have INV3/INV5 tag from MANTA output (INV converted)")

        return [mate1, mate2]

    return [rec]


def convert_vcffile(filename, outfile_name, source):

    vcf_in = VariantFile(filename)
    vcf_out = VariantFile(outfile_name, "w", header = vcf_in.header)

    for rec in vcf_in.fetch():
        if source == "svaba":
            rec = convert_indel(rec)
        for conv_rec in convert_to_bnd(rec, source):
            vcf_out.write(conv_rec)



convert_vcffile(input_vcf_path, output_vcf_path, source)
