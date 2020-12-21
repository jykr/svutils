# svutils
### Descriptions
* **alt_string_utils.py**: Extracts information from ALT column of VCF file.
* **assign_svtype.py**: Assigns SVTYPE to VCF files with no SVTYPE assignment (ex. Svaba output) based on the record and its mate strand orientation.
* **convert_svtype_to_BND.py**: Converts all SVTYPEs, including INS, DELs in indel encoding into adjacencies with pair of breakends.
* **find_matching_breakends.py**: Finds if each breakend in query vcf file has its match in target vcf file.

