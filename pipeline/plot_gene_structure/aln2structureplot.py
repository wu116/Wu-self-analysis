#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import logging
import argparse

def processGFF(GFFfile_path):
    gene_dict = {}
    trans_dict = {}
    exon_dict = {}
    CDS_dict = {}
    error = []
    with open(GFFfile_path, 'r') as GFFfile_path_file:
        for line in GFFfile_path_file:
            if line.startswith("#"):
                continue
            chr_id, _, feature, start, end, _, strand, _, attribute = line.strip().split('\t')
            start = int(start)
            end = int(end)
            if feature == "mRNA":
                RNA_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                mRNA_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                exon_count = CDS_count = 0
                if chr_id:
                    trans_dict[mRNA_id] = (RNA_parent_id, int(start), int(end), strand, chr_id)
                else:
                    error.append(line.strip())
            elif feature == "exon":
                exon_count += 1
                exon_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                try:
                    exon_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                except Exception as e:
                    exon_id = exon_count
                if strand not in ["+", "-"]:
                    error.append(line.strip())
                    continue
                if trans_dict[exon_parent_id]:
                    if not exon_parent_id in exon_dict:
                        exon_dict[exon_parent_id] = {}
                    exon_dict[exon_parent_id][str(exon_id)] = (int(start), int(end), strand)
                else:
                    error.append(line.strip())
            elif feature == "CDS":
                CDS_count += 1
                CDS_parent_id = re.search(r"Parent=([^\s,;]+)", attribute).group(1)
                try:
                    CDS_id = re.search(r"ID=([^\s,;]+)", attribute).group(1)
                except Exception as e:
                    CDS_id = CDS_count
                if strand not in ["+", "-"]:
                    error.append(line.strip())
                    continue
                if trans_dict[CDS_parent_id]:
                    if not CDS_parent_id in CDS_dict:
                        CDS_dict[CDS_parent_id] = {}
                    CDS_dict[CDS_parent_id][str(CDS_id)] = (int(start), int(end), strand)
                else:
                    error.append(line.strip())
    return trans_dict, exon_dict, CDS_dict, error

def processRpsbProc(RpsbProcout_path, domainDB):
    cdd_dict = {}
    with open(RpsbProcout_path, 'r') as RpsbProcout_path_file:
        session_ord = None
        for line in RpsbProcout_path_file:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("SESSION"):
                _, session_ord, _, _, _, _, _ = line.strip().split('\t')
            elif line.startswith("QUERY"):
                _, query_id, _, query_length, query_name = line.strip().split('\t')
                query_name = query_name.split(" ")[0]
                domain_id = 0
                continue
            elif line.startswith("DOMAINS"):
                domain_id += 1
                hsp_id = 0
                read_switch = 1
                continue
            elif line.startswith("ENDDOMAINS"):
                read_switch = 0
                continue
            elif not session_ord is None and read_switch == 1:
                if line.startswith(session_ord):
                    hsp_id += 1
                    sord, qid, hit_type, _, start, end, evalue, bitscore, accession, short_name, inC, _ = line.strip().split('\t')
                    if sord == session_ord and qid == query_id:
                        if not query_name in cdd_dict:
                            cdd_dict[query_name] = {}
                            cdd_dict[query_name][hsp_id] = (int(start), int(end), float(evalue), accession, short_name)
                        else:
                            cdd_dict_iteration = list(cdd_dict[query_name].items())
                            compare_result = True
                            old_id_list = []
                            for hsp, (hsp_start, hsp_end, hsp_evalue, hsp_accession, hsp_short_name) in cdd_dict_iteration:
                                compare = (int(end) >= int(hsp_start) + 15) and (int(start) <= int(hsp_end) - 15)
                                if compare:
                                    if hsp_short_name.upper() in short_name.upper() or short_name.upper() in hsp_short_name.upper():
                                        if domainDB in accession and domainDB in hsp_accession:
                                            if float(hsp_evalue) > float(evalue):
                                                old_id_list.append(hsp)
                                                continue
                                            else:
                                                compare_result = False
                                        elif domainDB in accession and not domainDB in hsp_accession:
                                            old_id_list.append(hsp)
                                            continue
                                        elif not domainDB in accession and domainDB in hsp_accession:
                                            compare_result = False
                                        else:
                                            if float(hsp_evalue) > float(evalue):
                                                old_id_list.append(hsp)
                                                continue
                                            else:
                                                compare_result = False
                                    else:
                                        continue
                                else:
                                    continue
                            if compare_result:
                                old_id_list_len = len(old_id_list)
                                if old_id_list_len == 1:
                                    hsp_replace = old_id_list[0]
                                    cdd_dict[query_name][hsp_replace] = (int(start), int(end), float(evalue), accession, short_name)
                                elif old_id_list_len > 1:
                                    for i in range(1, old_id_list_len + 1):
                                        old_id = old_id_list.pop()
                                        del cdd_dict[query_name][old_id]
                                    cdd_dict[query_name][old_id] = (int(start), int(end), float(evalue), accession, short_name)
                                else:
                                    cdd_dict[query_name][hsp_id] = (int(start), int(end), float(evalue), accession, short_name)

            else:
                continue
    return cdd_dict

def domainOrdChange(cdd_dict, trans_dict, CDS_dict):
    domain_dict = {}
    for query, hsp_data in cdd_dict.items():
        domain_dict[query] = {}
        for hsp, (hsp_start, hsp_end, hsp_evalue, hsp_accession, hsp_short_name) in hsp_data.items():
            hsp_length = (hsp_end - hsp_start + 1) * 3
            hsp_start = (hsp_start - 1) * 3 + 1
            hsp_end = hsp_end * 3
            (gene_id, trans_start, trans_end, trans_strand, chr_id) = trans_dict[query]
            CDS_list = []
            CDS_list_length = 0
            domain_list = []
            CDS_dict_sorted = sorted(CDS_dict[query].items(), key = lambda x:x[1][0])
            for CDS_id, (CDS_start, CDS_end, CDS_strand) in CDS_dict_sorted:
                CDS_list.append([CDS_start, CDS_end, CDS_strand])
                CDS_list_length += 1
            if trans_dict[query][3] == '+':
                for i in range(1, CDS_list_length + 1):
                    head_list = CDS_list.pop(0)
                    head_start, head_end, head_strand = head_list
                    if head_end - head_start + 1 < hsp_start:
                        hsp_start = hsp_start - (head_end - head_start + 1)
                        continue
                    if head_start + hsp_start - 1 + hsp_length - 1 > head_end:
                        hsp_length = (head_start + hsp_start - 1 + hsp_length - 1) - head_end + 1
                        domain_list.append([head_start + hsp_start - 1, head_end, head_strand, hsp_short_name])
                        hsp_start = 1
                    else:
                        domain_list.append([head_start + hsp_start - 1, head_start + hsp_start - 1 + hsp_length - 1, head_strand, hsp_short_name])
                        break
            elif trans_dict[query][3] == '-':
                for i in range(1, CDS_list_length + 1):
                    head_list = CDS_list.pop(-1)
                    head_start, head_end, head_strand = head_list
                    if head_end - head_start + 1 < hsp_start:
                        hsp_start = hsp_start - (head_end - head_start + 1)
                        continue
                    if head_end - hsp_start + 1 - hsp_length + 1 < head_start:
                        hsp_length = head_start - (head_end - hsp_start + 1 - hsp_length + 1)
                        domain_list.append([head_start, head_end - hsp_start + 1, head_strand, hsp_short_name])
                        hsp_start = 1
                    else:
                        domain_list.append([head_end - hsp_start + 1 - hsp_length + 1, head_end - hsp_start + 1, head_strand, hsp_short_name])
                        break

            domain_dict[query][hsp] = domain_list
    return domain_dict

def printoutput(trans_dict, exon_dict, CDS_dict, domain_dict):
    gff_lines = []
    trans_dict_sorted = sorted(trans_dict.items(), key = lambda x: (int(re.search(r"(\d+)", x[1][4]).group(1)), x[1][1]))
    for mRNA_id, (gene_id, mRNA_start, mRNA_end, mRNA_strand, chr_id) in trans_dict_sorted:
        gff_lines.append('\t'.join([chr_id, "CDD", "mRNA", str(mRNA_start), str(mRNA_end), ".", mRNA_strand, ".", mRNA_id, gene_id, "."]))

        exon_dict_sorted = sorted(exon_dict[mRNA_id].items(), key = lambda x:x[1][0])
        for exon_id, (exon_start, exon_end, exon_strand) in exon_dict_sorted:
            gff_lines.append('\t'.join([chr_id, "CDD", "exon", str(exon_start), str(exon_end), ".", exon_strand, ".", mRNA_id, gene_id, "."]))

        CDS_dict_sorted = sorted(CDS_dict[mRNA_id].items(), key=lambda x: x[1][0])
        for CDS_id, (CDS_start, CDS_end, CDS_strand) in CDS_dict_sorted:
            gff_lines.append('\t'.join([chr_id, "CDD", "CDS", str(CDS_start), str(CDS_end), ".", CDS_strand, ".", mRNA_id, gene_id, "."]))

        if mRNA_id in domain_dict:
            domain_dict_sorted = sorted(domain_dict[mRNA_id].items(), key=lambda x: x[1][0])
            for domain_id, domain_line in domain_dict_sorted:
                if mRNA_strand == "-":
                    domain_line = list(reversed(domain_line))
                for domain_start, domain_end, domain_strand, short_name in domain_line:
                    gff_lines.append('\t'.join([chr_id, "CDD", "domain", str(domain_start), str(domain_end), ".", domain_strand, ".", mRNA_id, gene_id, short_name]))
    return gff_lines

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="aln2structureplot")

    parser.add_argument('-g', '--gff3', required=True, help='the gff3 annotation of all gene in .aln')
    parser.add_argument('-c', '--cdd', required=True, help='CDD annotation from rpsblast and rpsbProc')
    parser.add_argument('-d', '--db', required=True, default = 'pfam', choices = ['CHL','COG','KOG','MTH',
                                                                                  'NF','PHA','PLN','PRK','PTZ','TIGR',
                                                                                  'cd','pfam','sd','smart'], help='use annotation from which database as priority')
    parser.add_argument('-o', '--out', required=True, help='ggtranscript plot gff file')
    args = parser.parse_args()

    gff3 = args.gff3
    cdd = args.cdd
    out = args.out
    domainDB = args.db

    trans_dict, exon_dict, CDS_dict, error = processGFF(gff3)
    cdd_dict = processRpsbProc(cdd, domainDB)
    domain_dict = domainOrdChange(cdd_dict, trans_dict, CDS_dict)
    gff_lines = printoutput(trans_dict, exon_dict, CDS_dict, domain_dict)

    with open(out, 'w') as f:
        for line in gff_lines:
            f.write(line + '\n')

    with open(out + "_errorlines.txt", 'w') as e:
        for line in error:
            e.write(line + '\n')

if __name__ == "__main__":
    main()