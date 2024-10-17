#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

def gff3_to_gtf(gff3_file, gtf_file):
    with open(gff3_file, 'r') as f:
        gff3_lines = f.readlines()

    gtf_lines = []
    for line in gff3_lines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) != 9:
            continue
        if 'RNA' in fields[2] or 'transcript' in fields[2]:
            attributes = fields[8].split(';')
            gene_id = ''
            transcript_id = ''
            for attr in attributes:
                if attr.startswith('Parent='):
                    gene_id = attr.replace('Parent=', '')
                elif attr.startswith('ID='):
                    transcript_id = attr.replace('ID=', '')
            if gene_id and transcript_id:
                gtf_line = '\t'.join([fields[0], fields[1], 'transcript', fields[3], fields[4], fields[5], fields[6], fields[7],
                                      'transcript_id "{}"; gene_id "{}";'.format(transcript_id, gene_id)]) + '\n'
                gtf_lines.append(gtf_line)
        elif fields[2] == 'exon':
            attributes = fields[8].split(';')
            exon_parent = ''
            for attr in attributes:
                if attr.startswith('Parent='):
                    exon_parent = attr.replace('Parent=', '')
            if exon_parent == transcript_id:
                gtf_line = '\t'.join([fields[0], fields[1], 'exon', fields[3], fields[4], fields[5], fields[6], fields[7],
                                      'transcript_id "{}"; gene_id "{}";'.format(exon_parent, gene_id)]) + '\n'
                gtf_lines.append(gtf_line)
        elif fields[2] == 'CDS':
            attributes = fields[8].split(';')
            CDS_parent = ''
            for attr in attributes:
                if attr.startswith('Parent='):
                    CDS_parent = attr.replace('Parent=', '')
            if CDS_parent == transcript_id:
                gtf_line = '\t'.join(
                    [fields[0], fields[1], 'CDS', fields[3], fields[4], fields[5], fields[6], fields[7],
                     'transcript_id "{}"; gene_id "{}";'.format(CDS_parent, gene_id)]) + '\n'
                gtf_lines.append(gtf_line)

    with open(gtf_file, 'w') as f:
        f.writelines(gtf_lines)


def gtf_to_gff3(gtf_file, gff3_file):
    with open(gtf_file, 'r') as f:
        gtf_lines = f.readlines()

    gff3_lines = []
    old_gene_id = ''
    for line in gtf_lines:
        fields = line.strip().split('\t')
        if len(fields) != 9:
            continue
        if fields[2] == 'transcript':
            attributes = fields[8].split(';')
            gene_id = ''
            transcript_id = ''
            for attr in attributes:
                if attr.strip().startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                elif attr.strip().startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
            if gene_id and transcript_id:
                if gene_id != old_gene_id:
                    gff3_line = '\t'.join([fields[0], fields[1], "gene", fields[3], fields[4], fields[5], fields[6], fields[7],
                         'ID={}'.format(gene_id)]) + '\n'
                    gff3_lines.append(gff3_line)
                gff3_line = '\t'.join([fields[0], fields[1], "transcript", fields[3], fields[4], fields[5], fields[6], fields[7],
                                       'ID={};Parent={}'.format(transcript_id, gene_id)]) + '\n'
                gff3_lines.append(gff3_line)
                old_gene_id = gene_id
        elif fields[2] == 'exon':
            attributes = fields[8].split(';')
            exon_parent = ''
            for attr in attributes:
                if attr.strip().startswith('transcript_id'):
                    exon_parent = attr.split('"')[1]
            if exon_parent == transcript_id:
                gff3_line = '\t'.join([fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7],
                                       'Parent={}'.format(exon_parent)]) + '\n'
                gff3_lines.append(gff3_line)
        elif fields[2] == 'CDS':
            attributes = fields[8].split(';')
            CDS_parent = ''
            for attr in attributes:
                if attr.strip().startswith('transcript_id'):
                    CDS_parent = attr.split('"')[1]
            if CDS_parent == transcript_id:
                gff3_line = '\t'.join([fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7],
                                       'Parent={}'.format(CDS_parent)]) + '\n'
                gff3_lines.append(gff3_line)

    with open(gff3_file, 'w') as f:
        f.writelines(gff3_lines)

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="GXF2GXF")

    parser.add_argument('-i', '--input', required=True, help='gtf file')
    parser.add_argument('-o', '--out', required=True, help='gff3 file')
    parser.add_argument('-m', '--mode', required=True, choices=['gff2gtf', 'gtf2gff'], help='choose change mode')
    args = parser.parse_args()

    if args.mode == 'gff2gtf':
        gff3_to_gtf(args.input, args.out)
    elif args.mode == 'gtf2gff':
        gtf_to_gff3(args.input, args.out)

if __name__ == "__main__":
    main()
