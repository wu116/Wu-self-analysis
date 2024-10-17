#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import logging
import argparse
import pandas as pd
import numpy as np


def read_gene_gff(file_path):
    rows = []
    with open(file_path, 'r') as file_path_file:
        for line in file_path_file:
            if line.startswith("#"):
                continue
            chr_id, _, feature, start, end, _, strand, _, attribute = line.strip().split()
            start = int(start) - 1
            end = int(end)
            if feature != "gene":
                continue
            gene_id = re.search(r"ID=(\S+)", attribute).group(1)
            rows.append((chr_id, start, end, strand, gene_id))
        gene_df = pd.DataFrame(rows, columns = ['chr_id', 'start', 'end', 'strand', 'gene_id'])
    return gene_df

def read_meth_tsv(file_path):
    meth_df = pd.read_csv(file_path,
                          sep = '\t',
                          names = ['chr_id', 'start', 'strand', 'feature', 'meth_reads', 'total_reads'],
                          comment = '#')
    meth_df['meth_radio'] = meth_df['meth_reads'] / meth_df['total_reads']
    return meth_df

def catch_meth_near_gene(gene_df, meth_df, l, r, o):
    gene_meth_df = pd.DataFrame(columns = meth_df.columns)
    for idx, row in gene_df.iterrows():
        if row['strand'] == '+':
            condition = ((meth_df['chr_id'] == row['chr_id']) &
                        (meth_df['start'] >= row['start'] - l) &
                        (meth_df['start'] <= row['end'] + r))
        elif row['strand'] == '-':
            condition = ((meth_df['chr_id'] == row['chr_id']) &
                        (meth_df['start'] >= row['start'] - r) &
                        (meth_df['start'] <= row['end'] + l))
        filter_df = meth_df[condition]
        filter_df['gene_id'] = row['gene_id']
        gene_meth_df = pd.concat([gene_meth_df, filter_df], ignore_index = True)
    gene_meth_df.to_csv(o, sep = '\t', index = False)

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Catch Meth base nearby genes")

    parser.add_argument('-G', '--gene', required=True, help='gene gff3 file')
    parser.add_argument('-M', '--meth', required=True, help='meth file')
    parser.add_argument('-l', '--left', type=int, required=True, help='left length bp')
    parser.add_argument('-r', '--right', type=int, required=True, help='right length bp')
    parser.add_argument('-o', '--out', required=True, help='output file')
    args = parser.parse_args()

    gene_df = read_gene_gff(args.gene)
    meth_df = read_meth_tsv(args.meth)
    catch_meth_near_gene(gene_df, meth_df, args.left, args.right, args.out)

if __name__ == "__main__":
    main()