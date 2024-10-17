#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import logging
import argparse
import pandas as pd
import numpy as np
import concurrent.futures

def read_gene_gff(file_path):
    rows = []
    with open(file_path, 'r') as file_path_file:
        for line in file_path_file:
            if line.startswith("#"):
                continue
            chr_id, _, feature, start, end, _, strand, _, attribute = line.strip().split()
            start = int(start)
            end = int(end)
            if feature != "gene":
                continue
            gene_id = re.search(r"ID=(\S+)", attribute).group(1)
            rows.append((chr_id, start, end, strand, gene_id))
        gene_df = pd.DataFrame(rows, columns = ['chr_id', 'start', 'end', 'strand', 'gene_id'])
        gene_df['start'] = pd.to_numeric(gene_df['start'], errors='coerce')
        gene_df['end'] = pd.to_numeric(gene_df['end'], errors='coerce')
    return gene_df

def read_meth_chunk(file_path, chunksize):
    chunk = pd.read_csv(file_path,
                        sep = '\t',
                        names = ['chr_id', 'start', 'strand', 'feature', 'meth_reads', 'total_reads'],
                        comment = '#',
                        chunksize = chunksize)
    return chunk

def process_chunk(gene_df, chunk, l, r):
    #print("processing",chunk.iloc[0,0:1])
    gene_chunk_df = pd.DataFrame()
    for idx, row in gene_df.iterrows():
        if row['strand'] == '+':
            condition = ((chunk['chr_id'] == row['chr_id']) &
                         (chunk['start'] >= row['start'] - l) &
                         (chunk['start'] <= row['end'] + r))
        elif row['strand'] == '-':
            condition = ((chunk['chr_id'] == row['chr_id']) &
                         (chunk['start'] >= row['start'] - r) &
                         (chunk['start'] <= row['end'] + l))
        filter_df = chunk[condition].copy()
        if not filter_df.empty:
            if row['strand'] == '+':
                condition_up = (filter_df['start'] < row['start'])
                condition_body = ((filter_df['start'] >= row['start']) & (filter_df['start'] <= row['end']))
                condition_down = (filter_df['start'] > row['end'])
                filter_df = filter_df.assign(distance=np.select(
                    [condition_up, condition_body, condition_down],
                    [filter_df['start'] - (row['start'] - l) + 1, 0, filter_df['start'] - row['end']],
                    default='Unknown'))
            elif row['strand'] == '-':
                condition_up = (filter_df['start'] > row['end'])
                condition_body = ((filter_df['start'] >= row['start']) & (filter_df['start'] <= row['end']))
                condition_down = (filter_df['start'] < row['start'])
                filter_df = filter_df.assign(distance=np.select(
                    [condition_up, condition_body, condition_down],
                    [(row['end'] + l) - filter_df['start'] + 1, 0, row['start'] - filter_df['start']],
                    default='Unknown'))
            filter_df = filter_df.assign(meth_radio = filter_df['meth_reads'] / filter_df['total_reads'],
                                         gene_id = row['gene_id'])
            filter_df = filter_df.assign(postiton = np.select(
                                             [condition_up, condition_body, condition_down],
                                             ['up', 'body', 'down'],
                                             default = 'Unknown'))
            gene_chunk_df = pd.concat([gene_chunk_df, filter_df], ignore_index=True)
    return gene_chunk_df

def catch_meth_near_gene(gene_df, meth_file, chunksize, l, r, cpu):
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        futures = []
        chunks = read_meth_chunk(meth_file, chunksize)
        for chunk in chunks:
            futures.append(executor.submit(process_chunk, gene_df, chunk, l, r))
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                if not result.empty:
                    results.append(result)
            except Exception as e:
                print(f"Exception while processing future: {e}")
    combined_df = pd.concat(results, ignore_index=True)
    return combined_df

def sort_df(df):
    print("Start sort result")
    def extract_number(text):
        match = re.search(r'\d+', text)
        return int(match.group()) if match else -float('inf')
    df['first_col_number'] = df['chr_id'].apply(extract_number)
    df_sorted = df.sort_values(by=['first_col_number', 'start'])
    df_sorted = df_sorted.drop(columns=['first_col_number'])
    print("Finish sort result")
    return df_sorted

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Catch Meth base nearby genes")

    parser.add_argument('-G', '--gene', required=False, help='gene gff3 file')
    parser.add_argument('-M', '--meth', required=False, help='meth file')
    parser.add_argument('-l', '--left', type=int, required=False, help='left length bp')
    parser.add_argument('-r', '--right', type=int, required=False, help='right length bp')
    parser.add_argument('-o', '--out', required=False, help='output file')
    parser.add_argument('-c', '--chunksize', type=int, default=200000, help='Number of rows per chunk for meth file')
    parser.add_argument('-t', '--thread', type=int, default=8, help='max CPU for using')
    args = parser.parse_args()

    gene_df = read_gene_gff(args.gene)
    gene_meth_df = catch_meth_near_gene(gene_df, args.meth, args.chunksize, args.left, args.right, args.thread)
    gene_meth_sorted_df = sort_df(gene_meth_df)
    gene_meth_sorted_df.to_csv(args.out, sep='\t', index=False)

if __name__ == "__main__":
    main()