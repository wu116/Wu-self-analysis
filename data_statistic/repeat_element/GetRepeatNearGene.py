#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import logging
import argparse


def read_gene_gff(file_path):
    gene_dict = {}
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
            if chr_id not in gene_dict:
                gene_dict[chr_id] = {}
            gene_dict[chr_id][gene_id] = (start, end, strand)
    return gene_dict


def read_repeat_gff(file_path):
    repeat_dict = {}
    with open(file_path, 'r') as file_path_file:
        for line in file_path_file:
            if line.startswith("#"):
                continue
            chr_id, _, feature, start, end, _, strand, _, attribute = line.strip().split()
            start = int(start)
            end = int(end)
            if re.match(r"Parent=", attribute):
                continue
            te_id = re.search(r"ID=(\S+?);", attribute).group(1)
            classification = re.search(r"Classification=(\S+?);", attribute).group(1)
            if chr_id not in repeat_dict:
                repeat_dict[chr_id] = {}
            repeat_dict[chr_id][te_id] = (start, end, classification)
    return repeat_dict


def intersect(gene_dict, repeat_dict, l, r, output_path):
    with open(output_path, 'w') as output_file:
        for chr_id, gene_info in gene_dict.items():
            if chr_id in repeat_dict:
                for gene_id, (gene_start, gene_end, gene_strand) in gene_info.items():
                    if gene_strand == "+":
                        up_length, down_length = l, r
                        for te_id, (te_start, te_end, te_class) in repeat_dict[chr_id].items():
                            orient, distance = None, None
                            if te_end > gene_start - up_length and te_start < gene_end + down_length:
                                if te_end <= gene_start:
                                    orient = "upstream"
                                    distance = gene_start - te_end
                                    output_list = [gene_id, te_id, te_class, orient, distance]
                                    output_file.write('\t'.join(map(str, output_list)) + '\n')
                                elif te_start >= gene_end:
                                    orient = "downstream"
                                    distance = te_start - gene_end
                                    output_list = [gene_id, te_id, te_class, orient, distance]
                                    output_file.write('\t'.join(map(str, output_list)) + '\n')
                                else:
                                    dis1 = te_end - gene_start
                                    dis2 = gene_end - te_start
                                    dis3 = gene_start - te_start
                                    dis4 = te_end - gene_end
                                    compare = [dis1, dis2, dis3, dis4]
                                    compare_min = min(map(abs, compare), default=None)
                                    if compare_min == dis1:
                                        orient = "upstream"
                                        distance = -dis1
                                    elif compare_min == dis2:
                                        orient = "downstream"
                                        distance = -dis2
                                    elif compare_min == abs(dis3):
                                        orient = "upstream"
                                        distance = dis3
                                    elif compare_min == abs(dis4):
                                        orient = "downstream"
                                        distance = dis4
                                    output_list = [gene_id, te_id, te_class, orient, distance]
                                    output_file.write('\t'.join(map(str, output_list)) + '\n')

                    elif gene_strand == "-":
                        up_length, down_length = r, l
                        for te_id, (te_start, te_end, te_class) in repeat_dict[chr_id].items():
                            if te_end > gene_start - up_length and te_start < gene_end + down_length:
                                if te_end <= gene_start:
                                    orient = "downstream"
                                    distance = gene_start - te_end
                                    output_list = [gene_id, te_id, te_class, orient, distance]
                                    output_file.write('\t'.join(map(str, output_list)) + '\n')
                                elif te_start >= gene_end:
                                    orient = "upstream"
                                    distance = te_start - gene_end
                                    output_list = [gene_id, te_id, te_class, orient, distance]
                                    output_file.write('\t'.join(map(str, output_list)) + '\n')
                                else:
                                    dis1 = te_end - gene_start
                                    dis2 = gene_end - te_start
                                    dis3 = gene_start - te_start
                                    dis4 = te_end - gene_end
                                    compare = [dis1, dis2, dis3, dis4]
                                    compare_min = min(map(abs, compare), default=None)
                                    if compare_min == dis1:
                                        orient = "downstream"
                                        distance = -dis1
                                    elif compare_min == dis2:
                                        orient = "upstream"
                                        distance = -dis2
                                    elif compare_min == abs(dis3):
                                        orient = "downstream"
                                        distance = dis3
                                    elif compare_min == abs(dis4):
                                        orient = "upstream"
                                        distance = dis4
                                    output_list = [gene_id, te_id, te_class, orient, distance]
                                    output_file.write('\t'.join(map(str, output_list)) + '\n')


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Catch TE elements nearby genes")

    parser.add_argument('-G', '--gene', required=True, help='gene gff3 file')
    parser.add_argument('-R', '--repeat', required=True, help='repeat gff3 file')
    parser.add_argument('-l', '--left', type=int, required=True, help='left length bp')
    parser.add_argument('-r', '--right', type=int, required=True, help='right length bp')
    parser.add_argument('-o', '--out', required=True, help='output file')
    args = parser.parse_args()

    gene_dict = read_gene_gff(args.gene)
    repeat_dict = read_repeat_gff(args.repeat)
    intersect(gene_dict, repeat_dict, args.left, args.right, args.out)


if __name__ == "__main__":
    main()
