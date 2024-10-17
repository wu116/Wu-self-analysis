#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import sys
import argparse
import logging
import gffutils
import subprocess

def fix(args):
    new_contig_list_file = args.list
    gene_overlap_contigs_file = '%s.gene_overlap_contigs_file.gff3' % args.prefix
    new_gff3_file = '%s.fixed.gff3' % args.prefix
    old_gff3_file = args.gff
    new_chr_fasta_file = '%s.fixed.fa' % args.prefix
    old_chr_fasta_file = args.chr

    gff3_db = gffutils.create_db(old_gff3_file, dbfn = "gff3.db", force = True, keep_order = True)
    new_contig_data = []
    with open(new_contig_list_file, "r") as contig_file:
        for line in contig_file:
            columns = line.strip().split("\t")
            new_contig_data.append(columns)

    new_start = 1
    new_end = 1

    err = open(gene_overlap_contigs_file, "w")
    new_gff3 = open(new_gff3_file, "w")
    new_chr_fasta = open(new_chr_fasta_file, "w")

    seq_id = None

    for data_row in new_contig_data:
        chr_id = data_row[0]
        old_start = int(data_row[2])
        old_end = int(data_row[3])
        length = int(data_row[4])

        seqkit_command = "seqkit subseq -r {}:{} {}".format(old_start,old_end,old_chr_fasta_file)
        subseq = subprocess.run(seqkit_command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
        if subseq.returncode == 0:
            lines = subseq.split('\n')
            for line in lines:
                match = re.match(r'^>(.*)',line)
                if match and seq_id is None:
                    seq_id = match.group(1)
                    new_chr_fasta.write(">"+seq_id+"\n")
                    continue
                elif not match:
                    stripped_line = line.strip()
                    new_chr_fasta.write(stripped_line)
        else:
            print(subseq.stderr)

        overlap_element_db = gff3_db.region(seqid = chr_id, start = old_start, end = old_end, completely_within = False)

        for element in overlap_element_db:
            if element.start < old_start or element.end > old_end:
                err.write(str(element) + "\n")
                continue

            element.start = element.start - old_start + new_start
            element.end = element.end - old_start + new_start
            new_gff3.write(str(element) + "\n")

        new_start += length

    new_gff3.close()
    err.close()
    new_chr_fasta.close()

    subprocess.run("seqkit seq -w 60 {} > {}".format(new_chr_fasta_file, new_chr_fasta_file + ".temp"), shell=True)
    subprocess.run("mv {} {}".format(new_chr_fasta_file + ".temp", new_chr_fasta_file), shell=True)

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Fix gff3 file with reorder broken contigs, must run this script on a single chromosome once time")

    parser.add_argument('-g', '--gff', required=True, help='old gff3 file')
    parser.add_argument('-c', '--chr', required=True, help='old Chromosome fasta file')
    parser.add_argument('-l', '--list', required=True, help='a file with reorder broken contigs\n'
                                                            'chr\tb_contig\tstart\tend\tlength\tgap(W/U)\n'
                                                            'chr1\tb_contig1\t1\t5000\t5000\tW\n'
                                                            'chr1\tb_contig2\t5001\t5100\t100\tU\n')
    parser.add_argument('-p', '--prefix', default='result', help='prefix of output, final output contains: \n'
                                                                 '"prefix.gene_overlap_contigs_file.gff3" - a file show old genes with overlap on different contigs, especially overlap on gap\n'
                                                                 '"prefix.fixed.fa" - a file which fixed the Chromosome fasta file with reorder broken contigs\n'
                                                                 '"prefix.fixed.gff3" - a file which fixed the coordinate with reorder broken contigs')
    args = parser.parse_args()

    fix(args)

if __name__ == "__main__":
    main()
