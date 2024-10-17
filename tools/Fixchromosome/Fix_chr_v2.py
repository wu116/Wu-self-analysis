#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import sys
import argparse
import logging
import gffutils
import subprocess
#Author = Alter X
def Break(args):
    old_chr_fasta_file = args.chr
    broken_contig_file = args.out

    seq_id = None
    seq = ""
    seq_data = []

    with open(old_chr_fasta_file, "r") as fasta_input:
        for line in fasta_input:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line[1:]
            else:
                seq += line

    seq_data = [match.group() for match in re.finditer(r"[ACGT]+|N+", seq)]

    num = 1
    start = 0
    end = 0

    with open(broken_contig_file, "w") as broken_contig:
        for sub_seq in seq_data:
            start = end + 1
            length = len(sub_seq)
            end = start + length - 1
            if re.match(r"[ACGT]+", sub_seq):
                output = [seq_id, f"contig{num}", start, end, length, "W"]
                broken_contig.write("\t".join(map(str, output)) + "\n")
                num += 1
            elif re.match(r"N+", sub_seq):
                output = [seq_id, f"contig{num}", start, end, length, "U"]
                broken_contig.write("\t".join(map(str, output)) + "\n")
                num += 1
            else:
                sys.exit(f"Split seq failed")


def Reorder(args):
    new_contig_list_file = args.list
    gene_overlap_contigs_file = '%s.mRNA_overlap_contigs_file.txt' % args.prefix
    new_gff3_file = '%s.fixed.gff3' % args.prefix
    old_gff3_file = args.gff
    new_chr_fasta_file = '%s.fixed.fa' % args.prefix
    old_chr_fasta_file = args.chr
    gff3_db = gffutils.create_db(old_gff3_file, dbfn = ":memory:", force = True, keep_order = True)
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

    #seq = {}
    #current_seq_id = None

    #with open(old_chr_fasta_file, "r") as old_chr_fasta:
    #    for line in old_chr_fasta:
    #        line = line.strip()
    #        if line.startswith(">"):
    #            current_seq_id = line[1:]
    #            seq[current_seq_id] = ""
    #        else:
    #            seq[current_seq_id] += line

    for data_row in new_contig_data:
        chr_id = data_row[0]
        old_start = int(data_row[2])
        old_end = int(data_row[3])
        length = int(data_row[4])
        seqkit_command = "seqkit subseq -r {}:{} {}".format(old_start,old_end,old_chr_fasta_file)
        subseq = subprocess.run(seqkit_command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
        if subseq.returncode == 0:
            lines = subseq.stdout.split('\n')
            for line in lines:
                match = re.match(r'^>(.*)',line)
                if match and seq_id is None:
                    seq_id = match.group(1)
                    new_chr_fasta.write(">" + seq_id + "\n")
                    continue
                elif not match:
                    stripped_line = line.strip()
                    new_chr_fasta.write(stripped_line)
        else:
            print(subseq.stderr)
        #subseq = seq[chr_id][old_start-1:old_end]
        #if seq_id is None:
        #    seq_id = chr_id
        #    new_chr_fasta.write(">" + seq_id + "\n")
        #elif seq_id is not None:
        #    new_chr_fasta.write(subseq)

        overlap_element_db = gff3_db.region(seqid = chr_id, start = old_start, end = old_end, completely_within = False)
        for element in overlap_element_db:
            if element.start < old_start or element.end > old_end:
                new_gff3.write(str(element) + "\n")
                if element.featuretype == "gene":
                    if "ID" in element.attributes:
                        gene_id = element.attributes["ID"][0]
                        err.write(str(element) + "\n")
                        for child in gff3_db.children(gene_id, level = None):
                            err.write(str(child) + "\n")
			#最重要的一段，计算contig坐标的变换
			#每个元素的旧起始/终止位点减去contig的旧起始位点，就能得到每个元素相对于contig的偏移值
			#再加上contig的新起始位点，就能得到每个元素的新坐标
            else:
                element.start = element.start - old_start + new_start
                element.end = element.end - old_start + new_start
                new_gff3.write(str(element) + "\n")
		#基于contig的长度
        new_start += length
    new_gff3.close()
    err.close()
    new_chr_fasta.close()
    #subprocess.run("seqkit seq -w 60 {} > {}".format(new_chr_fasta_file, new_chr_fasta_file + ".temp"), shell=True)
    #subprocess.run("move /Y {} {}".format(new_chr_fasta_file + ".temp", new_chr_fasta_file), shell=True)
    #subprocess.run("sort /unique {} > {}".format(gene_overlap_contigs_file, gene_overlap_contigs_file + ".temp"), shell=True)
    #subprocess.run("move /Y {} {}".format(gene_overlap_contigs_file + ".temp", gene_overlap_contigs_file), shell=True)
    #subprocess.run("findstr /v /g:{} {} > {}".format(gene_overlap_contigs_file, new_gff3_file, new_gff3_file + ".temp"), shell=True)
    #subprocess.run("move /Y {} {}".format(new_gff3_file + ".temp", new_gff3_file), shell=True)

def wrap_fasta(args, width):
    fasta_file = '%s.fixed.fa' % args.prefix
    wrapped_seq = []
    lines = []
    seq_id = None

    with open(fasta_file, "r") as input_fasta:
        lines = input_fasta.readlines()

    for line in lines:
        stripped_line = line.strip()
        if stripped_line.startswith(">"):
            seq_id = stripped_line[1:]
        else:
            for i in range(0, len(stripped_line), width):
                wrapped_seq.append(stripped_line[i:i + width])

    with open(fasta_file, "w") as output_fasta:
        output_fasta.write(">" + seq_id + "\n")
        output_fasta.writelines(line + "\n" for line in wrapped_seq)

def sort_uniq_overlap_file(args):
    gene_overlap_contigs_file = '%s.mRNA_overlap_contigs_file.txt' % args.prefix
    lines = []

    with open(gene_overlap_contigs_file, "r") as input_gene_overlap_contigs:
        lines = input_gene_overlap_contigs.readlines()

    sort_uniq_lines = sorted(set(lines), key=lambda x: (int(x.split()[3]), sort_feature(x.split()[2]), int(x.split()[4])))

    with open(gene_overlap_contigs_file, "w") as output_gene_overlap_contigs:
        output_gene_overlap_contigs.writelines(sort_uniq_lines)

def sort_feature(element):
    order = {'gene' : 1, 'mRNA' : 2, 'exon' : 3 , 'CDS' : 4}
    type_str = element.strip()
    return order.get(type_str, 5)

def grep_vf(args):
    element_overlap_contigs_file = '%s.mRNA_overlap_contigs_file.txt' % args.prefix
    new_gff3_file = '%s.fixed.gff3' % args.prefix
    patterns = []
    lines = []

    with open(element_overlap_contigs_file, "r") as element_overlap_contigs:
        patterns = set(line.strip().split()[8] for line in element_overlap_contigs)

    with open(new_gff3_file, "r") as input_new_gff:
        lines = input_new_gff.readlines()

    with open(new_gff3_file, "w") as output_new_gff:
        for line in lines:
            strip_line = line.strip()
            if not any(pattern in strip_line for pattern in patterns):
                output_new_gff.write(strip_line + "\n")

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Fix gff3 file with reorder broken contigs, must run this script on a single chromosome once time")
    parser.add_argument('-m', '--mode', required=True, choices=['Break', 'Reorder'], help='Choose the mode,\n'
                                                                                          '"Break" for breaking the draft genome,\n'
                                                                                          '"Reorder" for reordering the draft genome with manul check contigs list\n')

    parser.add_argument('-g', '--gff', help='old gff3 file')
    parser.add_argument('-c', '--chr', help='old Chromosome fasta file')
    parser.add_argument('-o', '--out', help='out broken contigs list file')
    parser.add_argument('-l', '--list', help='a file with reorder broken contigs\n'
                                                            'chr\tb_contig\tstart\tend\tlength\tgap(W/U)\n'
                                                            'chr1\tb_contig1\t1\t5000\t5000\tW\n'
                                                            'chr1\tb_contig2\t5001\t5100\t100\tU\n')
    parser.add_argument('-p', '--prefix', default='result', help='prefix of output, final output contains: \n'
                                                                 '"prefix.element_overlap_contigs_file.txt" - a file show old genes with overlap on different contigs, especially overlap on gap\n'
                                                                 '"prefix.fixed.fa" - a file which fixed the Chromosome fasta file with reorder broken contigs\n'
                                                                 '"prefix.fixed.gff3" - a file which fixed the coordinate with reorder broken contigs')
    args = parser.parse_args()
    if args.mode == 'Break':
        if all([args.chr, args.out]):
            Break(args)
        else:
            parser.error("Break mode requires -c and -o arguments.")
    elif args.mode == 'Reorder':
        if all([args.gff, args.chr, args.list, args.prefix]):
            Reorder(args)
            wrap_fasta(args, 60)
            sort_uniq_overlap_file(args)
            grep_vf(args)
        else:
            parser.error("Reorder mode requires -g, -c, -l, and -p arguments.")

if __name__ == "__main__":
    main()