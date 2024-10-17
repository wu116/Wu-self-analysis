#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <input.gff3> <output.list>";
die($usage) unless (@ARGV == 2);

# 输入文件和输出文件
my $input_file = $ARGV[0];
my $output_file = $ARGV[1];

# 打开输入文件
open(ERR,">",'err_line.txt');
open(my $in,'<',$input_file) or die "无法打开输入文件 $input_file: $!";
# 初始化变量
my %gene_data;
my $gene_id;

# 逐行读取GFF3文件
while (my $line = <$in>) {
    next if $line =~ /^#/;  # 跳过注释行
    chomp $line;
    my @fields = split(/\t/, $line);
    my ($feature, $start, $end, $strand, $attributes) = @fields[2, 3, 4, 6, 8];

    # 提取基因信息
    if ($feature eq 'gene') {
        $gene_id = $1 if($attributes =~ /ID=([^;]+)/);
        $gene_data{$gene_id} = {};
    }

    # 提取转录本信息
    if ($feature eq 'mRNA') {
        my $mRNA_id = $1 if($attributes =~ /ID=([^;]+)/);
        my $Parent_id = $1 if($attributes =~ /Parent=([^;]+)/);
        print ERR "$line\n" and next unless(defined $Parent_id);
        $gene_data{$Parent_id}{$mRNA_id} = {};
    }

    # 提取CDS信息
    if ($feature eq 'CDS') {
        my $CDS_id = $1 if($attributes =~ /ID=([^;]+)/);
        my $Parent_id = $1 if($attributes =~ /Parent=([^;]+)/);
        print ERR "$line\n" and next unless(defined $Parent_id);
        print ERR "$line\n" and next unless(exists $gene_data{$gene_id}{$Parent_id});
        # 存储CDS信息
        $gene_data{$gene_id}{$Parent_id}{$CDS_id} = {
            'start' => $start,
            'end'   => $end
        }
    }
}
# 关闭输入文件

close($in);

# 打开输出文件
open(my $out,'>',$output_file) or die "无法打开输出文件 $output_file: $!";
# 遍历基因数据，提取最长转录本
foreach my $gene (keys(%gene_data)){
    my $longest_mRNA;
    my $longest_mRNA_length = 0;
    foreach my $mRNA (keys(%{$gene_data{$gene}})){
        my $mRNA_length = 0;
        # 计算最长转录本长度
        while(my($CDS_key,$CDS_value) = each(%{$gene_data{$gene}{$mRNA}})){
            my $CDS_length = $CDS_value->{end} - $CDS_value->{start} + 1;
            $mRNA_length += $CDS_length;
        }
        if($longest_mRNA_length < $mRNA_length){
            $longest_mRNA_length = $mRNA_length;
            $longest_mRNA = $mRNA;
        }
    }
    print ERR "noCDS: $gene\n" and next unless(defined $longest_mRNA);
    print $out "$longest_mRNA\n";
}
# 关闭输出文件
close($out);
close(ERR);
print "提取完成，结果保存在 $output_file 中。\n";
