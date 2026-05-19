#####python nanopore_bam_info_select.py SARS2-PRJNA675370-vero-24h-SRR13089344_mapped.bam SARS2-PRJNA675370-vero-24h-SRR13089344_nanopore_info.txt
#####参数1：bam；参数2：输出文件

import pysam
import sys

pyname, bam_file, outfile = sys.argv

cigar = []
start = []

with pysam.AlignmentFile(bam_file) as alnfile:
    for aln in alnfile:
        cigar.append(aln.cigarstring)
        start.append(aln.reference_start)

with open(outfile,"w")as f:
    for i in range(len(cigar)):
        f.write(cigar[i]+"\t"+str(start[i])+"\n")