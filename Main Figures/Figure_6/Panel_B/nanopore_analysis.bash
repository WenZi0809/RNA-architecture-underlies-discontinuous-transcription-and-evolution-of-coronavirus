# Single-junction
python sam_to_bed_linux.py NC_045512.2 SARS2-PRJNA675370-vero-24h-SRR13089344_mapped.sam SARS2-PRJNA675370-vero-24h-SRR13089344_dyx.bed
python counts_gene_linux_no_plot.py NC_045512.2 SARS2-PRJNA675370-vero-24h-SRR13089344_dyx.bed SARS2-PRJNA675370-vero-24h-SRR13089344_dyx
python select_loc_from_nanopore_bed.py SARS2-PRJNA675370-vero-24h-SRR13089344_dyx-2.bed SARS2-PRJNA675370-vero-24h-SRR13089344_mapped.sam SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.sam
samtools view -@ 10 -b -S -o SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.bam SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.sam
samtools sort -@ 10 -o SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.sort.bam SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.bam
bedtools genomecov -ibam SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.sort.bam -dz -split > SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton_coverage.txt
python nanopore_bam_info_select_cigar_start.py SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.bam SARS2-PRJNA675370-vero-24h-SRR13089344_orf1aton.info

# Double-junction
grep "N-UTR_3" SARS2-PRJNA675370-vero-24h-SRR13089344_dyx-3.bed > SARS2-PRJNA675370-vero-24h-SRR13089344_2more_N-UTR_3.info
grep "leader-leader" SARS2-PRJNA675370-vero-24h-SRR13089344_2more_N-UTR_3.info >SARS2-PRJNA675370-vero-24h-SRR13089344_2more_N-UTR_3_leader.info
python get_id_from_nanopore_bed.py SARS2-PRJNA675370-vero-24h-SRR13089344_2more_N-UTR_3_leader.info SARS2-PRJNA675370-vero-24h-SRR13089344_2more_id.txt

python get_nanopore_sam_from_id.py SARS2-PRJNA675370-vero-24h-SRR13089344_mapped.sam SARS2-PRJNA675370-vero-24h-SRR13089344_2more_id.txt  SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.sam
samtools view -@ 10 -b -S -o SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.bam SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.sam
samtools sort -@ 10 -o SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.sort.bam SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.bam
bedtools genomecov -ibam SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.sort.bam -dz -split > SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10_coverage.txt
python nanopore_bam_info_select_cigar_start.py SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.sort.bam SARS2-PRJNA675370-vero-24h-SRR13089344_2more_orf10.info

samtools view -@ 10 -b -S -o SARS2-PRJNA675370-vero-24h-SRR13089344_mapped_dyx.bam SARS2-PRJNA675370-vero-24h-SRR13089344_mapped.sam
samtools sort -@ 10 -o SARS2-PRJNA675370-vero-24h-SRR13089344_mapped_dyx_sort.bam SARS2-PRJNA675370-vero-24h-SRR13089344_mapped_dyx.bam
bedtools genomecov -ibam SARS2-PRJNA675370-vero-24h-SRR13089344_mapped_dyx_sort.bam -dz -split > SARS2-PRJNA675370-vero-24h-SRR13089344_mapped_dyx_coverage.txt