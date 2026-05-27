# Cross-Genus Analysis Reveals Architecturally Programmed sgRNA Synthesis Patterns in Coronaviruses

## Code packaging
The main analysis functions and plotting functions in the article have been packaged and released on PyPI. The installation command is: 
```python
pip install corona-rna==0.1.0
```
The processing pipeline Corona-mRNA for the third-generation sequencing data has been packaged as the nextflow framework.
usage：
```
nextflow run main.nf \
    --refseq ./data/reference/NC_045512.2.fasta \
    --fqdata ./data/fastq/SRR13089344.fastq \
    --v_chr NC_045512.2 \
    --title sample1 \
    --outdir ./results
```

## Project Description
Here, we use Jupyter to display how these functions are used and to help reproduce the main findings of the paper.

## Data available
Data can be  available on CODE (https://webofgroup.cn/dgwei/hngene).

## Related article
RNA Architecture Underlies Discontinuous Transcription and Evolution of Coronavirus

Zi Wen#, Lei Chen#, Dehua Luo#, Ju Sun#, Liangrong Guo, Yingxiang Deng, Zhiyuan Huang, Yuxiang Wang, Ke Pan, Fan Wang*, Shaobo Xiao*, Li Li* and Dengguo Wei*

Correspondence: fwang@mail.ccnu.edu.cn; vet@mail.hzau.edu.cn; li.li@mail.hzau.edu.cn; dgwei@mail.hzau.edu.cn

Coronaviruses employ discontinuous transcription to produce canonical subgenomic RNAs (sgRNAs) essential for gene expression. Although TRS-dependent template switching mechanism is proposed, its structural basis remains poorly defined, and the functional significance of abundant non-canonical sgRNAs persists as a critical gap since the discovery of discontinuous RNA synthesis. Here, we bridge this gap through the first cross-genus integrated analysis of coronavirus transcriptomes and RNA interactomes. We reveal that canonical sgRNA formation is associated with same-direction RNA-RNA interactions. In contrast, non-canonical sgRNAs form through distinct architectural mechanisms: short-range junctions mediated by stem-loop structures overlapping genomic deletion hotspots, and conserved long-range ORF1a-N interactions generating sgRNAs encoding immune-modulatory ORFs — a function never previously attributed to non-canonical transcription. These findings establish architecturally programmed discontinuous RNA synthesis and highlight a potential link between non-canonical sgRNAs, genomic plasticity, and immune modulation, which may have implications for coronavirus adaptation.
