# usage
```
nextflow run main.nf \
    --refseq ./data/reference/NC_045512.2.fasta \
    --fqdata ./data/fastq/SRR13089344.fastq \
    --v_chr NC_045512.2 \
    --title sample1 \
    --outdir ./results
```

# dependencies:
- minimap2
- samtools
- numpy
- matplotlib
