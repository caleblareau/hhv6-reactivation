# Allogene collaboration
CL notes for Allogene collaboration

### Get the reference genome

```
https://www.ncbi.nlm.nih.gov/nuccore/AF157706
```

### Download scRNA data
```
wget https://cg.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_fastqs.tar
wget https://cg.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_fastqs.tar
```

### Quantify using kallisto | bustools

```
kallisto bus -t 8 -i HHV6b_transcriptome.idx -o count_v3 -x 10xv3 fastq/pbmc_10k_v3_S1_L002_R1_001.fastq.gz fastq/pbmc_10k_v3_S1_L002_R2_001.fastq.gz
kallisto bus -t 8 -i HHV6b_transcriptome.idx -o count_v5p -x 10xv2 fastq/vdj_v1_hs_pbmc2_5gex_S1_L002_R1_001.fastq.gz fastq/vdj_v1_hs_pbmc2_5gex_S1_L002_R2_001.fastq.gz
```
outcome:
```
[v3_3prime quant] processed 317,964,164 reads, 164 reads pseudoaligned
[5prime quant] processed 215,257,851 reads, 366 reads pseudoaligned
[allo_38low] processed 6,006,377 reads, 10 reads pseudoaligned
[allo_34high] processed 20,581,694 reads, 1,558 reads pseudoaligned
```

v3: 0.5 reads / million
5p: 1.7 reads / million


<br><br>
