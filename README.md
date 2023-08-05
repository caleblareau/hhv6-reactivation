# HHV6 Reactivation
Repository for reproducing analyses related to HHV-6 reactivation in (CAR) T cells manuscript:

Lareau et al. 2023 "Latent human herpesvirus 6 is reactivated in chimeric antigen receptor T cells." _Nature. In press._

## Reproducibility

### Raw sequencing data
Our .fastq is made freely available on the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210063).

### Figure 1
- [This repository](https://github.com/caleblareau/serratus-reactivation-screen) contains files to reproduce Serratus screening analyses for all viruses (via ViralZone) known to infect humans (see panels b,c). Note: due to the large number of `.soft` files annotating the GSM samples, it made sense to make this a separate repository.
- Code to reproduce the followup and supplemental analyses panels are contained in [this subfolder](https://github.com/caleblareau/hhv6-reactivation/tree/main/public_bulk_data).

### Figure 2
- All code to reproduce the analysis figure panels are contained in [this subfolder](https://github.com/caleblareau/hhv6-reactivation/tree/main/single_cell_data/stanford-invitro/code).

### Figure 3
- All code to reproduce the analysis figure panels are contained in [this subfolder](https://github.com/caleblareau/hhv6-reactivation/tree/main/single_cell_data/broad-invivo/code).

### Figure 4
- All code to reproduce the analysis figure panels are contained in [this subfolder](https://github.com/caleblareau/hhv6-reactivation/tree/main/single_cell_data/cd7-cart/code).



### HHV6B reference genome

For both ChIP-seq and RNA-seq followup analyses, we use the genome / transcriptome of HHV-6B [available here](https://www.ncbi.nlm.nih.gov/nuccore/AF157706)

## File tree and annotaitons

```
tree -L 2     
.
├── README.md
├── hhv6-reference # folder relevant to HHV-6B annotations / references used in this manuscript
│   ├── AF157706_sequence.txt
│   ├── HHV6b_expression_annotations.tsv 
│   ├── HHV6b_only.ec.txt
│   ├── HHV6b_only.index.txt
│   └── HHV6b_transcriptome.idx
├── public_bulk_data
│   ├── calderon-stim-immunecells-expression
│   ├── lamere-chipseq
│   ├── qu-atacseq
│   ├── richards-endothelial-CRS
│   ├── serratus-rnaseq-followup
│   ├── tcga # Preliminary analysis of OX40 expression in cancers
│   └── treg
└── single_cell_data # Data/analysis repositories for Fig. 2-4
    ├── broad-invivo # Folder for Fig 3 / cohort 2
    ├── cd7-cart # Folder for Fig 4
    ├── stanford-invitro # Folder for Fig 2
    ├── stjude-invitro # Folder for Fig 3 / in vitro reactivation
    └── stjude-invivo # Folder for Fig 3 / cohort 3
```

*Contact:* Raise a GitHub issue or email [Caleb](mailto:clareau@stanford.edu).

<br>
