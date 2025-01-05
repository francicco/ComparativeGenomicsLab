# Comparative Genomics Lab - Workshop on Genomics, Cesky Krumlov

## Class Overview
In this part of the workshop, we will focus on certain aspects of comparative genomics. Specifically, our attention will be directed toward whole genome alignment (WGA) as a tool for investigating the changes occurring in different regions of the genome that interest us. Additionally, we will explore identifying potential regulatory elements in intergenic regions. All of this will be done in the context of the radiation of the neotropical *Heliconius* butterflies.


## Dependencies
Below is a breakdown of the required dependencies.
- [Kent toolkit](https://github.com/ucscGenomeBrowser/kent). All of the binaries required are available pre-compiled on the utility page. The tools needed for this tutorial are `wigToBigWig`, `gtfToGenePred`, `genePredToBed`.
- [Bedtools](https://bedtools.readthedocs.io/en/latest/) & [Samtools](http://www.htslib.org/)
- [Progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). With it, you should also find precompiled binaries for [Hal Tools](https://github.com/ComparativeGenomicsToolkit/hal).
- `gffread` from the [Cufflink](http://cole-trapnell-lab.github.io/cufflinks/) package.
- Phylogenetic Analysis with Space/Time Models ([PHAST](http://compgen.cshl.edu/phast/)) package.
- [IGV](https://software.broadinstitute.org/software/igv/) or any other genome browser you like, to have a look at all the tracks you will generate.

![Phylogeny](https://github.com/user-attachments/assets/1c6becb6-bd86-4fcf-91e6-ab15b6f078d7)
*Figure. The dated species phylogeny built from the concatenated single-copy orthologous groups (scOGs) from all sequenced Heliconiinae and outgroups, using a combination of Maximum Likelihood and Bayesian Inference. The branch color represents the number of substitutions per site per 100 Mya of that specific branch. Species names in bold indicate the species with chromosome- or sub-chromosome-level assemblies, asterisks indicate genomes assembled in this study, C curated assemblies; ii) genome assembly size, in red the TE fractions; iii) BUSCO profiles for each species. Blue indicates the fraction of complete single-copy genes; iv) bar plots show total gene counts partitioned according to their orthology profiles, from Nymphalids to lineage-restricted and clade-specific genes. From [Cicconardi et al. (2023)](https://www.nature.com/articles/s41467-023-41412-5).*

## The tutorial
The tutorial is divided into three sections

1. [Part I: Whole Genome Alignment](https://github.com/francicco/ComparativeGenomicsLab/blob/main/PartI/WholeGenomeAlignment.md)
2. [Part II: HAL tools and Alignment Manipulation](https://github.com/francicco/ComparativeGenomicsLab/blob/main/PartII/AlignmentManipulation.md)
3. [Part III: Identification of Conserved Regions](https://github.com/francicco/ComparativeGenomicsLab/blob/main/PartIII/IndentificationConservedElements.md)

## All data is available at this [link](https://uob-my.sharepoint.com/:f:/g/personal/tk19812_bristol_ac_uk/El4csr5H5jpHvVBhL3OVNZIB63COCfld3kpyB3FHzeAR_g?e=eAf5d3).

## I hope this will be useful, Have fun!
