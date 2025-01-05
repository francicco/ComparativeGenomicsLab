# Comparative Genomics Lab - Workshop on Genomics, Cesky Krumlov

## Class Overview
In this part of the workshop, we will focus on certain aspects of comparative genomics. Specifically, our attention will be directed toward whole genome alignment (WGA) as a tool for investigating the changes occurring in different regions of the genome that interest us. Additionally, we will explore identifying potential regulatory elements in intergenic regions. All of this will be done in the context of the radiation of the neotropical Heliconius butterflies.


## Dependencies
Below is a breakdown of the required dependencies.
- [Kent toolkit](https://github.com/ucscGenomeBrowser/kent). All of the binaries required are available pre-compiled on the utility page. The tools needed for this tutorial are `wigToBigWig`, `gtfToGenePred`, `genePredToBed`.
- [Bedtools](https://bedtools.readthedocs.io/en/latest/) & [Samtools](http://www.htslib.org/)
- [Progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). With it, you should also find precompiled binaries for [Hal Tools](https://github.com/ComparativeGenomicsToolkit/hal).
- `gffread` from the [Cufflink](http://cole-trapnell-lab.github.io/cufflinks/) package.
- Phylogenetic Analysis with Space/Time Models ([PHAST](http://compgen.cshl.edu/phast/)) package.
- [IGV](https://software.broadinstitute.org/software/igv/) or any other genome browser you like, to have a look at all the tracks you will generate
Make sure you put the newly created ~/bin/$MACHTYPE directory on your path. 
![Phylogeny](https://github.com/user-attachments/assets/1c6becb6-bd86-4fcf-91e6-ab15b6f078d7)

## The workshop
The workshop is divided into four section

1. [RNAseq mapping on the reference genome (Short-reads & Iso-Seq)](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/1.Mapping/1.MappingStep.md)
2. [Homology and evidence-based prediction of protein coding genes (PCGs) using BRAKER2](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/2.Prediction/BRAKER.md)
3. [*De Novo* annotation using Short-reads RNAseq using Trinity](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/3.DeNovoAnnotation/DeNovoTrinity.md)
4. *Ab Initio* annotation using [Short-reads](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/4.AbInitioAnnotation/1.ShortReadAnnotation.md) RNAseq & [Long-reads RNA-Seq](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/4.AbInitioAnnotation/2.LongReadAnnotation.md)
5. Metrics to evaluate annotations, [Splice-site filtering](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/5.SpliceJunctionFiltering/5.PortcullisRun.md), and [Annotation Consensus using Mikado](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/6.Consensus/ConsensusAnnotationMikado.md)

All data is available at this [link](https://drive.google.com/drive/folders/1IreMRHaOa1kvOomyjoEm8xFw1fmOR-oK?usp=drive_link), but don't forget to [set up your environment](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/0.VariableSetting.md)!!!

## I hope this will be useful, Have fun!
