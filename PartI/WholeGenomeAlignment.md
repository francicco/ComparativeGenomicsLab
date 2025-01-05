# Part I: Whole Genome Alignemnt
As a first step in this analysis, we will generate the WGA using the [*Progressive Cactus*](https://www.nature.com/articles/s41586-020-2871-y) aligner. *Cactus* is a highly accurate reference-free multiple genome alignment program.
*Progressive Cactus* enables the alignment of tens to thousands of large genomes without the use of a reference while maintaining high alignment quality and scalability.
*Cactus* solves large alignment problems by breaking them down into smaller subproblems using an input guide tree. Each subproblem involves comparing a set of ingroup genomes (the children of the internal node to be reconstructed) against each other, as well as a sample of outgroup genomes (non-descendants of the internal node in question).

**Dataset:** For my study, I utilized a total of 63 genomes. Among them, 44 genomes belong to *Heliconius* species, 14 genomes belong to other Heliconiini species, and 5 genomes are used as outgroup references. Due to time constraints, I selected 7 species from the complete dataset: *H. melpomene*, *H. doris*, *H. erato*, *H. demophoon*, and *H. charitonia* as ingroup genomes; *Eueides isabella*, *Dryadula phaetusa*, and *Speyeria mormonia* as outgroup genomes.
Furthermore, I narrowed down the dataset to a single chromosome. The ancestral chromosome 21, which in *Heliconius* species was renamed chromosome 2.

The documentation for Progressive Cactus can be found on the [GitHub repository](https://github.com/ComparativeGenomicsToolkit/cactus). In this workshop, we will align the 7 selected species, but first, we have to instruct *Cactus* on how to do it. Meaning that we have to tell the location of the fasta files and the phylogenetic tree that cactus will use as a guide to progressively decompose the alignment by iteratively aligning sibling genomes to estimate their parent genomes in a bottom-up fashion.
Cactus utilizes the predicted branch lengths from the tree to determine appropriate pairwise alignment parameters, enabling quicker alignment for closely related species without sacrificing accuracy.

The configuration file needs to be prepared by the user and it's a simple plan text file that should be formatted as follows:
```
NEWICK tree
name1 path1
name2 path2
...
nameN pathN
```

This is an example that you can find in the *Cactus* repository for some primate species
```
  # Sequence data for progressive alignment of 4 genomes
  # human, chimp and gorilla are flagged as good assemblies.  
  # since orang isn't, it will not be used as an outgroup species.
 (((human:0.006,chimp:0.006667):0.0022,gorilla:0.008825):0.0096,orang:0.01831);
 *human /data/genomes/human/human.fa
 *chimp /data/genomes/chimp/
 *gorilla /data/genomes/gorilla/gorilla.fa
 orang /cluster/home/data/orang/
```

Note: The asteriscs indicate genomes with good contiguity that will be used by *Cactus* as "better references"

I also provide you with a phylogenetic tree `7SpeciesPhylogeny.nex.treefile` (in newick format) and seven fasta files. You can use [`FigTree`](http://tree.bio.ed.ac.uk/software/figtree/) or [`TreeViewer`](https://treeviewer.org/) to have a look at the tree (optional).

![PhylogenyTutorial](https://github.com/user-attachments/assets/b4b53126-ebc0-4644-b851-ded8718924c1)
*Figure. Maximum likelihood phylogenetic tree of the seven selected species. Leaf names correspond to species tags that need to be reported in the cactus configuration file/*

Now:
1. Prepare the config file
2. Activate *Cactus*
3. Execute *Cactus*


*Solution:*
```bash
source $HOME/software/.source/cactus-bin-v2.5.1/cactus_env/bin/activate
cactus jobStore Cactus.Chr2.config Chr2.hal --workDir . --maxMemory 20G --binariesMode local --stats --logFile Cactus.log --maxCores 2 --defaultMemory 5G
```


________________________________________________________________________________________________________________________________________________________________________________
Previous section: [Overview](https://github.com/francicco/ComparativeGenomicsLab/edit/main/README.md)

Next section: [Part II: HAL tools and Alignment Manipulation](https://github.com/francicco/ComparativeGenomicsLab/blob/main/PartII/AlignmentManipulation.md)
