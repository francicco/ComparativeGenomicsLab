# Part II: HAL tools and Alignemnt Manupulation
*Cactus* can take more than an hour to align this dataset on your instance, and it is strongly dependent on the allocated CPUs and memory for the job. I have provided the alignment in case you want to skip this step. Look for [HelicChr2.hal] among the files you downloaded.

The output of *Cactus* is a `HAL` ([Hierarchical Alignment Format](https://academic.oup.com/bioinformatics/article/29/10/1341/256598)) file. As the name suggests, the file is hierarchical and designed to store genome-scale multiple sequence alignments.
A `HAL` can represent a complete tree of genomes, including their sequence data, annotations, and evolutionary relationships. You can also think of a `HAL` file as some kind of database in binary format, meaning that you cannot directly access with a simple text editor. The `HAL` file is optimized for random access, as well as supporting incremental updates to the alignment. If needed, you can use the `Haltools` to manipulate and convert the alignment into other formats.


![bioinformatics_29_10_1341_f2](https://github.com/user-attachments/assets/1b64dec2-a786-4a46-9347-94bca92eecd3)
Figure. (A) A single genome as represented in `HAL`. Two sequences are stored in an array of DNA characters and are segmented with respect to its parent (top segments) and children (bottom segments). (B) The same genome in the context of HAL graph of five genomes. The dashed edge corresponds to an inversion event. Hickey G, Paten B, Earl D, Zerbino D, Haussler D (2013). [HAL: a hierarchical format for storing and analyzing multiple genome alignments](https://academic.oup.com/bioinformatics/article/29/10/1341/256598). Bioinformatics

In this section, you can utilize halSummarizeMutations to extract various useful information, such as mutations (insertions, deletions, inversions, duplications, transpositions, gap insertions, and gap deletions) in each branch of the alignment. The job is a single CPU process and may require ~25 minutes to complete. You can try redirecting the output to a log file and running it in the background using the '&' symbol.

1. What is the predominant mutation occurring at the branch of *Heliconius*?
2. And how many bases does it involve?

*Solution*
```bash
halSummarizeMutations HelicChr2.hal [optional --maxNFraction]
```

*NOTE:* In `halSummarizeMutations` you can specify `--targetGenomes` or `--rootGenome` option, `--maxNFraction 0` will instead prevent rearrangements with missing data as being identified as such. More generally, if an insertion of length 50 contains c N-characters, it will be labeled as missing data (rather than an insertion) if c/N > maxNFraction.

*Solution 1.* The predominant mutation should be *Insertions*, with over 1.6 Mb.

*NOTE:* This command can take a while, if you don't want to wait there's a hidden file named `.HelicChr2.SummarizeMutations.csv` you can run the `mv` command to make it visible:
```bash
mv .HelicChr2.SummarizeMutations.csv HelicChr2.SummarizeMutations.csv
```

You can also compute the alignment coverage using `halAlignmentDepth`, which returns the number of genomes aligned to a reference species of your choice. It will output data to the `Stdout` in wiggle format. You can redirect the output to a file and convert it into a compressed binary file (`BigWig` format) using `wigToBigWig` that can be easily visualized on a genome browser such as IGV.

1. Use the tool to compute the overall coverage using *H. melpomene* (`Hmel`) as a reference
```bash
halAlignmentDepth --noAncestors HelicChr2.hal Hmel > Hmel.Cov.wig
```

3. convert the wig file into a compress binary
```bash
wigToBigWig Hmel.Cov.wig Hmel.Chr2.fasta.fai Hmel.Cov.bw
```

4. Compute the coverage of only Heliconius species and convert the output
```bash
halAlignmentDepth --rootGenome Heliconius --noAncestors HelicChr2.hal Hmel > Hmel.HelOnly.Cov.wig
wigToBigWig Hmel.HelOnly.Cov.wig Hmel.Chr2.fasta.fai Hmel.HelOnly.Cov.bw
```

5. Compute the coverage of only non-Heliconius species
```bash
halAlignmentDepth --targetGenomes Eisa,Dpha,Smor --noAncestors HelicChr2.hal Hmel > Hmel.NonHelOnly.Cov.wig
wigToBigWig Hmel.NonHelOnly.Cov.wig Hmel.Chr2.fasta.fai Hmel.NonHelOnly.Cov.bw
```

6. Now load the file into IGV, can you spot any region of the chromosome that is more specific to Heliconius species?
   
