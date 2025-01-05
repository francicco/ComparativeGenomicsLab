# Part II: HAL tools and Alignemnt Manupulation
*Cactus* can take more than an hour to align this dataset on your instance, and it is strongly dependent on the allocated CPUs and memory for the job. I have provided the alignment in case you want to skip this step. Look for [HelicChr2.hal] among the files you downloaded.

The output of *Cactus* is a `HAL` ([Hierarchical Alignment Format](https://academic.oup.com/bioinformatics/article/29/10/1341/256598)) file. As the name suggests, the file is hierarchical and designed to store genome-scale multiple sequence alignments.
A `HAL` can represent a complete tree of genomes, including their sequence data, annotations, and evolutionary relationships. You can also think of a `HAL` file as some kind of database in binary format, meaning that you cannot directly access with a simple text editor. The `HAL` file is optimized for random access, as well as supporting incremental updates to the alignment. If needed, you can use the `Haltools` to manipulate and convert the alignment into other formats.


![bioinformatics_29_10_1341_f2](https://github.com/user-attachments/assets/1b64dec2-a786-4a46-9347-94bca92eecd3)
Figure. (A) A single genome as represented in `HAL`. Two sequences are stored in an array of DNA characters and are segmented with respect to its parent (top segments) and children (bottom segments). (B) The same genome in the context of HAL graph of five genomes. The dashed edge corresponds to an inversion event. Hickey G, Paten B, Earl D, Zerbino D, Haussler D (2013). [HAL: a hierarchical format for storing and analyzing multiple genome alignments](https://academic.oup.com/bioinformatics/article/29/10/1341/256598). Bioinformatics

In this section, you can utilize halSummarizeMutations to extract various useful information, such as mutations (insertions, deletions, inversions, duplications, transpositions, gap insertions, and gap deletions) in each branch of the alignment. The job is a single CPU process and may require ~25 minutes to complete. You can try redirecting the output to a log file and running it in the background using the '&' symbol.

1. What is the predominant mutation occurring at the branch of Heliconius?
2. And how many bases does it involve?


Create a folder where to run all the analyses:
```bash
mkdir -p $DIR && cd $DIR
```
