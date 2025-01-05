# Part IIIa: Analysis of conserved regions and putative regulatory elements
Coding regions have been extensively analyzed for positive, neutral, or relaxing selection, primarily through the study of nonsynonymous vs synonymous nucleotide substitutions (omega = dN/dS) rates.
Non-coding regions, on the other hand, have often been neglected and not thoroughly explored. Whole genome alignment (WGA) provides a unique opportunity to investigate these genomic regions in various ways.
In this section, we will compute two metrics: the `PhyloP` **CONACC** and the identification of **Conserved Elements** (CEs) via `PhastCons` scores using a reference genome.
These metrics aim to identify functional elements - candidate *cis*-regulatory elements. Additionally, I am supplementing you with previously identified [ATAC (Assay for Transposase-Accessible Chromatin)](https://www.nature.com/articles/s41596-022-00692-9) peaks, for the study of chromatin accessibility, which in turn provides insights into the regulatory landscape of genomes.

The `PhyloP` **CONACC** metric employs the [likelihood ratio test (LRT)](https://en.wikipedia.org/wiki/Likelihood-ratio_test#:~:text=In%20statistics%2C%20the%20likelihood%2Dratio,the%20ratio%20of%20their%20likelihoods.) method to calculate conservation scores for each site in the alignment. The scores are then outputted in the fixed-step wiggle format. Positive values (log p-values or scores) indicate conservation, while negative values indicate acceleration.
On the other hand, `PhastCons` scores are computed using a [hidden Markov model-based method](https://en.wikipedia.org/wiki/Hidden_Markov_model) to estimate the probability of each nucleotide belonging to a conserved element (CE), based on the multiple alignments. This method takes into account not only the alignment of each nucleotide but also its flanking regions. Therefore, the PhastCons score represents the posterior probability of purifying selection and ranges between 0 and 1.

Both metrics rely on the [nucleotide substitution model](https://en.wikipedia.org/wiki/Substitution_model), which needs to be degenerate. Also known as *models of sequence evolution*, are models that describe how nucleotides or amino acids change over evolutionary time.
These substitution models are generally used to calculate the likelihood of phylogenetic trees using multiple sequence alignment data. Thus, substitution models are central to maximum likelihood estimation of phylogeny as well as Bayesian inference in phylogeny.

In our specific case, identifying genomic regions that vary under neutrality we will first estimate a neutral model and subsequently use it to identify genomic regions that deviate from neutrality and are evolving under acceleration or purifying selection (conservation).
There are various ways to identify regions evolving under neutrality, one way to compote this is considering the 4-fold degenerate sites; sites in which changes in the DNA sequence do not alter the amino acid sequence of the coded protein. These sites are more likely to be "neutral" or "silent" and less likely to have a functional impact on the protein.
Due to their neutrality, we can utilize these sites to estimate the rate of neutral mutations.

I provided a gene annotation of *H. melpomene* (`Hmel.Chr2.annotation.gff3`) in [gff3 format](https://www.ensembl.org/info/website/upload/gff3.html), our reference species, to extract at first coding regions and subsequently the 4-fold degenerative sites using `hal4dExtract` from HAL tools.

1. Get the coding regions using a combination of `gffread` and `grep -P`.
```bash
gffread --no-pseudo -C -T -o Hmel.CDS.gtf Hmel.Chr2.annotation.gff3
grep -P "\tCDS\t" Hmel.CDS.gtf > Hmel_filt.CDS.gtf
```

`gffread` will extract only coding transcripts in [gtf format](https://www.ensembl.org/info/website/upload/gff.html); while `grep` will extract only the coding exons.

2. Now we need to do some format conversion using `gtfToGenePred`...
```bash
gtfToGenePred Hmel_filt.CDS.gtf Hmel.CDS.gp
```

3. ...and finally into a `bed` file format with `genePredToBed`
```bash
genePredToBed Hmel.CDS.gp Hmel.CDS.bed
```

4. Now that we have a `bed` file we can use `hal4dExtract` from HAL tools to extract the site evolving under neutrality
```bash
hal4dExtract HelicChr2.hal Hmel Hmel.CDS.bed Hmel.4d.bed
```

Once we extracted the coordinates of those regions, we had to extract the relative alignment in a converted `MAF` format using `hal2maf`. `MAF` is a format for storing multiple sequence alignments.
It represents an alignment as a series of blocks, each of which contains one or more sequences and their alignment positions.

5. Extract those regions removing duplicated and ancestral genomes and using only ortholog sequences.
```bash
hal2maf --onlyOrthologs --noDupes --refGenome Hmel --refTargets Hmel.4d.bed HelicChr2.hal HelicChr2.4d.maf
```
--

6. Now, using the phylogenetic tree and the extracted regions, build the neutral model using `phyloFit`. We can specify a few options, such as the substitution model. Usually between `REV` and `SSREV` models work well enough.
(Hint: You can also try the `python` wrapper `halPhyloPTrain.py` which should be faster as parallelize the job).
```bash
phyloFit --tree 7SpeciesPhylogeny.nex.treefile --subst-mod SSREV --sym-freqs --out-root neutralModel.4d HelicChr2.4d.maf
```

```bash
halPhyloPTrain.py --substMod REV --noAncestors HelicChr2.hal Hmel Hmel.4d.bed neutralModel.mod --numProc 10
```


7. Finally using the neutral model, we can run `phyloP` on the whole alignment converted in `MAF` format. We can select whichever reference genome; internal branches are also accepted. Here you can use *H. melpomene* or *H. erato*.
Remember, you can compress the output converting the `wig` file into a `bigwig`.
```bash
hal2maf --refGenome Hmel --noAncestors --noDupes HelicChr2.hal HelicChr2.maf
```

(Alternative: `cactus-hal2maf ./js HelicChr2.hal HelicChr2.maf --refGenome Hmel --chunkSize 100000 --noAncestors`)

```bash
phyloP --mode CONACC --wig-scores --method LRT neutralModel.4d.mod HelicChr2.maf | sed 's/HelicChr2/Hmel202001o/' > Hmel.CONACC.phyloP.wig
```

```bash
wigToBigWig Hmel.CONACC.phyloP.wig Hmel.Chr2.fasta.fai Hmel.CONACC.phyloP.bw
```

You can use `halPhyloPMP.py`, a Python wrapper to leverage parallelisation.
```bash
halPhyloPMP.py --numProc 12 HelicChr2.hal Hmel neutralModel.4d.mod Hmel.phyloP.wig
```

In the same fashion we can also generate `PhastCons` scores and the relative CEs using `phastCons`. But before that, instead of a neutral model phastCons needs a model of the "conserved" state and another for "non-conserved" state. The model for conserved regions is optional. If it is not given, then this model is defined indirectly by "scaling" the model for non-conserved regions by a factor called *rho* (see the `--estimate-rho` option in `phastCons`).
The basic idea of the program is to scan along the alignment for regions that are better "explained" by the conserved model than by the non-conserved model; such regions will be output as CEs, and the probability that each base is in such a region will be output as the conservation score for that base.

8. Compute *rho*
```bash
phastCons HelicChr2.maf neutralModel.4d.mod --estimate-rho chr2 --expected-length 12 --target-coverage 0.25 > HelicChr2.phastCons.log
```

9. Use the two models obtained to generate a wiggle file and a bed file for the CEs (hint: use the option --most-conserved).
```bash
phastCons --score --most-conserved HelicChr2.mostcons.bed HelicChr2.maf chr2.cons.mod,chr2.noncons.mod | sed 's/HelicChr2/Hmel202001o/' > Chr2.PhastConsScores.wig
```

```bash
sed -i.BK 's/HelicChr2/Hmel202001o/' HelicChr2.mostcons.bed
```

```bash
wigToBigWig Chr2.PhastConsScores.wig Hmel.Chr2.fasta.fai Chr2.PhastConsScores.bw
```

These regions can be very close to each other, sometimes the distance is of just a few nucleotides. You can choose to merge these regions, let's say 5nt apart. Since the format of `PhastCons` is a `bed` file check if `bedtools` can be of any help.
```bash
bedtools merge -o first,mean -c 4,5 -d 5 -i HelicChr2.mostcons.bed > HelicChr2.5b.Merged.mostcons.bed
```

- How many regions did you merge?

Now load all the generated data (`Chr2.PhastConsScores.bw`|`Hmel.phyloP.wig`|`HelicChr2.5b.Merged.mostcons.bed` ) on `IGV` and have a look. Do you see any pattern? Any region where you see lesser CEs in *H. melpomene*?

![Screenshot 2023-03-28 at 12 13 23 PM](https://github.com/user-attachments/assets/7e25c5ec-af96-4912-822b-c95eaec5a13a)



# Part IIIc: ATAC peak enrichment
If you examine the chromosome closely, you will notice that there is a significant amount of overlap between the conserved elements (CEs), phastCons scores, and PhyloP scores with coding regions. This is expected since coding exons are typically the most conserved regions of the genome. However, upon closer inspection, you will also find CEs within intergenic and intronic regions. In this particular case, these CEs span approximately 40 million years within the phylogenetic framework. It is important to note that most of these intergenic regions do not possess a valid open reading frame (ORF). Currently, it is unknown whether these regions represent functional elements such as enhancers or promoters. However, we can perform an assay for transposase-accessible chromatin (ATAC) to investigate whether there is an enrichment of chromatin accessibility with CEs. Indicating possible role as *cis*-regulatory elements (CREs).
important "switches" or "control switches" that influence when and how genes are turned on or off and, therefore expressed.

In this section, you can explore one possible application CE identification, checking if for example these regions correspond to promoters in your tissue of interest. With this goal, I have provided you with two annotations for ATAC peaks: one for *H. melpomene* (`melp_ros.MACS2.peaks.all.Chr2.counts`) and another for *H. erato demophoon* (dem_hyd.MACS2.peaks.all.Chr2.counts). The data was derived from the 5th instar heads and 5th instar, day 1, and day 2 pupal wings ([Ruggeri et al. 2022](https://genome.cshlp.org/content/32/10/1862.short)). The data contains the coordinates (start and end positions) of all ATAC-seq peaks plus the number of reads within those peaks (similar to gene expression counts in RNA-seq) for each sample/tissue.

![F5 large](https://github.com/user-attachments/assets/894bbd0e-f20e-47c0-a43f-170317e333f3)
*Figure. Example intervals of the genome alignment of *H. charithonia* (pink), *H. erato* (blue), and *H. melpomene* (orange) with the alignment of lineage-specific genome sequences, transposable element (TE) annotations, and ATAC-seq profiles. The plots show an illustrative interval of the genome assembly near the gene chiffon (chif) and sloppy paired 2 (slp2). For more details look at [Ruggeri et al. 2022](https://genome.cshlp.org/content/32/10/1862.short).*

Let's check how many CEs overlap with the ATAC peaks, meaning that correspond to open chromatin regions (putative *cis*-regulatory elements), most likely promoters. 
In the first step, we need to Convert the ATAC peaks files into a valid bed file





