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

6. Now, using the phylogenetic tree and the extracted regions, build the neutral model using `phyloFit`. We can specify a few options, such as the substitution model. Usually between `REV` and `SSREV` models work well enough.
```bash
phyloFit --tree 7SpeciesPhylogeny.nex.treefile --subst-mod REV --sym-freqs --out-root neutralModel.4d HelicChr2.4d.maf
```
Have a look at the model:
```
ALPHABET: A C G T 
ORDER: 0
SUBST_MOD: REV
TRAINING_LNL: -862084.094730
BACKGROUND: 0.273505 0.226495 0.226495 0.273505 
RATE_MAT:
  -0.959160    0.165999    0.476776    0.316385 
   0.200452   -1.050496    0.280580    0.569464 
   0.575732    0.280580   -1.055584    0.199273 
   0.316385    0.471585    0.165022   -0.952993 
TREE: ((Smor:0.531082,(((Hdor:0.119252,Hmel:0.111786):0.049647,(Herd:0.0679488,Hcha:0.0799231):0.0520345):0.105919,Eisa:0.199437):0.138494):0.137492,Dpha:0.137492);
```
You can tell the model was done using nucleotides (`ALPHABET: A C G T `), with a substitution model `REV`, the likelihood of the training: `-862084.094730` , the rate matrix (`RATE_MAT`), and finally an adjusted phylogenetic tree with modified branch lengths.
<img width="1333" alt="ModelTree" src="https://github.com/user-attachments/assets/087aff22-1e9e-4405-9bed-0997d1c53df2" />



And this is the original tree:
<img width="1355" alt="OriginalTree" src="https://github.com/user-attachments/assets/365b2573-bbde-4b66-8ec7-d126475f3750" />



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

```bash
wc -l HelicChr2.mostcons.bed
16236 HelicChr2.mostcons.bed

wc -l HelicChr2.5b.Merged.mostcons.bed
15980 HelicChr2.5b.Merged.mostcons.bed

bc <<< 16236-15980
```

Now load all the generated data (`Chr2.PhastConsScores.bw`|`Hmel.phyloP.wig`|`HelicChr2.5b.Merged.mostcons.bed` ) on `IGV` and have a look. Do you see any pattern? Any region where you see lesser CEs in *H. melpomene*?

![Screenshot 2023-03-28 at 12 13 23 PM](https://github.com/user-attachments/assets/7e25c5ec-af96-4912-822b-c95eaec5a13a)



# Part IIIc: ATAC peak enrichment
If you examine the chromosome closely, you will notice that there is a significant amount of overlap between the conserved elements (CEs), phastCons scores, and PhyloP scores with coding regions. This is expected since coding exons are typically the most conserved regions of the genome. However, upon closer inspection, you will also find CEs within intergenic and intronic regions. In this particular case, these CEs span approximately 40 million years within the phylogenetic framework. It is important to note that most of these intergenic regions do not possess a valid open reading frame (ORF). Currently, it is unknown whether these regions represent functional elements such as enhancers or promoters. However, we can perform an assay for transposase-accessible chromatin (ATAC) to investigate whether there is an enrichment of chromatin accessibility with CEs. Indicating possible role as *cis*-regulatory elements (CREs).
important "switches" or "control switches" that influence when and how genes are turned on or off and, therefore expressed.

In this section, you can explore one possible application CE identification, checking for example if these regions correspond to promoters in your tissue of interest. With this goal, I have provided you with two annotations for ATAC peaks: one for *H. melpomene* (`melp_ros.MACS2.peaks.all.Chr2.counts`) and another for *H. erato demophoon* (dem_hyd.MACS2.peaks.all.Chr2.counts). The data was derived from the 5th instar heads and 5th instar, day 1, and day 2 pupal wings ([Ruggeri et al. 2022](https://genome.cshlp.org/content/32/10/1862.short)). The data contains the coordinates (start and end positions) of all ATAC-seq peaks plus the number of reads within those peaks (similar to gene expression counts in RNA-seq) for each sample/tissue.

![F5 large](https://github.com/user-attachments/assets/894bbd0e-f20e-47c0-a43f-170317e333f3)
*Figure. Example intervals of the genome alignment of *H. charithonia* (pink), *H. erato* (blue), and *H. melpomene* (orange) with the alignment of lineage-specific genome sequences, transposable element (TE) annotations, and ATAC-seq profiles. The plots show an illustrative interval of the genome assembly near the gene chiffon (chif) and sloppy paired 2 (slp2). For more details look at [Ruggeri et al. 2022](https://genome.cshlp.org/content/32/10/1862.short).*

To determine if the CEs correspond to promoters in these tissues, we could examine whether any of the CEs overlap with open chromatin regions, which are likely putative cis-regulatory elements, such as promoters.

1. In the first step, we need to prepare the data and convert the ATAC peaks files extracting peaks only relative to our chromosome of interest. We can easily do it with `grep`:
```bash
grep Hmel202001o melp_ros.MACS2.peaks.all.Chr2.counts > Hmel.Chr2.ATACpeaks.bed
```

2. Get the Conserved Non-Exonic Elements (CNEEs) and ATAC peaks that sit in intergenic regions [or in non-Exonic regions].
To do that you need to break the problems into smaller ones using `bedtools intersect`.

- Extract intergenic regions first:
```bash
grep -w gene Hmel.Chr2.annotation.gff3 | grep -v transcript | bedtools complement -i - -g Hmel.Chr2.fasta.fai > Hmel.IntergenicReg.bed
```

- Extract the ATAC peaks that are located in those regions...
```bash
bedtools intersect -u -wa -a Hmel.Chr2.ATACpeaks.bed -b Hmel.IntergenicReg.bed > Hmel.Chr2.ATACpeaks.InterGenic.bed
```

- ...and Do the same for intergenic CNEEs
```bash
bedtools intersect -u -wa -a HelicChr2.5b.Merged.mostcons.bed -b Hmel.IntergenicReg.bed > HelicChr2.5b.Merged.intergCEs.bed
```

3. Now check how many intergenic ATAC peaks overlap with CNEEs.
```bash
bedtools intersect -u -wb -a Hmel.Chr2.ATACpeaks.InterGenic.bed -b HelicChr2.5b.Merged.intergCEs.bed | wc -l
```

- What is their proportion?
```
Of the 1292 intergenic ATAC peaks 871 overlap with CEs, the ~67%
```

________________________________________________________________________________________________________________________________________________________________________________
## Extra assignment

Is this overlap significant in any way? That is, is there a significant enrichment in CNEEs where the chromatin is opened?

### - How could you test if there is a significant enrichment in CEs?

*Solution:*
You could perform a permutation test by reshuffling the ATAC peaks *n* times to estimate the expected overlap under the null hypothesis (*H₀*), which assumes that chromatin accessibility is randomly distributed in intergenic regions. By generating a distribution of overlaps across these permutations, you can evaluate where the observed overlap falls. If the observed overlap is within the lowest 5% of the null distribution, it would indicate a *p*-value less than 0.05, suggesting significant enrichment. Alternatively, you could use a binomial test to assess the probability of observing such an overlap under *H₀*.
 
`bedtools` have an implemented operator, `shuffle`, that randomly replaces the elements' coordinates of a bed file.

- Have a look at it and see if you can come up with a configuration that suits our needs.
- For each replica compute the overlap and store it in a file.
- Plot the results and what is the enrichment, if there is one.

*Solution:*
```bash
mkdir Permutation

echo -e "Perm\tOvlFrac" > Overlaps.dat

for REP in $(seq 1 1000); do
  echo -e "Rep: $REP"
  bedtools shuffle -excl Hmel.CDS.bed -chrom -chromFirst -noOverlapping -i Hmel.Chr2.ATACpeaks.InterGenic.bed \
    -g Hmel.Chr2.fasta.fai | bedtools sort -i - > Permutation/Hmel.Chr2.ATACpeaks.ReShuf.$REP.bed

  Nline=`wc -l Permutation/Hmel.Chr2.ATACpeaks.ReShuf.$REP.bed | cut -f1 -d' '`
  bedtools intersect -u -wa -a Permutation/Hmel.Chr2.ATACpeaks.ReShuf.$REP.bed -b HelicChr2.5b.Merged.mostcons.bed \
    | wc -l | awk '{ print "'"$REP"'""\t"$0/"'"$Nline"'"}' >> Overlaps.dat
done
```

You can run this script in R to plot the data and compute the binomial test
```Rscript
# Load the file
FracOvl <- read.table("Overlaps.dat", header=TRUE, sep="\t")

head(FracOvl)

# Statistics on the distribution
summary(FracOvl$OvlFrac)


# Calculate median and density
median <- summary(FracOvl$OvlFrac)[[3]]
density_default <- density(FracOvl$OvlFrac, adjust = 0.1)
density_smoothed <- density(FracOvl$OvlFrac, adjust = 2)

# Perform binomial test
x <- 871
n <- 1292
obsfreq <- 871/1292
binom_res <- binom.test(x, n, p = median)

# Extend x-axis range
x_range <- range(c(density_default$x, density_smoothed$x)) # Find the range of x-values
x_range <- c(x_range[1], obsfreq) # Extend the range by 1 on both sides


# Save the plot to a file
png("density_plot.png", width = 1200, height = 800)

# Create plot
plot(
  density_default, 
  main = "Density Plot with Smoothed Curve", 
  xlab = "Values", 
  ylab = "Density", 
  col = "blue", 
  lwd = 2,
  xlim = x_range # Extended x-axis range
)
lines(density_smoothed, col = "green", lwd = 2)
abline(v = median, col = "red", lty = 2, lwd = 2)
abline(v = obsfreq, col = "green", lty = 2, lwd = 2)
legend(
  "topright", 
  legend = c("Default Density (adjust = 0.1)", "Smoothed Density (adjust = 2)", "Expected freq.", "Observed freq."), 
  col = c("blue", "green", "red","green"), 
  lwd = 2, 
  lty = c(1, 1, 2, 2)
)

# Close the device
dev.off()
```

You should get a plot similar to this:
![density_plot](https://github.com/user-attachments/assets/ed217100-4b0e-420a-a1dc-96f0d7a79629)





________________________________________________________________________________________________________________________________________________________________________________
[Overview](https://github.com/francicco/ComparativeGenomicsLab/edit/main/README.md)

[Part I: Whole Genome Alignment](https://github.com/francicco/ComparativeGenomicsLab/blob/main/PartI/WholeGenomeAlignment.md)

Previous section: [Part II: HAL tools and Alignment Manipulation](https://github.com/francicco/ComparativeGenomicsLab/blob/main/PartII/AlignmentManipulation.md)
