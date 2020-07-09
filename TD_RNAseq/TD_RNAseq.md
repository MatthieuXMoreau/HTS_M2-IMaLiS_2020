# Practical training – RNA-Seq data analysis


- [1. Introduction](#intro)
- [2. Connect to the server](#logging)
- [3. Sequence Quality Controls](#QC)
- [4. Mapping of Reads on the Reference Genome](#mapping)
- [5. Alignments Visualization with a Genome Browser](#genome_browser)
- [6. Search for Differentially Expressed Genes](#DEtest)

#

## Introdution <a name="intro"></a>

#

**Aim**: Getting started with bioinformatics tools and statistical approaches applied to analyze results coming from RNA-Seq technology. Data used in these practical were collected from the following publication:

>Guida, A., Lindstädt, C., Maguire, S. L., Ding, C., Higgins, D. G., Corton, N. J., Berriman, M., et al. (2011). Using RNA-seq to determine the transcriptional landscape and the hypoxic response of the pathogenic yeast Candida parapsilosis. BMC genomics
>[Guida* et al*. BMC Genomics 2011 ](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-628)

#

## Connect to the server <a name="logging"></a>

#

#### 1 - Sign in on the server
  * On Windows using [MobaXterm](https://mobaxterm.mobatek.net/)
  
> Session : ssh
> Host : core.cluster.france-bioinformatique.fr
> Specify username : ticked and filled in
> Advanced SSH settings : X11-Forwarding

  * On MacOS and Linux
```bash
ssh -XY <login>@core.cluster.france-bioinformatique.fr
```

#### 2 - Set up your working environment
1. Go to your project directory
```bash
cd /shared/projects/ens_HTseq_2020/<your login>/
```
2. Create a directory that will contain all results of the upcoming analyses.
```bash
mkdir RNAseq_Practical
```
3. Go to the newly created directory
```bash
cd RNAseq_Practical
```

4. Check your are in the right directory using `pwd`:

```bash
pwd

/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical
```

#

## Sequence Quality Controls <a name="QC"></a>

#

**FASTQ** files are raw results from RNA-Seq experiments. These files comprise all the **sequences** (or reads) obtained with the sequencer machine (Illumina technology here), together with **base quality scores** (PHRED scores).

Two different files will be analyzed in this practical :
- ***O2rep2_SRR352263.fastq*** refereed to a transcriptome analysis of yeasts C. parapsilosis under normoxic condition.
- ***noO2rep3_SRR352271.fastq*** refereed to a transcriptome analysis of yeasts C. parapsilosis under hypoxic condition (see [Guida* et al*.](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-628) for more information).

In a first step, quality controls will be perform on each FASTQ files in order to evaluate the quality of the sequences and identify potential problems that could impact the following analyses. Dedicated JAVA software will be used [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) . Note that other software exists.

FastQC is a quality control application for high throughput sequence data. It provides an interactive application to review the results of several different quality control checks, and create an HTML based report. The main functions of FastQC are:

>  - **Import of data** from BAM, SAM or FastQ files (any variant)
>   - Providing a **quick overview** to tell you in which areas there may be problems
>   - **Summary graphs** and tables to quickly assess your data
>   - Export of results to an **HTML based permanent report**
>   - **Offline operation** to allow automated generation of reports without running the interactive application


####  To do: Use FASTQC to evaluate the quality of sequences in each FASTQ files.


1. Create a new directory to store the output of fastqc

```bash
mkdir 1-QualityControl
```
Using the `tree` command, your directory should look like this :

```bash
/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical
│
└───1-QualityControl
```

2. Go to this directory

```bash
cd 1-QualityControl
```
3. Get Fastqc available in your environment

```bash
module add fastqc/0.11.8
```

4. Check the help page of the programme to see its usage and parameters 

```bash
srun fastqc --help
```

5. Run fastqc on each experiment files

- /shared/projects/ens_HTseq_2020/RNAseq/Fastqc/O2rep2_SRR352263.fastq : **absolute path** to the first file
- -o: creates all output files in the specified output directory. '.' means current directory.

```bash
# O2 condition reads
srun fastqc /shared/projects/ens_HTseq_2020/RNAseq/Fastqc/O2rep2_SRR352263.fastq -o .
```
```bash
# noO2 condition reads
srun fastqc /shared/projects/ens_HTseq_2020/RNAseq/Fastqc/noO2rep2_SRR352263.fastq -o .
```
At this point you should see the two new files in your directory using the `tree` command

```bash
/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical
│
└───1-QualityControl
	│
	└─── O2rep2_SRR352263.fastqc.zip
	│
	└───noO2rep2_SRR352263.fastqc.zip
```
6. Download the HTML file SRR576933_fastqc.html on your local machine

```bash
### OPEN A NEW TERMINAL
## Create a directory where to put generated files on your computer
mkdir ~/Desktop/RNAseq_Practical/

## Go to the location on your computer, where you want to put the data, for example:
cd ~/Desktop/RNAseq_Practical/

## Download the first file
scp <login>@core.cluster.france-bioinformatique.fr:/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical/1-QualityControl/O2rep2_SRR352263.fastqc.zip .
# Enter your password

## Download the second
scp <login>@core.cluster.france-bioinformatique.fr:/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical/1-QualityControl/O2rep2_SRR352263.fastqc.zip .
# Enter your password
```
7. Open a new shell and Unzip the files in your personal computer with `uzip`

```bash
unzip *.zip
```

8. Open the *.html* report with firefox

Using information from the [Fastqc help page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help) as well as exemples of [good](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)  and [bad](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) illumina data as references, compare results between the two FASTQ files.
**Is there any concern related to the following analyses?**

#

## Mapping of Reads on the Reference Genome <a name="mapping"></a>

#

Once data quality is verified, reads will be mapped onto the reference genome of yeast C. parapsilosis. The **genome sequence** of C. parapsilosis and its **annotations** (locations of ORFs) were retrieved from the [CGD database](http://www.candidagenome.org/ "Candidat genome database").

Different aligner and algorithms for RNA-Seq analysis exist. We will use [Bowtie1.2.2](http://bowtie-bio.sourceforge.net/index.shtml)  an ultrafast (memory-efficient) short read aligner. As an input, BOWTIE uses a FASTQ file (with reads to be aligned) and **“pre-built indexes”** of the reference genome. These indexes are named ***“C_parapsilosis.1.ebwt”***, ***“C_parapsilosis.2.ebwt”***, etc. They will allow boosting the alignment process.  
As an output, BOWTIE provides a SAM file. SAM (Sequence Alignment/Map) is a generic format for storing large nucleotide sequence alignments.

#### To do: Run sequence alignments with Bowtie using the two FASTQ files.

1. Create a new directory to store the output of fastqc

```bash
#Go to the parental directory "RNAseq_Practical"
cd ../

#Create a new directory to store results of the alignment
mkdir 2-Mapping
```
Your directory should now look like this :

```bash
/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical
│
└───1-QualityControl
	│
	└─── O2rep2_SRR352263.fastqc.zip
	│
	└───noO2rep2_SRR352263.fastqc.zip
│
└─── 2-Mapping
```

2. Go to this directory

```bash
cd 2-Mapping
```
3. Load Bowtie into your environment

```bash
module add bowtie/1.2.2
```

3. Map the reads to the reference genome

>- **-S** will output the result in SAM format
>- **/shared/projects/ens_HTseq_2020/RNAseq/bowtie_indexes/C_parapsilosis** specify the location and the **prefix (C_parapsilosis)** of the bowtie's index files
>- **/shared/projects/ens_HTseq_2020/RNAseq/Fastqc/O2rep2_SRR352263.fastq** location of the input fastq
>- **2>** will print some statistic about the aligment (#of reads mapped, etc...)
>- **>** redirects the mapping output into a .sam file

```bash
# Map the aerobic condition reads
srun bowtie -S /shared/projects/ens_HTseq_2020/RNAseq/bowtie_indexes/C_parapsilosis \
	/shared/projects/ens_HTseq_2020/RNAseq/Fastqc/O2rep2_SRR352263.fastq 2> O2rep2_SRR352263_bowtie_mapping.out > O2rep2_SRR352263_bowtie_mapping.sam
```

```bash
# Map the hypoxic condition reads
srun bowtie -S /shared/projects/ens_HTseq_2020/RNAseq/bowtie_indexes/C_parapsilosis \
 	O2rep2_SRR352263.fastq 2> noO2rep2_SRR352263_bowtie_mapping.out > noO2rep2_SRR352263_bowtie_mapping.sam
```

Your directory should now look like this :

```bash
/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical
│
└───1-QualityControl
	│
	└─── O2rep2_SRR352263.fastqc.zip
	│
	└─── noO2rep2_SRR352263.fastqc.zip
│
└─── 2-Mapping
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.sam
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.out
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.sam
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.out
```

**Take a look at [Bowtie documentation](http://bowtie-bio.sourceforge.net/manual.shtml) and describe the specified options (-m in particular).**
**What is the proportion of reads aligned on the reference genome?**

#

## Alignments Visualization with a Genome Browser <a name="genome_browser"></a>

#

The [Integrative Genomics Viewer](http://software.broadinstitute.org/software/igv/home) (IGV) is a high-performance **visualization tool** for interactive exploration of large, integrated genomic datasets. It supports a wide variety of data types, including array-based, next-generation sequence data and genomic annotations. In this practical, we will use IGV to visualize mapping results (see previous section). For that, **SAM files** has to be converted into **BAM files** (a binary version of SAM) and “sorted” according to the genomic sequence. We will use programs available in the [SAMTOOLS](http://samtools.sourceforge.net/) suite.

1. Sort and Converte *.sam* into *.bam* files

>- **samtools sort** Sort alignments by genomic coordinates
> - **|** "pipe" the output of samtools sort to the next programme i.e. samtools view
> - **samtools view** will convert sam into bam
> - **-b** specify the output to be in BAM format
>- **>** write the output in the bam file

```bash
module add samtools/1.9
```

```bash
# Sort and convert O2 condition
srun samtools sort O2rep2_SRR352263_bowtie_mapping.sam | srun samtools view -b  > O2rep2_SRR352263_bowtie_sorted.bam

# Sort and convert noO2 condition
srun samtools sort noO2rep2_SRR352263_bowtie_mapping.sam | srun samtools view -b  > noO2rep2_SRR352263_bowtie_sorted.bam
```
2. Create an index of the tow bam files

> IGV requieres the to have an index of the bam file. The index of a bam file is name ***.bam.bai***

```bash
#Index the O2 condition
srun samtools index O2rep2_SRR352263_bowtie_sorted.bam

#Index the noO2 condition
srun samtools index noO2rep2_SRR352263_bowtie_sorted.bam
```

3. Visualize mapping results with IGV

Once the IGV program is launched, it is necessary to **import the reference genome** “Genomes/Greate .genome File...” (see below). Select the FASTA file with the genomic sequence of C. parapsilosis “C_parapsilosis_CGD.fasta” (“Browse / FASTA file”) and **enter information regarding ORFs positions**, GFF file “C_parapsilosis_ORFs.gff” (“Browse Gene file”).  
Finally, give a name to your genome (“Unique identifier” and “Descriptive name”) and press “OK”. Save the genome file in your home.  
**Warning!** In order to IGV to create a index of your genome, you need to copy the reference genome FASTA file in writable directory.


<p align="center">

<img src="./images/IGV_genome.png" width="70%">

</p>

*C. parapsilosis* genome is now loaded into IGV and can be selected from the top/left menu (see 1 below). The genomic sequence can be therefore explored, choosing for instance, a particular chromosome (see 2 below) or a genomic region (see 3). Note that gene annotations (ORF positions) are shown on the bottom of the window (see 4, blue lines) and you can obtain a more detailed view of the sequence using the cursor located on the top/right of the window, see 5).

<p align="center">

<img src="./images/IGV_tools.png" width="50%">

</p>

Mapping results (“.sorted.bam” files) can now be imported (“File / Load from File”):

<p align="center">

<img src="./images/Reads_visualization.png" width="50%">

</p>


**Compare your results with those presented in the original publication :** 
**Did the authors use stranded-specific protocols ?**
**Can you observe differences between hypoxic and normoxic conditions ?**


#

## Search for Differentially Expressed Genes <a name="DEtest"></a>

#

To identify genes whose expression is different between hypoxic and normoxic conditions, we will **count and compare the number of reads mapped to each ORF**. A program available in the BEDTOOLS suite will be used.

#### Calculate for each ORF the number of reads that were aligned (normoxic and hypoxic conditions).

1. Create a new directory to store the ORF count matrix

```bash
#Go to the parental directory "RNAseq_Practical"
cd ../

#Create a new directory to store results of the alignment
mkdir 3-ORF_reads_count
```
Your directory should now look like this :

```bash
/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical
│
└───1-QualityControl
	│
	└─── O2rep2_SRR352263.fastqc.zip
	│
	└───noO2rep2_SRR352263.fastqc.zip
│
└─── 2-Mapping
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.sam
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.out
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.sam
	│
	└─── noO2rep2_SRR352263_bowtie_mapping.out
│
└─── 3-ORF_reads_count
```

2. Go to this directory

```bash
cd  3-ORF_reads_count
```

3. Calculate for each ORF the number of reads that were aligned to it

```bash
module add bedtools/2.27.1
```

```bash
srun bedtools multicov -bams ../O2rep2_SRR352263_bowtie_sorted.bam \
-bed /shared/projects/ens_HTseq_2020/RNAseq/C_parapsilosis_ORFs.gff > O2rep2_SRR352263_gene_counts.gff

srun sed 's/^.*ID=//' O2rep2_SRR352263_gene_counts.gff > O2rep2_SRR352263_gene_counts.tab
```

```bash
srun bedtools multicov -bams ../noO2rep2_SRR352263_bowtie_sorted.bam \
-bed /shared/projects/ens_HTseq_2020/RNAseq/C_parapsilosis_ORFs.gff > noO2rep3_SRR352271_gene_counts.gff

srun sed 's/^.*ID=//' noO2rep3_SRR352271_gene_counts.gff > noO2rep3_SRR352271_gene_counts.tab
```

4. Unload the tools you used

```bash
module rm samtools/1.9 bowtie/1.2.2 bedtools/2.27.1
```

#### Statistical analysis using DEseq2 R package.
In their article (Guida et al., 2011), the authors repeated the experiment 6 times for normoxic condition (with O2) and 4 times for hypoxic conditions (without O2). Results obtained for all experiments are combined in the file “/shared/projects/ens_HTseq_2020/RNAseq/R/count_data_diffAnalysis.txt”. This file will be used to search for differentially expressed genes using the **DESeq2** ([Love *et al*. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)) R package. The [DESeq package](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) provides methods to test for differential expression by use of the negative binonial distribution and a shrinkage estimator for the distribution’s variance.

1. Connect to Rstudio serveur of the IFB

In a web browser, connect to https://rstudio.cluster.france-bioinformatique.fr/auth-sign-in and log in using your user name and pasword (same as for ssh connection)

![Rstudio_login](./images/Rstudio.png)

You will reached the familiar Rstudio environment :

<p align="center">

<img src="./images/RstudioScreen.png" width="50%">

</p>

2. Save the working notebook in your personal environment

In *File > Open File...* enter the path ***/shared/projects/ens_HTseq_2020/RNAseq/R/DEseq2.Rmd*** to open the notbook containing all the code needed for the practical.
Save it into your personal folder using *File > Save As* and the path ***/shared/projects/ens_HTseq_2020/<your login>/RNAseq_Practical/DEseq2.Rmd***

3. Conduct statistical analysis in R

At this point you can use the notebook to run each chunk of code and conserve all the output along with the script


----

```r
# library loading
library(DESeq2)
```

```r
# data reading
countData = read.table("/shared/projects/ens_HTseq_2020/RNAseq/R/count_data_diffAnalysis.txt", row.names = 1, header = TRUE)

# general information
row.names(countData)
head(countData)
```

```r
# Loading metadata for the experiment
# N = Normoxic condition (with O2)
# H = Hypoxic condition  (without O2)
colData = read.table("/shared/projects/ens_HTseq_2020/RNAseq/R/design.txt", row.names = 1, header = TRUE)

```

```r
# DESeqDataSet object creation
dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
head(dds)

```

```r
# normalization of counts
# calculation of sizeFactors
dds = estimateSizeFactors(dds)
sizeFactors(dds)
```

```r
# bar plots
barplot(colSums(counts(dds)))
```

```r
# boxplots
boxplot(log2(counts(dds)+1))
boxplot(log2(counts(dds,normalized=TRUE)+1))
```

```r
# variance estimations
dds = estimateDispersions(dds)
plotDispEsts(dds)
```

```r
# differential analysis
dds = nbinomWaldTest(dds)
res = results(dds)
mcols(res,use.names=TRUE)
resultsNames(dds)
res <- as.data.frame(res)
```


```r
# diagnostic plots
attach(res)

#volcano plot
plot(log2FoldChange, -log10(padj),pch =21, main = "Normoxic VS Hypoxic")
abline(v = c(-2,2), col ="red")
abline(h = 1.3, col= "red", lty =2)
points(log2FoldChange[log2FoldChange >2 & padj < 0.05], -1 * log10(padj[log2FoldChange >2  & padj < 0.05]), col ="green")
points(log2FoldChange[log2FoldChange <(-2) & padj < 0.05], -1 * log10(padj[log2FoldChange <(-2)  & padj < 0.05]), col ="red")
legend("topright", c("107 genes with LogFc  > 2 in Hypoxic VS normoxic", " 85 genes with LogFc  > 2 in normoxic VS hypoxic"),pch = 21, col = c("red", "green"), bty ="n", cex =.9)
```

```r
FoldFc_sup2 <- res[res[,2] < (-2) & res[,5] <0.05,]
nrow(FoldFc_sup2)
write.table(FoldFc_sup2, row.names = T, quote = F, sep = "\t", file = "../results/RNAseq_FoldFc_sup2.txt")

```

```r
FoldFc_inf2 <- res[res[,2] <(-2),]
nrow(FoldFc_inf2)
```

```r
# MA plot
plotMA(res)
```

```r
# Histogram of the p values
hist(res$pvalue,breaks=20,col="grey")
# PCA plots
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
p = plotPCA(vsd)
p = update(p, panel = function(x, y, ...) {lattice::panel.xyplot(x, y, ...);lattice::ltext(x=x, y=y, labels=rownames(colData(vsd)), pos=1, offset=1, cex=0.8)})
print(p)
```

```r
# genes are sorted according to adjusted p-values
res = res[order(res$padj),]
dim(res)
res[rownames(res) == "CPAR2_212440",]

write.table(res, "res.txt", row.name=T, quote=F, sep='\t')
```

**Search for differentially expressed genes using DESeq R package.**
**How many genes are selected with different p-value thresholds (5%, 1%, etc.) ?**
**Check your results with IGV and use GOtermFinder (see practical on microarrays) to analyse the function of the selected genes.**
