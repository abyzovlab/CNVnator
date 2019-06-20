# README 

## Quick start guide

```
# Extract read mapping
$ ./cnvnator -root file.root -tree file.bam

# Generate histogram
$ ./cnvnator -root file.root -his 1000 -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -d dir_with_genome_fa/
  OR
$ ./cnvnator -root file.root -his 1000 -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y -fasta file_genome.fa.gz
  OR
$ ./cnvnator -root file.root -his 1000 -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -fasta file_genome.fa.gz

# Calculate statistics
$ ./cnvnator -root file.root -stat 1000 

# Partition
$ ./cnvnator -root file.root -partition 1000

# Call CNVs
$ ./cnvnator -root file.root -call 1000
```

## 1. Compilation

### Dependencies

You must install [ROOT package](http://root.cern.ch) and set up `$ROOTSYS` variable (see ROOT documentation [here](https://root.cern.ch/root/html534/guides/users-guide/GettingStarted.html)).

Also, a link to the samtools binary should be present in your CNVnator directory.

If compilation is not completed but the file libbam.a has been created, you can continue.

### Installation from github

```
git clone https://github.com/abyzovlab/CNVnator.git

cd CNVnator

ln -s /path/to/src/samtools samtools

make
```

If make doesn't work, try `"make OMP=no"` which will disable parallel support.

### Installing with Yeppp support

[Yeppp](http://www.yeppp.info/) is a library which provides high-performance implementations of math functions.

To install with Yeppp support, download Yeppp from [here](http://bitbucket.org/MDukhan/yeppp/downloads/yeppp-1.0.0.tar.bz2)
and extract it to a location of your choice. Set `YEPPPLIBDIR` and `YEPPPINCLUDEDIR` directories appropriately. 

Typically, for Linux-based systems on x86-64, `YEPPPLIBDIR` will be yeppp-1.0.0/binaries/linux/x86_64/ and `YEPPPINCLUDEDIR` will be
yeppp-1.0.0/library/headers. 

To build, type  
`make YEPPPLIBDIR=... YEPPPINCLUDEDIR=...`  

To disable OpenMP, add `OMP=no` to the make command.

## 2. Predicting CNV regions

Running CNVnator involves a few steps outlined below. Chromosome names and lengths are
parsed from the input sam/bam file header. 

### 2.1 EXTRACTING READ MAPPING FROM BAM/SAM FILES

```
$ ./cnvnator -root out.root [-chrom name1 ...] -tree [file1.bam ...] [-lite]
```
where,


-root out.root  -- specifies output ROOT file. See ROOT package documentation.  
-chrom name1 ... -- specifies chromosome name(s).  
-tree file1.bam ...  -- specifies bam file(s) names.
-lite -- use this option to produce a "lighter" (smaller) root file.


Chromosome names must be specified the same way as they are described in the sam/bam
header, e.g., chrX or X. One can specify multiple chromosomes separated by
space. If no chromosome is specified, read mapping is extracted for all chromosomes
in the sam/bam file. Note that this would require machines with a large physical
memory of at least 7Gb. Extracting read mapping for subsets of chromosomes is a way
around this issue. Also note that the root file is not being overwritten.


Example:

```
./cnvnator -root NA12878.root -chrom 1 2 3  -tree NA12878_ali.bam
```

for bam files with a header like this:  
@HD VN:1.4    GO:none  SO:coordinate  
@SQ SN:1      LN:249250621  
@SQ SN:2      LN:243199373  
@SQ SN:3      LN:198022430  
...  

or

```
./cnvnator -root NA12878.root -chrom chr1 chr2 chr3 -tree NA12878_ali.bam
```
for bam files with a header like this:  
@HD VN:1.4    GO:none  SO:coordinate  
@SQ SN:chr1   LN:249250621  
@SQ SN:chr2   LN:243199373  
@SQ SN:chr3   LN:198022430  
...  

Example:

```
./cnvnator -root NA12878.root -chrom 4 5 6 -tree NA12878_ali.bam
./cnvnator -root NA12878.root -chrom 7 8 9 -tree NA12878_ali.bam
```

is equivalent to

```
./cnvnator -root NA12878.root -chrom 4 5 6 7 8 9 -tree NA12878_ali.bam
```


### 2.2 GENERATING A READ DEPTH HISTOGRAM

```
$ ./cnvnator -root file.root [-chrom name1 ...] -his bin_size [-d dir]
```

This step is not memory consuming and so can be done for all chromosomes at once. 
It can also be carried for a subset of chromosomes. 
Files with individual chromosome sequences (.fa) are required and should reside in the 
current directory or in the directory specified by the -d option. 
Files should be named as: chr1.fa, chr2.fa, etc.


### 2.3 CALCULATING STATISTICS

```
$ ./cnvnator -root file.root [-chrom name1 ...] -stat bin_size
```

This step must be completed before proceeding to partitioning and CNV calling.


### 2.4 RD SIGNAL PARTITIONING

```
$ ./cnvnator -root file.root [-chrom name1 ...] -partition bin_size [-ngc]
```

Option `-ngc` specifies not to use GC corrected RD signal. Partitioning is the most time consuming step.

### 2.5 CNV CALLING

```
$ ./cnvnator -root file.root [-chrom name1 ...] -call bin_size [-ngc]
```

Calls are printed to STDOUT by default. You may redirect them to a file using the redirect operator >

The output columns are as follows:

CNV\_type coordinates CNV\_size normalized\_RD e-val1 e-val2 e-val3 e-val4 q0

where,

normalized_RD -- read depth normalized to 1.  
e-val1        -- is calculated using t-test statistics.  
e-val2        -- is from the probability of RD values within the region to be in
the tails of a gaussian distribution describing frequencies of RD values in bins.  
e-val3        -- same as e-val1 but for the middle of CNV  
e-val4        -- same as e-val2 but for the middle of CNV  
q0            -- fraction of reads mapped with q0 quality


### 2.6 REPORTING READ SUPPORT

To find and report read support for deletions and duplications by abnormal read pairs, use the -pe option as below:

```
./cnvnator -pe file1.bam ... -qual val(20) -over val(0.8) [-f file]
```

Once prompted, enter a genomic region and the CNV type, e.g.,

```
>12:11396601-11436500 del
or
>chr12:11396601-11436500 del
```

Please note that the bin size should be equal to a whole number of 100 bases (e.g., 2500, 3700,…)

### 2.7 MERGING ROOT FILES

```
./cnvnator -root out.root [-chrom name1 ...] -merge file1.root ...
```
Merging can be used when combining read mappings extracted from multiple files.  
Note: histogram generation, statistics calculation, signal partitioning, and
CNV calling should be completed/redone after merging.


## 3. Importing VCF data

To import variant data from VCF file use following option:

```
./cnvnator -root file.root [-chrom name1 ...] [-rmchr | -addchr] -vcf file.vcf.gz
```

If chromosome names are not specified, data for all chromosomes from file.vcf.gz will be imported. If 
you would like to add or remove the "chr" prefix from your chromosome names, use options `-addchr` or `-rmchr` respectively. 
It is important that chromosome names in the vcf file and the SAM/BAM file match. 

To mark known SNPs from the SNP database:

```
./cnvnator -root file.root [-chrom name1 ...] [-rmchr | -addchr] -idvar databasefile.vcf.gz
```

On running the above line, each SNP will be associated with a binary flag which equals 1 if it's in the database.


To mark variants based on genome accessibility using mask file from the 1000 Genomes Project:

```
./cnvnator -root file.root [-chrom name1 ...] [-rmchr | -addchr] -mask maskfile.fa.gz
```

On running the above line, each SNP will be associated with a binary flag which equals 1 if it's in the P-region


## 4. Genotyping genomic regions and visualization

For efficient genotype calculations, we recommend that you sort the list of regions by
chromosomes.

```
./cnvnator -root file.root -genotype bin_size [-ngc]
```

Once prompted enter a genomic region, e.g., 

```
>12:11396601-11436500
 or
>chr12:11396601-11436500
 or 
>12 11396601 11436500
 or
>chr12 11396601 11436500
```

One can also perform instant visualization by adding the word 'view', e.g.,

```
>12:11396601-11436500 view
 or
>chr12:11396601-11436500 view
 or
>12 11396601 11436500 view
 or
>chr12 11396601 11436500 view
```

### Additional notes

For genotyping of multiple regions one can use input piping, e.g.,
```
./cnvnator -root NA12878.root -genotype 100 << EOF
12:11396601-11436500
22:20999401-21300400
exit
EOF
```

Another example:
```
awk '{ print $2 } END { print "exit" }' calls.cnvnator | ./cnvnator -root NA12878.root -genotype 100
```


### 4.1 Visualizing specified regions

```
./cnvnator -root file.root [-chrom name1 ...] -view bin_size [-ngc]
```

Once prompted, enter a genomic region, e.g.,

```
>12:11396601-11436500
 or
>chr12:11396601-11436500
 or 
>12 11396601 11436500
 or
>chr12 11396601 11436500
```

Additionally, one can specify the length of flanking regions (default is 10 kb) to
be displayed as well, e.g.,

```
>12:11396601-11436500 100000
 or
>chr12:11396601-11436500 100000
 or
>12 11396601 11436500 100000
 or
>chr12 11396601 11436500 100000
```

One can also perform instant genotyping by adding the word 'genotype', e.g.,

```
>12:11396601-11436500 genotype
 or
>chr12:11396601-11436500 genotype
 or
>12 11396601 11436500 genotype
 or
>chr12 11396601 11436500 genotype
```

### 4.2 Plotting B-allele frequency (BAF)

To plot BAF data along RD use baf option in view mode:

```
./cnvnator -root file.root -view bin_size

>1:1-200000000 baf
```

The resulting output plot has two panels. On the uper panel, black line corresponds to binned RD
signal, green to segmentation, and red to calls. On the bottom panel each dot corresponds to BAF 
value of the SNPs. Colors represent following:
 
* black - homozygous (1/1 or 1|1) SNPs in P-region of the strict mask,
* grey - homozygous (1/1 or 1|1) SNPs out of P-region of the strict mask,
* blue - heterozygous (0/1 or 0|1) SNPs in P-region of the strict mask,
* cyan - heterozygous (0/1 or 0|1) SNPs out of P-region of the strict mask,
* red - heterozygous (1|0) SNPs in P-region of the strict mask,
* orange - heterozygous (1|0) SNPs out of P-region of the strict mask.

Plot BAF data with python tool plotbaf.py (requires numpy, matplotlib installed):

```
./plotbaf.py [-h] [-bs BINSIZE] [-res RESOLUTION] [-o SAVE_FILE] [-t TITLE]
             [-nomask] [-useid] root_file region
```

required arguments:
root\_file: cnvnator root file name
region: chromosomal coordinates in the format chr:start-end


optional arguments: 
size of bins: -bs BINSIZE, --binsize BINSIZE
likelihood function resolution: -res RESOLUTION, --resolution RESOLUTION
save plot to file: -o SAVE_FILE, --save_file SAVE_FILE
plot title: -t TITLE, --title TITLE
do calculations without mask: -nomask
do calculations using idvar filter: -useid

Output plot consists of four panels. Starting from the top one, they are:

* BAF value for heterozygous SNPs.
* Likelihood function. Light dots on the imagemap represent the most likely value of BAF at each bin.
* Red line represents a distance between maxima positions in likelihood function that is equivalent
  to twice the absolute difference between most likely BAF value and 0.5. Blue dots represent the ratio
  between the value of the likelihood function at 0.5 and its maximum value.
* Green dots and blue error-bars correspond to mean MAF and standard deviation per bin,
  respectively. Bin size is 100k base pairs.


## 5. Exporting CNV calls as VCFs

In order to export your CNV calls as a VCF file, use the script `cnvnator2VCF.pl` as

```
cnvnator2VCF.pl -prefix study1 -reference GRCh37 sample1.cnvnator.out /path/to/individual/fasta_files
```

where, 

-prefix specifies a prefix string you want to append to the ID field in your output VCF. For e.g., if you set your -prefix as "study1", then your resulting ID column will be study1_CNVnator_del_1, study1_CNVnator_del_2 etc.

-reference stands for the name of reference genome you used, for e.g., GRCh37, hg19 etc.

file.calls is your CNVnator output file with the CNV calls

genome_dir is the directory containing your individual reference fasta files such as 1.fa, 2.fa etc. (or chr1.fa, chr2.fa etc.)


## Contact Us

Please send your comments and suggestions to abyzov.alexej@mayo.edu
