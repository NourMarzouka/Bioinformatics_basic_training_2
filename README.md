# Understanding VCF Files and Focus on Specific Regions

This guide provides a simple introduction to Variant Call Format (VCF) files and how they are used in analyzing genetic variations, with a focus on specific regions of chromosomes. It's designed for school students interested in bioinformatics and genomics.

## What are VCF Files?

VCF (Variant Call Format) files are text files that store information about genetic variations found in DNA sequences. Think of them as a detailed report card for a genome, highlighting all the differences from a reference sequence.

Each line in a VCF file typically represents a single genetic variation, such as a Single Nucleotide Polymorphism (SNP) or an Insertion/Deletion (INDEL).

### Key Information in a VCF File:

*   **CHROM:** The chromosome or contig where the variation is located.
*   **POS:** The position on the chromosome where the variation starts.
*   **ID:** An identifier for the variation (e.g., rsID for SNPs).
*   **REF:** The reference allele (the DNA base(s) found in the reference genome).
*   **ALT:** The alternative allele(s) (the DNA base(s) found in the sample that differ from the reference).
*   **QUAL:** Quality score for the variation call (higher is better).
*   **FILTER:** Indicates if the variation passed certain quality filters.
*   **INFO:** Additional information about the variation (e.g., allele frequency, depth of coverage).
*   **FORMAT:** Describes the format of the genotype information for each sample.
*   **SAMPLE(s):** Genotype information for each individual sample.

## Analyzing VCF Files: Tools and Concepts

Analyzing VCF files helps scientists understand genetic differences between individuals or populations, which can be important for studying diseases, ancestry, and evolution.

### 1. `bcftools`

`bcftools` is a powerful command-line tool used for manipulating and analyzing VCF and BCF (binary VCF) files. It can perform many tasks, including:

*   **Filtering:** Selecting variations based on quality, type, or location.
*   **Merging:** Combining VCF files from different samples or experiments.
*   **Stats:** Generating summary statistics about variations (e.g., number of SNPs, INDELs, quality distributions).

**Example `bcftools` subfunctions and usage with `SRR2584866_variants.vcf.gz`:**

#### Step 1: Get VCF Statistics

To get a summary of statistics from your VCF file:

```bash
bgzip SRR2584866_variants.vcf
tabix SRR2584866_variants.vcf.gz
bcftools stats SRR2584866_variants.vcf.gz
```

#### Step 2: Filter VCF by Quality

To filter variants with a quality score greater than 30:

```bash
bcftools view -i "QUAL>30" SRR2584866_variants.vcf.gz > high_quality_variants.vcf
```

#### Step 3: Focus on specific region

```bash
tabix SRR2584866_variants.vcf.gz CP000819.1:1500-2000
```
or we can use this also

```bash
bcftools view SRR2584866_variants.vcf.gz -r CP000819.1:1500-2000 > CP000819.1_1500_2000_variants.vcf
```

### 2. `samtools`

While `samtools` is primarily used for Sequence Alignment/Map (SAM) and Binary Alignment/Map (BAM) files (which store aligned sequencing reads), it's often used in conjunction with `bcftools` because VCF files are generated from BAM files. `samtools` can provide statistics on the alignment process, which helps in understanding the quality of the data used to call variants.

**Example `samtools` subfunctions and usage (assuming `SRR2584866.aligned.sorted.bam` exists):**

#### Step 4: View BAM Header

To view the header of a BAM file:

```bash
samtools view -H SRR2584866.aligned.sorted.bam
```

#### Step 5: Count Reads in BAM

To count the number of reads in a BAM file:

```bash
samtools view -c SRR2584866.aligned.sorted.bam
```

#### Step 6: Subset BAM by Region

To extract aligned reads for a specific region (e.g., `CP000819.1:1500-2000`) from a BAM file:

```bash
samtools view -b SRR2584866.aligned.sorted.bam CP000819.1:1500-2000 > SRR2584866_subset.bam
```

##  Focusing on Specific Regions

Instead of analyzing an entire genome, scientists often focus on specific chromosomes or even smaller regions of interest. This is particularly useful when studying diseases linked to particular genes or chromosomal abnormalities.

### Why focus on specific regions?

*   **Targeted Analysis:** Allows for deeper investigation of areas known or suspected to be involved in a particular trait or disease.
*   **Reduced Computational Load:** Analyzing smaller datasets is faster and requires less computing power.
*   **Clinical Relevance:** Many genetic disorders are associated with specific chromosomal regions (e.g., Down syndrome with chromosome 21, certain cancers with specific genes on particular chromosomes).

#### Step 7: Subset FASTQ files
# Subset the first 1000 sequences (4000 lines) from paired-end FASTQ files
```bash
head -n 4000 SRR2584866_1.trim.sub.fastq > 1000_reads_1.fastq

head -n 4000 SRR2584866_2.trim.sub.fastq > 1000_reads_2.fastq
```
