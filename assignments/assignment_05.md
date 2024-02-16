---
title: 'Homework 6: bedtools and GenomicRanges'
author: "Your name here"
date: "2024-02-13"
output:
  html_document:
    df_print: paged
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


# Homework 6

For each of these questions, show the commands you would use to answer the question, If the answer can be shown in less than 5 lines, show the answer as well.

### Question #1
Create a BED file representing all of the intervals in the genome that are NOT exonic.

***
### Question #2

What is the average distance from GWAS SNPs to the closest exon? (Hint - have a look at the `closest` tool.)

***
### Question #3
Count how many exons occur in each 500kb interval (“window”) in the human genome, what is the average value? (Hint - have a look at the `makewindows` tool.)

***
### Question #4
Are there any exons that are completely overlapped by an enhancer? If so, how many?

***
### Question #5
What fraction of the GWAS SNPs are exonic?

***
### Question #6
What fraction of the GWAS SNPs lie in either enhancers or promoters in the hESC data we have?

***
### Question #7
Create intervals representing the canonical 2bp splice sites on either side of each exon (don’t worry about excluding splice sites at the first or last exon). (Hint - have a look at the flank tool.)

***
### Question #8
Which hESC ChromHMM state (e.g., `11_Weak_Txn`, `10_Txn_Elongation`) represents the most number of base pairs in each of chromosome 19 and chromosome 8?

***

### Question #9
Correct the following command so that it completes:
awk '(NR >= 100 && NR <= 200) || ($1 == "chr18)' cpg.bed

***

### Question #10
Correct the following command so that it completes:
egrep --color=always "ATC[ATGC+TGA" Mov10_oe_1.subset.fq | head

***

### Question #11
Correct the following command to copy the genome.txt file into the current directory:
cp ~/data/genome.txt

***

### Question #12
Correct the following command to copy the directory ~/data into the directory ~/lecture6:
cp ~/data ~/lecture6

***
