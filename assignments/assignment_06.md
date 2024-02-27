---
title: 'Homework 6: NGS File Formats'
author: "Your name here"
date: "2024-02-27"
output:
  html_document:
    df_print: paged
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

1. Recall that a FASTQ file contains 4 lines per entry.

For the following entries, what position in the first read is most confidently called? what is it's Quality score (Q) and Probability of incorrect call (P)?

For the following entries, what position in the second read is least confidently called? What is it's Quality score (Q) and Probability of incorrect call (P)?

```
@SRR794330.2 HWI-ST434:134117522:C1N85ACXX:8:1101:1493:2158/1
TGGCTTTGAAGAAGGAGGATGGGGCCACCAGCCAAGGAATGCAGGGAGCCTCTAGAAATTAGAAAAGGCAAGGCAACAGATTCTCCCCTAAAGCCTCCAG
+
;=?ADDD?<4C4;CEFEGE@@7@+@C;EDB)CDDE9B@>?D9??@FA<@@FAE@=D@EEE.)=?BB75=@@1;2(5;<9?>@@B>A:@B?@@########
@SRR794330.3 HWI-ST434:134117522:C1N85ACXX:8:1101:1684:2048/1
ATACAAAAATTAGCTGGGCATGGTGGTGTGCACCTGTAATCCCAGCTACTTGGGAAGCTGAGGCAGGAGAATCGCTTGAACCTGGGAGGTAGAGGTTGCA
+
<@<BBADABHHFFIJIIJG>FHGHGHCFFGIIJJJIAGADFGGIFHD@DDGHEICGA@FFAHGECC>CD?@;>AC@A??AABC2???@@2>@CC?:?CB@
```


2. Let's say you have aligned some RNA-Seq, and now you have a BAM file. You would like to know how many reads you have for each gene in your sample's genome.

What is one tool that can be used to count gene abundance, and what type of file would you need to describe the the position and structure of those genes.


3. With this BAM file, you want see all the reads that you have that are paired, mapped in a proper pair, and are duplicates, what SAM flags will you filter for?

4. Decode the following CIGAR strings:

`76H130M`
`15M8D4I7M1X7M2X5M1X158M`
