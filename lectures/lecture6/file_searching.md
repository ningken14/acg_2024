<style>
.small-code pre code {
  font-size: 1 em;
}
.reveal small {
  position: absolute;
  bottom: 0;
  font-size: 50%;
}
.reveal .slides{
    width: 90% !important;  
    height: 90% !important;
}
.reveal p{
  line-height: 1.5 !important;
}
.reveal code{
  background-color: #F8F8F8;
  color: #1E90FF;
}
.reveal code.bad {
  color: #B73239;
  background-color: #F4F4F4;
}
.reveal pre code {
  white-space: pre;
  overflow: auto;
  color: #00008B;
}
</style>

Lecture 6 - Working with data on the command line (and a bit in R too!)
========================================================
author: Simon Coetzee
date: 02/13/2024
width: 1440
height: 900
transition: fade


## **Genome Arithmetic with bedtools + advanced grep and awk**

<small>This work, "Working with data on the command line (and a bit in R too!)", is a derivative of part of ["Applied Computational Genomics Course at UU"](https://github.com/quinlan-lab/applied-computational-genomics) by quinlan-lab, and ["bedtools Tutorial"](http://quinlanlab.org/tutorials/bedtools/bedtools.html) by quinlan-lab used under [CC BY-SA](https://creativecommons.org/licenses/by-sa/4.0/). "Working with data on the command line (and a bit in R too!)" is licensed under [CC BY](https://creativecommons.org/licenses/by/4.0/) by Simon Coetzee.</small>

Goals (bedtools and GenomicRanges)
========================================================

- Learn to use bedtools to accomplish genome arithmetic.
- Learn similar skills in R with `GenomicRanges`


Genome Arithmetic with bedtools
========================================================
type: sub-section

<center>
![](bedtools.swiss.png)
</center>

bedtools allows one to intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF.

GenomicRanges
========================================================

The `GenomicRanges` package serves as the foundation for representing genomic locations within the Bioconductor project.
This package lays a foundation for genomic analysis by introducing three classes (`GRanges`, `GPos`, and `GRangesList`), which are used to represent genomic ranges, genomic positions, and groups of genomic ranges.

HelloRanges
========================================================
`HelloRanges` is a tool that will help us transfer our bedtools knowledge to R.
We will install it first.

```r
BiocManager::install("HelloRanges")
library(HelloRanges)
```


Genome Arithmetic with bedtools
========================================================

bedtools is a command line tool with many subcommands to perform specific tasks.
```
bedtools is a powerful toolset for genome arithmetic.

Version:   v2.30.0
About:     developed in the quinlanlab.org and by many contributors worldwide.
Docs:      http://bedtools.readthedocs.io/
Code:      https://github.com/arq5x/bedtools2
Mail:      https://groups.google.com/forum/#!forum/bedtools-discuss

Usage:     bedtools <subcommand> [options]

The bedtools sub-commands include:
```
Genome Arithmetic with bedtools
========================================================

bedtools is a command line tool with many subcommands to perform specific tasks.
```
[ Genome arithmetic ]
    intersect     Find overlapping intervals in various ways.
    window        Find overlapping intervals within a window around an interval.
    closest       Find the closest, potentially non-overlapping interval.
    coverage      Compute the coverage over defined intervals.
    map           Apply a function to a column for each overlapping interval.
    genomecov     Compute the coverage over an entire genome.
    merge         Combine overlapping/nearby intervals into a single interval.
    cluster       Cluster (but don't merge) overlapping/nearby intervals.
    complement    Extract intervals _not_ represented by an interval file.
    shift         Adjust the position of intervals.
    subtract      Remove intervals based on overlaps b/w two files.
    slop          Adjust the size of intervals.
    flank         Create new intervals from the flanks of existing intervals.
    sort          Order the intervals in a file.
    random        Generate random intervals in a genome.
    shuffle       Randomly redistribute intervals in a genome.
    sample        Sample random records from file using reservoir sampling.
    spacing       Report the gap lengths between intervals in a file.
    annotate      Annotate coverage of features from multiple files.
```

Genome Arithmetic with bedtools
========================================================

bedtools is a command line tool with many subcommands to perform specific tasks.

each tool is called in the form `bedtools <subcommand>` as follows:
```
$ bedtools intersect

Tool:    bedtools intersect (aka intersectBed)
Version: v2.30.0
Summary: Report overlaps between two feature files.

Usage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>

	Note: -b may be followed with multiple databases and/or
	wildcard (*) character(s).
Options:
	-wa	Write the original entry in A for each overlap.

	-wb	Write the original entry in B for each overlap.
		- Useful for knowing _what_ A overlaps. Restricted by -f and -r.
...		
```

Some example files
========================================================
type: sub-section

- `cpg.bed` - represents CpG islands in the human genome
- `exons.bed` - represents RefSeq exons from human genes
- `gwas.bed` - represents human disease-associated SNPs that were identified in genome-wide association studies
- `hesc.chromHmm.bed` - represents the predicted function (by chromHMM) of each interval in the genome of a human embryonic stem cell based upon ChIP-seq experiments from ENCODE

I'll take a look at these quickly in IGV^1 - an easy to use local genome browser

<small>
1. https://software.broadinstitute.org/software/igv/
</small>

bedtools - intersect
========================================================

This is probably the command you will use the most often.

It compares 2 or more bed files (or bam, or vcf, or gtf) and identifies all the
regions in the genome where the features overlap.

Want to know where one feature lies relative to another, this is your tool.

bedtools - intersect
========================================================

It compares 2 or more bed files (or bam, or vcf, or gtf) and identifies all the
regions in the genome where the features overlap.

<center>
![](intersect.png)
</center>

bedtools - intersect
========================================================
### Default behavior
By default, intersect reports the intervals that represent overlaps between your
two files. To demonstrate, let’s identify all of the CpG islands that overlap exons.

**these are not the original cpg intervals, *only the portion that overlaps* with the exons**
```
$ cp ~/data/*.bed ~/lecture6
$ cd ~/lecture6
$ bedtools intersect -a cpg.bed -b exons.bed | head -5
chr1    29320   29370   CpG:_116
chr1    135124  135563  CpG:_30
chr1    327790  328229  CpG:_29
chr1    327790  328229  CpG:_29
chr1    327790  328229  CpG:_29
$ █
```

GRanges - intersect
========================================================
### Default behavior
By default, intersect reports the intervals that represent overlaps between your
two files. To demonstrate, let’s identify all of the CpG islands that overlap exons.

**these are not the original cpg intervals, *only the portion that overlaps* with the exons**

```r
setwd("~/lecture6")
code <- bedtools_intersect("-a cpg.bed -b exons.bed")
code
```

```
{
    genome <- Seqinfo(genome = NA_character_)
    gr_a <- import("cpg.bed", genome = genome)
    gr_b <- import("exons.bed", genome = genome)
    pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    ans <- pintersect(pairs, ignore.strand = TRUE)
    ans
}
```

GRanges - intersect
========================================================
### Default behavior
By default, intersect reports the intervals that represent overlaps between your
two files. To demonstrate, let’s identify all of the CpG islands that overlap exons.

**these are not the original cpg intervals, *only the portion that overlaps* with the exons**

```r
eval(code)
```

```
GRanges object with 45500 ranges and 2 metadata columns:
          seqnames            ranges strand |        name       hit
             <Rle>         <IRanges>  <Rle> | <character> <logical>
      [1]     chr1       29321-29370      * |    CpG:_116      TRUE
      [2]     chr1     135125-135563      * |     CpG:_30      TRUE
      [3]     chr1     327791-328229      * |     CpG:_29      TRUE
      [4]     chr1     327791-328229      * |     CpG:_29      TRUE
      [5]     chr1     327791-328229      * |     CpG:_29      TRUE
      ...      ...               ...    ... .         ...       ...
  [45496]     chrY 59213949-59214117      * |     CpG:_36      TRUE
  [45497]     chrY 59213949-59214117      * |     CpG:_36      TRUE
  [45498]     chrY 59213949-59214117      * |     CpG:_36      TRUE
  [45499]     chrY 59213949-59214117      * |     CpG:_36      TRUE
  [45500]     chrY 59213949-59214117      * |     CpG:_36      TRUE
  -------
  seqinfo: 69 sequences from an unspecified genome; no seqlengths
```

bedtools - intersect
========================================================
### Report the original interval in each file.

The `-wa` (write A) and `-wb` (write B) options allow one to see the original 
records from the A and B files that overlapped. As such, instead of not only 
showing you where the intersections occurred, it shows you what intersected.

```
bedtools intersect -a cpg.bed -b exons.bed -wa -wb | head -5
chr1    28735   29810   CpG:_116    chr1    29320   29370   NR_024540_exon_10_0_chr1_29321_r        -
chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_039983_exon_0_0_chr1_134773_r    0   -
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028322_exon_2_0_chr1_324439_f    0   +
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028325_exon_2_0_chr1_324439_f    0   +
chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_028327_exon_3_0_chr1_327036_f    0   +
$ █
```

GRanges - intersect
========================================================
### Report the original interval in each file.

The `-wa` (write A) and `-wb` (write B) options allow one to see the original 
records from the A and B files that overlapped. As such, instead of not only 
showing you where the intersections occurred, it shows you what intersected.


```r
code <- bedtools_intersect("-a cpg.bed -b exons.bed -wa -wb")
code
```

```
{
    genome <- Seqinfo(genome = NA_character_)
    gr_a <- import("cpg.bed", genome = genome)
    gr_b <- import("exons.bed", genome = genome)
    pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    ans <- pairs
    ans
}
```

GRanges - intersect
========================================================
### Report the original interval in each file.

The `-wa` (write A) and `-wb` (write B) options allow one to see the original 
records from the A and B files that overlapped. As such, instead of not only 
showing you where the intersections occurred, it shows you what intersected.


```r
eval(code)
```

```
Pairs object with 45500 pairs and 0 metadata columns:
                           first                   second
                       <GRanges>                <GRanges>
      [1]       chr1:28736-29810       chr1:29321-29370:-
      [2]     chr1:135125-135563     chr1:134773-139696:-
      [3]     chr1:327791-328229     chr1:324439-328581:+
      [4]     chr1:327791-328229     chr1:324439-328581:+
      [5]     chr1:327791-328229     chr1:327036-328581:+
      ...                    ...                      ...
  [45496] chrY:59213795-59214183 chrY:59213949-59214117:+
  [45497] chrY:59213795-59214183 chrY:59213949-59214117:+
  [45498] chrY:59213795-59214183 chrY:59213949-59214117:+
  [45499] chrY:59213795-59214183 chrY:59213949-59214117:+
  [45500] chrY:59213795-59214183 chrY:59213949-59214117:+
```

bedtools - intersect
========================================================
### How many base pairs of overlap were there?

The `-wo` (write overlap) option allows one to also report the number of base 
pairs of overlap between the features that overlap between each of the files.

```
$ bedtools intersect -a cpg.bed -b exons.bed -wo | head -5
chr1    28735   29810   CpG:_116    chr1    29320   29370   NR_024540_exon_10_0_chr1_29321_r        -   50
chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_039983_exon_0_0_chr1_134773_r    0       439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028322_exon_2_0_chr1_324439_f    0       439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028325_exon_2_0_chr1_324439_f    0       439
chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_028327_exon_3_0_chr1_327036_f    0       439
$ █
```

GRanges - intersect
========================================================
### How many base pairs of overlap were there?

The `-wo` (write overlap) option allows one to also report the number of base 
pairs of overlap between the features that overlap between each of the files.


```r
code <- bedtools_intersect("-a cpg.bed -b exons.bed -wo")
code
```

```
{
    genome <- Seqinfo(genome = NA_character_)
    gr_a <- import("cpg.bed", genome = genome)
    gr_b <- import("exons.bed", genome = genome)
    pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    ans <- pairs
    mcols(ans)$overlap_width <- width(pintersect(ans, ignore.strand = TRUE))
    ans
}
```

GRanges - intersect
========================================================
### How many base pairs of overlap were there?

The `-wo` (write overlap) option allows one to also report the number of base 
pairs of overlap between the features that overlap between each of the files.


```r
eval(code)
```

```
Pairs object with 45500 pairs and 1 metadata column:
                           first                   second | overlap_width
                       <GRanges>                <GRanges> |     <integer>
      [1]       chr1:28736-29810       chr1:29321-29370:- |            50
      [2]     chr1:135125-135563     chr1:134773-139696:- |           439
      [3]     chr1:327791-328229     chr1:324439-328581:+ |           439
      [4]     chr1:327791-328229     chr1:324439-328581:+ |           439
      [5]     chr1:327791-328229     chr1:327036-328581:+ |           439
      ...                    ...                      ... .           ...
  [45496] chrY:59213795-59214183 chrY:59213949-59214117:+ |           169
  [45497] chrY:59213795-59214183 chrY:59213949-59214117:+ |           169
  [45498] chrY:59213795-59214183 chrY:59213949-59214117:+ |           169
  [45499] chrY:59213795-59214183 chrY:59213949-59214117:+ |           169
  [45500] chrY:59213795-59214183 chrY:59213949-59214117:+ |           169
```

bedtools - intersect
========================================================
### Counting the number of overlapping features.

We can also count, for each feature in the “A” file, the number of overlapping 
features in the “B” file. This is handled with the `-c` option.

```
$ bedtools intersect -a cpg.bed -b exons.bed -c | head -5
chr1    28735   29810   CpG:_116    1
chr1    135124  135563  CpG:_30 1
chr1    327790  328229  CpG:_29 3
chr1    437151  438164  CpG:_84 0
chr1    449273  450544  CpG:_99 0
$ █
```

GRanges - intersect
========================================================
### Counting the number of overlapping features.

We can also count, for each feature in the “A” file, the number of overlapping 
features in the “B” file. This is handled with the `-c` option.


```r
code <- bedtools_intersect("-a cpg.bed -b exons.bed -c")
code
```

```
{
    genome <- Seqinfo(genome = NA_character_)
    gr_a <- import("cpg.bed", genome = genome)
    gr_b <- import("exons.bed", genome = genome)
    ans <- gr_a
    mcols(ans)$overlap_count <- countOverlaps(gr_a, gr_b, ignore.strand = TRUE)
    ans
}
```

GRanges - intersect
========================================================
### Counting the number of overlapping features.

We can also count, for each feature in the “A” file, the number of overlapping 
features in the “B” file. This is handled with the `-c` option.


```r
eval(code)
```

```
GRanges object with 28691 ranges and 2 metadata columns:
          seqnames            ranges strand |        name overlap_count
             <Rle>         <IRanges>  <Rle> | <character>     <integer>
      [1]     chr1       28736-29810      * |    CpG:_116             1
      [2]     chr1     135125-135563      * |     CpG:_30             1
      [3]     chr1     327791-328229      * |     CpG:_29             3
      [4]     chr1     437152-438164      * |     CpG:_84             0
      [5]     chr1     449274-450544      * |     CpG:_99             0
      ...      ...               ...    ... .         ...           ...
  [28687]     chrY 27610116-27611088      * |     CpG:_76             0
  [28688]     chrY 28555536-28555932      * |     CpG:_32             0
  [28689]     chrY 28773316-28773544      * |     CpG:_25             0
  [28690]     chrY 59213795-59214183      * |     CpG:_36             5
  [28691]     chrY 59349267-59349574      * |     CpG:_29             0
  -------
  seqinfo: 69 sequences from an unspecified genome; no seqlengths
```


bedtools - intersect
========================================================
### Find features that DO NOT overlap
remember `grep -v` ?!

Often we want to identify those features in our A file that do not overlap features in the B file.

```
$ bedtools intersect -a cpg.bed -b exons.bed -v | head -5
chr1    437151  438164  CpG:_84
chr1    449273  450544  CpG:_99
chr1    533219  534114  CpG:_94
chr1    544738  546649  CpG:_171
chr1    801975  802338  CpG:_24
$ █
```

GRanges - intersect
========================================================
### Find features that DO NOT overlap
remember `grep -v` ?!

Often we want to identify those features in our A file that do not overlap features in the B file.


```r
code <- bedtools_intersect("-a cpg.bed -b exons.bed -v")
code
```

```
{
    genome <- Seqinfo(genome = NA_character_)
    gr_a <- import("cpg.bed", genome = genome)
    gr_b <- import("exons.bed", genome = genome)
    subsetByOverlaps(gr_a, gr_b, invert = TRUE, ignore.strand = TRUE)
}
```

GRanges - intersect
========================================================
### Find features that DO NOT overlap
remember `grep -v` ?!

Often we want to identify those features in our A file that do not overlap features in the B file.


```r
eval(code)
```

```
GRanges object with 9846 ranges and 1 metadata column:
         seqnames            ranges strand |        name
            <Rle>         <IRanges>  <Rle> | <character>
     [1]     chr1     437152-438164      * |     CpG:_84
     [2]     chr1     449274-450544      * |     CpG:_99
     [3]     chr1     533220-534114      * |     CpG:_94
     [4]     chr1     544739-546649      * |    CpG:_171
     [5]     chr1     801976-802338      * |     CpG:_24
     ...      ...               ...    ... .         ...
  [9842]     chrY 26351344-26352316      * |     CpG:_76
  [9843]     chrY 27610116-27611088      * |     CpG:_76
  [9844]     chrY 28555536-28555932      * |     CpG:_32
  [9845]     chrY 28773316-28773544      * |     CpG:_25
  [9846]     chrY 59349267-59349574      * |     CpG:_29
  -------
  seqinfo: 69 sequences from an unspecified genome; no seqlengths
```

bedtools - intersect
========================================================
### Require a minimal fraction of overlap.

Recall that the default is to report overlaps between features in A and B so 
long as at least one base pair of overlap exists. However, the `-f` option 
allows you to specify what fraction of each feature in A should be overlapped by
a feature in B before it is reported.

Let’s be more strict and require 50% of overlap.

```
$ bedtools intersect -a cpg.bed -b exons.bed -wo -f 0.50 | head -5
chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_039983_exon_0_0_chr1_134773_r    0       439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028322_exon_2_0_chr1_324439_f    0       439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028325_exon_2_0_chr1_324439_f    0       439
chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_028327_exon_3_0_chr1_327036_f    0       439
chr1    788863  789211  CpG:_28 chr1    788770  794826  NR_047519_exon_5_0_chr1_788771_f    0       348
$ █
```

GRanges - intersect
========================================================
### Require a minimal fraction of overlap.

Recall that the default is to report overlaps between features in A and B so 
long as at least one base pair of overlap exists. However, the `-f` option 
allows you to specify what fraction of each feature in A should be overlapped by
a feature in B before it is reported.

Let’s be more strict and require 50% of overlap.


```r
code <- bedtools_intersect("-a cpg.bed -b exons.bed -wo -f 0.50")
code
```

```
{
    genome <- Seqinfo(genome = NA_character_)
    gr_a <- import("cpg.bed", genome = genome)
    gr_b <- import("exons.bed", genome = genome)
    pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    olap <- pintersect(pairs, ignore.strand = TRUE)
    keep <- width(olap)/width(first(pairs)) >= 0.5
    pairs <- pairs[keep]
    ans <- pairs
    mcols(ans)$overlap_width <- width(olap)[keep]
    ans
}
```

GRanges - intersect
========================================================
### Require a minimal fraction of overlap.

Recall that the default is to report overlaps between features in A and B so 
long as at least one base pair of overlap exists. However, the `-f` option 
allows you to specify what fraction of each feature in A should be overlapped by
a feature in B before it is reported.

Let’s be more strict and require 50% of overlap.


```r
eval(code)
```

```
Pairs object with 12669 pairs and 1 metadata column:
                           first                   second | overlap_width
                       <GRanges>                <GRanges> |     <integer>
      [1]     chr1:135125-135563     chr1:134773-139696:- |           439
      [2]     chr1:327791-328229     chr1:324439-328581:+ |           439
      [3]     chr1:327791-328229     chr1:324439-328581:+ |           439
      [4]     chr1:327791-328229     chr1:327036-328581:+ |           439
      [5]     chr1:788864-789211     chr1:788771-794826:+ |           348
      ...                    ...                      ... .           ...
  [12665] chrY:26979890-26980116 chrY:26979967-26980276:+ |           150
  [12666] chrY:26979890-26980116 chrY:26979967-26980276:+ |           150
  [12667] chrY:26979890-26980116 chrY:26979967-26980276:+ |           150
  [12668] chrY:26979890-26980116 chrY:26979990-26980276:+ |           127
  [12669] chrY:26979890-26980116 chrY:26979990-26980276:+ |           127
```

bedtools - intersect
========================================================
### Bedtools is faster when your files are sorted.

We can sort with the tools we already know (these files are all pre-sorted already)
```
$ sort -k1,1V -k2,2n gwas.bed > gwas.sorted.bed
$ █
```
<center>
![](speed-compare.png)
</center>

bedtools - intersect
========================================================
### Intersecting multiple files at once
bedtools is able to intersect an “A” file against one or more “B” files. This 
greatly simplifies analyses involving multiple datasets relevant to a given
experiment. For example, let’s intersect exons with CpG islands, GWAS SNPs, 
and the ChromHMM annotations.

```
$ bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted | head -5
chr1    11873   11937   NR_046018_exon_0_0_chr1_11874_f 0   +
chr1    11937   12137   NR_046018_exon_0_0_chr1_11874_f 0   +
chr1    12137   12227   NR_046018_exon_0_0_chr1_11874_f 0   +
chr1    12612   12721   NR_046018_exon_1_0_chr1_12613_f 0   +
chr1    13220   14137   NR_046018_exon_2_0_chr1_13221_f 0   +
$ █
```

bedtools - intersect
========================================================
Now by default, this isn’t incredibly informative as we can’t tell which of 
the three **“B”** files yielded the intersection with each exon. However, if we use 
the `-wa` and `-wb` options, we can see from which file number (following the order 
of the files given on the command line) the intersection came. In this case, 
the 7th column reflects this file number.

```
bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted -wa -wb \
  | head -10000 \
  | tail -5
chr1    27632676    27635124    NM_015023_exon_15_0_chr1_27632677_f 0   +   3   chr1    27635013    27635413    7_Weak_Enhancer
chr1    27648635    27648882    NM_032125_exon_0_0_chr1_27648636_f  0   +   1   chr1    27648453    27649006    CpG:_63
chr1    27648635    27648882    NM_032125_exon_0_0_chr1_27648636_f  0   +   3   chr1    27648613    27649413    1_Active_Promoter
chr1    27648635    27648882    NR_037576_exon_0_0_chr1_27648636_f  0   +   1   chr1    27648453    27649006    CpG:_63
chr1    27648635    27648882    NR_037576_exon_0_0_chr1_27648636_f  0   +   3   chr1    27648613    27649413    1_Active_Promoter
$ █
```

bedtools - intersect
========================================================
Additionally, one can use file “labels” instead of file numbers to facilitate 
interpretation, especially when there are many files involved.

```
$ bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted -wa -wb -names cpg gwas chromhmm \
  | head -10000 \
  | tail -5
chr1    27632676    27635124    NM_015023_exon_15_0_chr1_27632677_f 0   +   chromhmm    chr1    27635013    27635413    7_Weak_Enhancer
chr1    27648635    27648882    NM_032125_exon_0_0_chr1_27648636_f  0   +   cpg chr1    27648453    27649006    CpG:_63
chr1    27648635    27648882    NM_032125_exon_0_0_chr1_27648636_f  0   +   chromhmm    chr1    27648613    27649413    1_Active_Promoter
chr1    27648635    27648882    NR_037576_exon_0_0_chr1_27648636_f  0   +   cpg chr1    27648453    27649006    CpG:_63
chr1    27648635    27648882    NR_037576_exon_0_0_chr1_27648636_f  0   +   chromhmm    chr1    27648613    27649413    1_Active_Promoter
$ █
```

GRanges - intersect
========================================================
Additionally, one can use file “labels” instead of file numbers to facilitate 
interpretation, especially when there are many files involved.


```r
code <- {
  genome <- Seqinfo(genome = NA_character_)
  gr_a <- import("exons.bed", genome = genome)
  b <- c("cpg.bed", "gwas.bed", "hesc.chromHmm.bed")
  names(b) <- c("cpg", "gwas", "chromhmm")
  bl <- List(lapply(b, import, genome = genome))
  gr_b <- stack(bl, "b")
  pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
  ans <- pairs
}
eval(code)
```

```
Pairs object with 541210 pairs and 0 metadata columns:
                              first                 second
                          <GRanges>              <GRanges>
       [1]       chr1:11874-12227:+       chr1:11538-11937
       [2]       chr1:11874-12227:+       chr1:11938-12137
       [3]       chr1:11874-12227:+       chr1:12138-14137
       [4]       chr1:12613-12721:+       chr1:12138-14137
       [5]       chr1:13221-14409:+       chr1:12138-14137
       ...                      ...                    ...
  [541206] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541207] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541208] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541209] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541210] chrY:59213949-59214117:+ chrY:59213795-59214183
```

GRanges - intersect
========================================================
Additionally, one can use file “labels” instead of file numbers to facilitate 
interpretation, especially when there are many files involved.


```r
eval(code)
```

```
Pairs object with 541210 pairs and 0 metadata columns:
                              first                 second
                          <GRanges>              <GRanges>
       [1]       chr1:11874-12227:+       chr1:11538-11937
       [2]       chr1:11874-12227:+       chr1:11938-12137
       [3]       chr1:11874-12227:+       chr1:12138-14137
       [4]       chr1:12613-12721:+       chr1:12138-14137
       [5]       chr1:13221-14409:+       chr1:12138-14137
       ...                      ...                    ...
  [541206] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541207] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541208] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541209] chrY:59213949-59214117:+ chrY:59213795-59214183
  [541210] chrY:59213949-59214117:+ chrY:59213795-59214183
```

GRanges - intersect
========================================================
Additionally, one can use file “labels” instead of file numbers to facilitate 
interpretation, especially when there are many files involved.


```r
head(first(eval(code)))
```

```
GRanges object with 6 ranges and 2 metadata columns:
      seqnames      ranges strand |                   name     score
         <Rle>   <IRanges>  <Rle> |            <character> <numeric>
  [1]     chr1 11874-12227      + | NR_046018_exon_0_0_c..         0
  [2]     chr1 11874-12227      + | NR_046018_exon_0_0_c..         0
  [3]     chr1 11874-12227      + | NR_046018_exon_0_0_c..         0
  [4]     chr1 12613-12721      + | NR_046018_exon_1_0_c..         0
  [5]     chr1 13221-14409      + | NR_046018_exon_2_0_c..         0
  [6]     chr1 13221-14409      + | NR_046018_exon_2_0_c..         0
  -------
  seqinfo: 49 sequences from an unspecified genome; no seqlengths
```

```r
head(second(eval(code)))
```

```
GRanges object with 6 ranges and 2 metadata columns:
      seqnames      ranges strand |        b              name
         <Rle>   <IRanges>  <Rle> |    <Rle>       <character>
  [1]     chr1 11538-11937      * | chromhmm       11_Weak_Txn
  [2]     chr1 11938-12137      * | chromhmm 14_Repetitive/CNV
  [3]     chr1 12138-14137      * | chromhmm       11_Weak_Txn
  [4]     chr1 12138-14137      * | chromhmm       11_Weak_Txn
  [5]     chr1 12138-14137      * | chromhmm       11_Weak_Txn
  [6]     chr1 14138-27537      * | chromhmm  9_Txn_Transition
  -------
  seqinfo: 70 sequences from an unspecified genome; no seqlengths
```

bedtools - merge
========================================================
Many datasets of genomic features have many individual features that overlap one 
another (e.g. aligments from a ChiP seq experiment). It is often useful to just 
combine the overlapping into a single, contiguous interval.

<center>
![](merge.png)
</center>

bedtools - merge
========================================================
### Input must be sorted

The merge tool requires that the input file is sorted by chromosome, then by start position. This allows the merging algorithm to work very quickly without requiring any RAM.

If your files are unsorted, the merge tool will raise an error.

We can use bedtools sort to sort by a specific chromosome order - for example, exactly the one your colleague who aligned your data used with their reference genome.

```
$ cp ~/data/genome.txt ~/lecture6
$ bedtools sort -g genome.txt -i exons.bed > exons.sorted.bed
$ █
```

bedtools - merge
========================================================
### Merge intervals.
Merging results in a new set of intervals representing the merged set of 
intervals in the input. That is, if a base pair in the genome is covered by 
ten features, it will now only be represented once in the output file.

```
$ bedtools merge -i exons.bed | head -n 10
chr1    11873   12227
chr1    12612   12721
chr1    13220   14829
chr1    14969   15038
chr1    15795   15947
chr1    16606   16765
chr1    16857   17055
chr1    17232   17368
chr1    17605   17742
chr1    17914   18061
$ █
```

bedtools - merge
========================================================
### Count the number of overlapping intervals.
A more sophisticated approach would be to not only merge overlapping intervals,
but also report the number of intervals that were integrated into the new, 
merged interval.

The `-c` option allows one to specify a column or columns in the input that you wish to summarize

The `-o` option defines the operation(s) that you wish to apply to each column listed for the -c option
```
$ bedtools merge -i exons.bed -c 1 -o count | head -n 10
chr1    11873   12227   1
chr1    12612   12721   1
chr1    13220   14829   2
chr1    14969   15038   1
chr1    15795   15947   1
chr1    16606   16765   1
chr1    16857   17055   1
chr1    17232   17368   1
chr1    17605   17742   1
chr1    17914   18061   1
$ █
```

bedtools - merge
========================================================
### Merge operations
```
-o	Specify the operation that should be applied to -c.
		Valid operations:
		    sum, min, max, absmin, absmax,
		    mean, median, mode, antimode
		    stdev, sstdev
		    collapse (i.e., print a delimited list (duplicates allowed)),
		    distinct (i.e., print a delimited list (NO duplicates allowed)),
		    distinct_sort_num (as distinct, sorted numerically, ascending),
		    distinct_sort_num_desc (as distinct, sorted numerically, desscending),
		    distinct_only (delimited list of only unique values),
		    count
		    count_distinct (i.e., a count of the unique values in the column),
		    first (i.e., just the first value in the column),
		    last (i.e., just the last value in the column),
		Default: sum
		Multiple operations can be specified in a comma-delimited list.
```

bedtools - merge
========================================================
### Merging features that are close to one another.

With the `-d` (distance) option, one can also merge intervals that do not
overlap, yet are close to one another. For example, to merge features that are
no more than 1000bp apart, one would run:

```
$ bedtools merge -i exons.bed -d 1000 -c 1 -o count | head -10
chr1    11873   18366   12
chr1    24737   24891   1
chr1    29320   29370   1
chr1    34610   36081   6
chr1    69090   70008   1
chr1    134772  140566  3
chr1    323891  328581  10
chr1    367658  368597  3
chr1    621095  622034  3
chr1    661138  665731  3
$ █
```

bedtools - merge
========================================================
### Listing the name of each of the exons that were merged.

Many times you want to keep track of the details of exactly which intervals were merged. One way to do this is to create a list of the names of each feature. We can do with with the **collapse** operation available via the `-o` argument.

```
$ bedtools merge -i exons.bed -d 90 -c 1,4 -o count,collapse | head -20
chr1    11873   12227   1   NR_046018_exon_0_0_chr1_11874_f
chr1    12612   12721   1   NR_046018_exon_1_0_chr1_12613_f
chr1    13220   14829   2   NR_046018_exon_2_0_chr1_13221_f,NR_024540_exon_0_0_chr1_14362_r
chr1    14969   15038   1   NR_024540_exon_1_0_chr1_14970_r
chr1    15795   15947   1   NR_024540_exon_2_0_chr1_15796_r
chr1    16606   16765   1   NR_024540_exon_3_0_chr1_16607_r
chr1    16857   17055   1   NR_024540_exon_4_0_chr1_16858_r
chr1    17232   17368   1   NR_024540_exon_5_0_chr1_17233_r
chr1    17605   17742   1   NR_024540_exon_6_0_chr1_17606_r
chr1    17914   18061   1   NR_024540_exon_7_0_chr1_17915_r
$ █
```

GRanges - merge
========================================================
### Listing the name of each of the exons that were merged.

Many times you want to keep track of the details of exactly which intervals were merged. One way to do this is to create a list of the names of each feature. We can do with with the **collapse** operation available via the `-o` argument.


```r
code <- bedtools_merge("-i exons.bed -d 90 -c 1,4 -o count,collapse")
code
```

```
{
    genome <- Seqinfo(genome = NA_character_)
    gr_a <- import("exons.bed", genome = genome)
    ans <- reduce(gr_a, ignore.strand = TRUE, with.revmap = TRUE, 
        min.gapwidth = 91L)
    mcols(ans) <- aggregate(gr_a, mcols(ans)$revmap, seqnames.count = lengths(seqnames), 
        name.collapse = unstrsplit(name, ","), drop = FALSE)
    ans
}
```

GRanges - merge
========================================================
### Listing the name of each of the exons that were merged.

Many times you want to keep track of the details of exactly which intervals were merged. One way to do this is to create a list of the names of each feature. We can do with with the **collapse** operation available via the `-o` argument.


```r
eval(code)
```

```
GRanges object with 221100 ranges and 3 metadata columns:
           seqnames            ranges strand |             grouping
              <Rle>         <IRanges>  <Rle> | <ManyToManyGrouping>
       [1]     chr1       11874-12227      * |                    1
       [2]     chr1       12613-12721      * |                    2
       [3]     chr1       13221-14829      * |                  3,4
       [4]     chr1       14970-15038      * |                    5
       [5]     chr1       15796-15947      * |                    6
       ...      ...               ...    ... .                  ...
  [221096]     chrY 59337091-59337236      * |        459744,459745
  [221097]     chrY 59337949-59338150      * |        459746,459747
  [221098]     chrY 59338754-59338859      * |        459748,459749
  [221099]     chrY 59340194-59340278      * |               459750
  [221100]     chrY 59342487-59343488      * |        459751,459752
           seqnames.count          name.collapse
                <integer>            <character>
       [1]              1 NR_046018_exon_0_0_c..
       [2]              1 NR_046018_exon_1_0_c..
       [3]              2 NR_046018_exon_2_0_c..
       [4]              1 NR_024540_exon_1_0_c..
       [5]              1 NR_024540_exon_2_0_c..
       ...            ...                    ...
  [221096]              2 NM_002186_exon_4_0_c..
  [221097]              2 NM_002186_exon_5_0_c..
  [221098]              2 NM_002186_exon_6_0_c..
  [221099]              1 NM_002186_exon_7_0_c..
  [221100]              2 NM_002186_exon_8_0_c..
  -------
  seqinfo: 49 sequences from an unspecified genome; no seqlengths
```


bedtools - complement
========================================================
We often want to know which intervals of the genome are NOT “covered” by intervals in a given feature file. For example, if you have a set of ChIP-seq peaks, you may also want to know which regions of the genome are not bound by the factor you assayed.

<center>
![](complement.png)
</center>

bedtools - complement
========================================================
As an example, let’s find all of the non-exonic (i.e., intronic or intergenic) regions of the genome.

Here we see another use for the `genome.txt` file.
```
$ bedtools complement -i exons.bed -g genome.txt > non-exonic.bed
$ head non-exonic.bed
chr1    0   11873
chr1    12227   12612
chr1    12721   13220
chr1    14829   14969
chr1    15038   15795
chr1    15947   16606
chr1    16765   16857
chr1    17055   17232
chr1    17368   17605
chr1    17742   17914
$ █
```

GRanges - complement
========================================================
As an example, let’s find all of the non-exonic (i.e., intronic or intergenic) regions of the genome.

Here we see another use for the `genome.txt` file. To use with the `GRanges` tools, we have to rename it.
```
$ cp ~/lecture6/genome.txt ~/lecture6/hg38.genome
```

```r
code <- bedtools_complement("-i exons.bed -g hg38.genome")
code
```

```
{
    genome <- import("hg38.genome")
    gr_a <- import("exons.bed", genome = genome)
    ans <- setdiff(as(seqinfo(gr_a), "GRanges"), unstrand(gr_a))
    ans
}
```

GRanges - complement
========================================================
As an example, let’s find all of the non-exonic (i.e., intronic or intergenic) regions of the genome.

Here we see another use for the `genome.txt` file. To use with the `GRanges` tools, we have to rename it.

```r
eval(code)
```

```
GRanges object with 229334 ranges and 0 metadata columns:
           seqnames            ranges strand
              <Rle>         <IRanges>  <Rle>
       [1]     chr1           1-11873      *
       [2]     chr1       12228-12612      *
       [3]     chr1       12722-13220      *
       [4]     chr1       14830-14969      *
       [5]     chr1       15039-15795      *
       ...      ...               ...    ...
  [229330]     chrY 59337237-59337948      *
  [229331]     chrY 59338151-59338753      *
  [229332]     chrY 59338860-59340193      *
  [229333]     chrY 59340279-59342486      *
  [229334]     chrY 59343489-59373566      *
  -------
  seqinfo: 93 sequences from hg38 genome
```

bedtools - genomecov
========================================================
For many analyses, one wants to measure the genome wide coverage of a feature file. For example, we often want to know what fraction of the genome is covered by 1 feature, 2 features, 3 features, etc. This is frequently crucial when assessing the “uniformity” of coverage from whole-genome sequencing.

<center>
![](genomecov.png)
</center>

bedtools - genomecov
========================================================
As an example, let’s produce a histogram of coverage of the exons throughout the genome.

* chromosome (or entire genome)
* depth of coverage from features in input file
* number of bases on chromosome (or genome) with depth equal to column 2.
* size of chromosome (or entire genome) in base pairs
* fraction of bases on chromosome (or entire genome) with depth equal to column 2.

```
$ bedtools genomecov -i exons.bed -g genome.txt > exon.histogram
$ less exon.histogram
chr1	0	241996316	249250621	0.970896
chr1	1	4276763	249250621	0.0171585
chr1	2	1475526	249250621	0.00591985
chr1	3	710135	249250621	0.00284908
chr1	4	388193	249250621	0.00155744
chr1	5	179651	249250621	0.000720765
chr1	6	91005	249250621	0.000365114
...
```

bedtools - genomecov
========================================================
### Producing BEDGRAPH output

bedgraph^1 is a format we didn't talk about, it allows display of continuous-valued data in track format. This display type is useful for probability scores, transcriptome data, and coverage! Using the `-bg` option, one can also produce this output which represents the “depth” of feature coverage for each base pair in the genome.

```
$ bedtools genomecov -i exons.bed -g genome.txt -bg | head -10
chr1    11873   12227   1
chr1    12612   12721   1
chr1    13220   14361   1
chr1    14361   14409   2
chr1    14409   14829   1
chr1    14969   15038   1
chr1    15795   15947   1
chr1    16606   16765   1
chr1    16857   17055   1
chr1    17232   17368   1
$ █
```

<small>
1. https://genome.ucsc.edu/goldenPath/help/bedgraph.html
</small>


bedtools - genomecov
========================================================
### Producing BEDGRAPH output

bedgraph^1 is a format we didn't talk about, it allows display of continuous-valued data in track format. This display type is useful for probability scores, transcriptome data, and coverage! Using the `-bg` option, one can also produce this output which represents the “depth” of feature coverage for each base pair in the genome.


```r
code <- bedtools_genomecov("-i exons.bed -g hg38.genome -bg")
code
```

```
{
    genome <- import("hg38.genome")
    gr_a <- import("exons.bed", genome = genome)
    cov <- coverage(gr_a)
    ans <- GRanges(cov)
    ans <- subset(ans, score > 0)
    ans
}
```

bedtools - genomecov
========================================================
### Producing BEDGRAPH output

bedgraph^1 is a format we didn't talk about, it allows display of continuous-valued data in track format. This display type is useful for probability scores, transcriptome data, and coverage! Using the `-bg` option, one can also produce this output which represents the “depth” of feature coverage for each base pair in the genome.


```r
eval(code)
```

```
GRanges object with 245860 ranges and 1 metadata column:
           seqnames            ranges strand |     score
              <Rle>         <IRanges>  <Rle> | <integer>
       [1]     chr1       11874-12227      * |         1
       [2]     chr1       12613-12721      * |         1
       [3]     chr1       13221-14361      * |         1
       [4]     chr1       14362-14409      * |         2
       [5]     chr1       14410-14829      * |         1
       ...      ...               ...    ... .       ...
  [245856]     chrY 59337120-59337236      * |         2
  [245857]     chrY 59337949-59338150      * |         2
  [245858]     chrY 59338754-59338859      * |         2
  [245859]     chrY 59340194-59340278      * |         1
  [245860]     chrY 59342487-59343488      * |         2
  -------
  seqinfo: 93 sequences from an unspecified genome
```

<small>
1. https://genome.ucsc.edu/goldenPath/help/bedgraph.html
</small>

bedtools - a practical example:
========================================================

Let’s imagine you have a BED file of ChiP-seq peaks from two different experiments. You want to identify peaks that were observed in both experiments (requiring 50% reciprocal overlap) and for those peaks, you want to find to find the closest, non-overlapping gene. Such an analysis could be conducted with two, relatively simple bedtools commands.

intersect the peaks from both experiments.
`-f` 0.50 combined with `-r` requires 50% reciprocal overlap between the
peaks from each experiment.
```
$ bedtools intersect -a exp1.bed -b exp2.bed -f 0.50 -r > both.bed
```

find the closest, non-overlapping gene for each interval where
both experiments had a peak
`-io` ignores overlapping intervals and returns only the closest,
non-overlapping interval (in this case, genes)
```
$ bedtools closest -a both.bed -b genes.bed -io > both.nearest.genes.txt
```

Advanced GREP
========================================================
## Just so you know it's there - regular expressions can be even more powerful:
There are certain meta-characters that are reserved in regex:
```
^: matches pattern at start of string
$: matches pattern at end of string
.: matches any character except new lines
[]: matches any of enclose characters
[^]: matches any characters *except* ones enclosed (note: is different from ^)
\: "escapes" meta-characters, allows literal matching
```

Advanced GREP
========================================================
One can do some quite sophisticated searching with grep using regular expressions. 
For these kinds of searches one can use "extended" grep or `egrep`.
I'll use screenshots here to show what is being matched with color:

One can match a certain number of characters with `{x,y}`. Here we are looking 
for matches of any string of C’s and/or G’s that is six to twelve characters 
long using the pattern `[GC]{6,12}`.

<center>
![](grepGC.png)
</center>

Advanced GREP
========================================================
## `+` The most prominent repetition modifier
Instead of asking for a specific length of match, one can ask for at least one
match using the special symbol `+`

As an example we can search for a start codon `ATG` followed by any number of 
bases (at least one), and ending with a stop codon `TGA` using the pattern `ATG[ATGC]+TGA`

<center>
![](grepPlus.png)
</center>

Advanced GREP
========================================================
We can also search for repetitive strings, where we want to search for the whole group of characters, not just one or two.

Here we search for `CAG` repeats between six and twelve repeats long with the pattern `(CAG){6,12}`. You'll note that we used parentheses instead of square brackets here.

<center>
![](grepCAG.png)
</center>

AWK
========================================================

AWK is a programming language that is specifically designed for quickly manipulating space delimited data.
The name AWK comes from the initials of the inventors of the languguage at ATT Bell Labs: Alfred ***A***ho, Peter ***W***einberger, and Brian ***K***ernighan

An `awk` program tends to have the structure of

```
condition { action }
condition { action }
...
```
We will only be dealing with simple conditions here though.

AWK
========================================================
## Our example file - BED
`cpg.bed`

 - Chromosome Name
 - Chromosome Start
 - Chromosome End
 
Optional Fields:

Name, Score, Strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts

___

<center>
![](basic_bed.png)
</center>

AWK
========================================================
## Our example file - BED
`cpg.bed`

```
$ head cpg.bed
chr1    28735   29810   CpG:_116
chr1    135124  135563  CpG:_30
chr1    327790  328229  CpG:_29
chr1    437151  438164  CpG:_84
chr1    449273  450544  CpG:_99
chr1    533219  534114  CpG:_94
chr1    544738  546649  CpG:_171
chr1    713984  714547  CpG:_60
chr1    762416  763445  CpG:_115
chr1    788863  789211  CpG:_28
$ █ 
```

AWK
========================================================

`condition { action }`

The default action in awk is to print the line, so we take advantage of that frequently for filtering a big file.


```bad
columns are identifed by number
      ↓
```
```
awk '$1 == "chr1"' cpg.bed
```

```bad
    ↑            ↑
    single quotes indicate the start and end of the program.
```

Here we are only printing records from `cpg.bed` that are for chr1


AWK can find lines that match your criteria
========================================================

In this case we are searching for lines where the end coordinate is less than the start coordinate.

```
$ awk '$3 < $2' cpg.bed
$ 
```
This returns no lines on a well formed bed file, because the start should come before the end.


AWK knows about line numbers too
========================================================

`NR` is a special variable in awk that stands for record number.

```
$ awk 'NR == 100' cpg.bed
chr1    1417285 1417753 CpG:_42
$ █
```

The above command prints the 100th line in the file.

AWK knows about line numbers too
========================================================

It's possible to have multiple conditions at once too!
in awk, like many languages, the logical AND operator is `&&`

Here we return lines 100 - 200

```
$ awk 'NR >= 100 && NR <= 200' cpg.bed
chr1    1417285 1417753 CpG:_42
chr1    1420226 1420684 CpG:_34
chr1    1444588 1444849 CpG:_19
chr1    1446878 1448205 CpG:_151
chr1    1452701 1453199 CpG:_40
chr1    1455302 1455892 CpG:_48
chr1    1456829 1457074 CpG:_18
chr1    1457954 1459029 CpG:_77
...
```

AWK knows about line numbers too
========================================================

Going further we can add in a logical `OR` operator and group the sub-statements too

```
$ awk '(NR >= 100 && NR <= 200) || ($1 == "chr18)'
chr1    1417285 1417753 CpG:_42
chr1    1420226 1420684 CpG:_34
chr1    1444588 1444849 CpG:_19
chr1    1446878 1448205 CpG:_151
chr1    1452701 1453199 CpG:_40
chr1    1455302 1455892 CpG:_48
chr1    1456829 1457074 CpG:_18
chr1    1457954 1459029 CpG:_77
...
```

The {action} side of AWK
========================================================
`condition { action }`

Rather than just manipulating the condition in which to print lines, we can also
manipulate the output.

We are using the action between curly braces to print the BED record and the length of that interval


```bad
$0 indicates to print the whole line
               ↓
```
```
$ awk '{print $0, $3-$2}' cpg.bed
chr1    28735   29810   CpG:_116 1075
chr1    135124  135563  CpG:_30 439
chr1    327790  328229  CpG:_29 439
...
```

Column separators in AWK
========================================================

By default awk uses a space to separate columns. Most genomics tasks prefer tabs, and keeping them all consistent is better than mixing them.

This also introduces another special variable for AWK, `OFS` which is the Output field separator.

```
$ awk -vOFS='\t' '{print $0, $3-$2}' cpg.bed
chr1    28735   29810   CpG:_116        1075
chr1    135124  135563  CpG:_30 439
chr1    327790  328229  CpG:_29 439
...
```

Column separators in AWK
========================================================

When one sets the `OFS` variable, it only changes for rows that are altered. The annoying sideffect of this is that a command like.
```
awk -vOFS='\t' '{print $0}'
```
does not change the delimiters between all columns in the file like you might expect - since nothing is being changed.


A quick trick to change all the column delimiters at once, is just to reassign a column to itself.
```
awk -vOFS='\t' '$1=$1'
```
This alters every row, and replaces all the column delimiters with tabs

Using AWK as an improved cut command
========================================================

We can use awk to print columns in any order, whereas `cut` prints them only in the original order.

```
$ cut -f4,1 cpg.bed
chr1    CpG:_116
chr1    CpG:_30
chr1    CpG:_29
...
```

With awk we print in arbitrary order, depending on how we list the columns in the print statement.
```
$ awk -vOFS='\t' '{print $4, $1}'
CpG:_116        chr1
CpG:_30 chr1
CpG:_29 chr1
...
```

How many columns on each line?
========================================================

One more helpful built in variable we have in awk is `NF`. This represents the number of fields/columns are detected for each line. By default, these are just white-space delimited.

```
$ awk '{print NF}' cpg.bed
4
4
4
...
```

Assigning variables with AWK
========================================================

Sometimes when ones computation get's a little more complex, its nice to be able to assign your own variables, and not just use the ones built into awk.


This is a simple example that still just prints the length of the element after printing the line.


```bad
len is a variable that is assigned the output of $3-$2
                    ↓
```
```
$ awk -vOFS='\t' '{len=($3-$2); print $0, len}' cpg.bed
chr1    28735   29810   CpG:_116        1075
chr1    135124  135563  CpG:_30 439
chr1    327790  328229  CpG:_29 439
```

Computing across a whole file
========================================================

We can compute the total number of base pairs represented by the CpG islands in two parts.

```
$ awk '{tot_len += $3-$2}END{print tot_len}' cpg.bed
21842742
$ █
```

First, we start with a variable `tot_len` that starts at 0 but increases by the length of each CpG island as awk proceeds down the file.

We then have an `END` statement that indicates that that the following part occurs after all the line by line processing is finished.

Second, we print the contents the `tot_len` variable

Computing across a whole file
========================================================

We can use one of the built in variables in our calculation too.

Here we calculate the average number of base pairs represented by a CpG island.

```
$ awk '{tot_len += $3-$2}END{print tot_len/NR}' cpg.bed
761.31
$ █
```

Computing across a whole file
========================================================

We can even - using [arrays in AWK](https://www.gnu.org/software/gawk/manual/html_node/Arrays.html) - 
compute something per unique value in a column.

Here we calculate the number of base pairs per chromatin annotation across the file.
This is advanced usage, but can be useful to add to ones toolbox.


```bad
array called "a" indexed by the values in column 4
      ↓        
```
```
awk '{a[$4] += ($3-$2)} END {for (i in a) print i, a[i]}' hesc.chromHmm.bed | sort -k2,2nr
```

```bad
                  ↑
we increment the value for each index by the length of each element of that index
```

Computing across a whole file
========================================================

We can even - using [arrays in AWK](https://www.gnu.org/software/gawk/manual/html_node/Arrays.html) - 
compute something per unique value in a column.

Here we calculate the number of base pairs per chromatin annotation across the file.
This is advanced usage, but can be useful to add to ones toolbox.

```
$ awk '{a[$4] += ($3-$2)} END {for (i in a) print i, a[i]}' hesc.chromHmm.bed | sort -k2,2nr
13_Heterochrom/lo 1992618522
11_Weak_Txn 503587684
10_Txn_Elongation 82134508
7_Weak_Enhancer 66789380
12_Repressed 38190271
6_Weak_Enhancer 35705628
9_Txn_Transition 25674539
8_Insulator 22467670
2_Weak_Promoter 18737021
3_Poised_Promoter 18724570
1_Active_Promoter 9908594
5_Strong_Enhancer 7139201
14_Repetitive/CNV 4079101
4_Strong_Enhancer 2648615
15_Repetitive/CNV 2342560
$ █
```
