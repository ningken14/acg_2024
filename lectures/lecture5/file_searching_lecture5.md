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

Lecture 5 - Working with data on the command line
========================================================
author: Simon Coetzee
date: 02/06/2024
width: 1440
height: 900
transition: fade


## **Searching NGS File Formats**

<small>This work, "Searching NGS File Formats", is a derivative of ["Intro to Shell"](https://hbctraining.github.io/Intro-to-Shell/) by hbctraining, used under [CC BY](https://creativecommons.org/licenses/by/4.0/) and is a derivative of part of ["Applied Computational Genomics Course at UU"](https://github.com/quinlan-lab/applied-computational-genomics) by quinlan-lab. "Searching NGS File Formats" is licensed under [CC BY](https://creativecommons.org/licenses/by/4.0/) by Simon Coetzee.</small>

Goals (grep and redirection)
========================================================

- Search for characters or patterns in a text file using the grep command

Goals (grep and redirection)
========================================================

- Search for characters or patterns in a text file using the grep command

- Write to and append a file using output redirection

Goals (grep and redirection)
========================================================

- Search for characters or patterns in a text file using the grep command

- Write to and append a file using output redirection

- Use the pipe `|` character to chain together commands

GREP
========================================================

`grep` is a command-line utility for searching plain-text data sets for lines 
that match a *regular expression*. Its name comes from the `ed` command `g/re/p` 
(globally search for a regular expression and print matching lines), 
which has the same effect

`grep` syntax looks like
```
grep search_term filename
```

FASTQ
========================================================

Get some FASTQ data:

```
$ wget https://cloud.coetzee.me/s/zsYsQEyjk84aCYt/download/data.tar.gz
$ tar xvf data.tar.gz
$ mkdir ~/lecture5
$ cp ~/data/*.fq ~/lecture5
$ ls ~/lecture5
Mov10_oe_1.subset.fq  Mov10_oe_2.subset.fq  Mov10_oe_3.subset.fq
$ head ~/lecture5/Mov10_oe_1.subset.fq
@HWI-ST330:304:H045HADXX:1:1101:1162:2055
NAGAACTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTT
+
#1=DDFFFHHHGHIJJJJIJJJGEGGAFGBHHEHGFBFFDEDECDDA==CB@BDDDDD?;B-<CBDDD>BBBBDDB5<@DDDCDDB@-9ACDDDDB?B<?
@HWI-ST330:304:H045HADXX:2:2111:20110:84312
GTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAG
+
@@<FFFFDFFH>DEGFEGIJGJIJD9;CFCG;@;9?DDCD8AHGEF@84ADB?CD>3@CAACBBBDD@@@??90))5055(22-95<-5(:<ACBB@?8?
@HWI-ST330:304:H045HADXX:1:1214:9417:35291
CTCCAGACTCCGATCGTACAGCTTGAACTTCACATCTGAGGGCAGCAACGAGACCCCACGGGAGGCCACAGGAAAAAGCATGGGCCATAGCACCCAGCGC
```

grep pattern FASTQ
========================================================

Let's say we consider "bad" reads from our file as those that contain 10
consecutive Ns. We can `grep` for that.


```bash
$ grep NNNNNNNNNN Mov10_oe_1.subset.fq
```

This will return many lines of text!

But only the actual read, what if we wanted to see the whole record, all four lines?

grep pattern FASTQ
========================================================

To see what command line arguments we can use to modify `grep`s behavior, we can enter the following in your open bash shell

```{}
$ man grep
```
```{}
GREP(1)                User Commands                GREP(1)

NAME
       grep,  egrep,  fgrep, rgrep - print lines that match
       patterns

SYNOPSIS
       grep [OPTION...] PATTERNS [FILE...]
       grep [OPTION...] -e PATTERNS ... [FILE...]
       grep [OPTION...] -f PATTERN_FILE ... [FILE...]

DESCRIPTION
       grep searches for PATTERNS in each  FILE.   PATTERNS
```

grep pattern FASTQ
========================================================

Given that we are looking for four lines (one record), we want the line before the read, and two lines after the read.

```{}
-A NUM, --after-context=NUM
  Print NUM lines  of  trailing  context  after
  matching  lines.   Places a line containing a
  group  separator  (--)   between   contiguous
  groups   of   matches.    With   the   -o  or
  --only-matching option, this  has  no  effect
  and a warning is given.
-B NUM, --before-context=NUM
  Print  NUM  lines  of  leading context before
  matching lines.  Places a line  containing  a
  group   separator   (--)  between  contiguous
  groups  of   matches.    With   the   -o   or
  --only-matching  option,  this  has no effect
```

The `-B` and `-A` arguments for `grep` will be useful to return the matched line plus one before `-B 1` and two lines after `-A 2`.

grep pattern FASTQ
========================================================

The `-B` and `-A` arguments for `grep` will be useful to return the matched line plus one before `-B 1` and two lines after `-A 2`.

```{}
$ grep -B 1 -A 2 NNNNNNNNNN Mov10_oe_1.subset.fq
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
--
@HWI-ST330:304:H045HADXX:1:1101:1106:89824
CACAAATCGGCTCAGGAGGCTTGTAGAAAAGCTCAGCTTGACANNNNNNNNNNNNNNNNNGNGNACGAAACNNNNGNNNNNNNNNNNNNNNNNNNGTTGG
+
?@@DDDDDB1@?:E?;3A:1?9?E9?<?DGCDGBBDBF@;8DF#########################################################
--
```

grep pattern FASTQ
========================================================

The `-B` and `-A` arguments for `grep` will be useful to return the matched line
plus one before `-B 1` and two lines after `-A 2`. The argument
`--no-group-separator` will make sure that the `--` group separator lines 
between results are removed, so the output appears in the same format as the input.

```{}
$ grep -B 1 -A 2  --no-group-separator NNNNNNNNNN Mov10_oe_1.subset.fq
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
@HWI-ST330:304:H045HADXX:1:1101:1106:89824
CACAAATCGGCTCAGGAGGCTTGTAGAAAAGCTCAGCTTGACANNNNNNNNNNNNNNNNNGNGNACGAAACNNNNGNNNNNNNNNNNNNNNNNNNGTTGG
+
?@@DDDDDB1@?:E?;3A:1?9?E9?<?DGCDGBBDBF@;8DF#########################################################
```

Invert grep
========================================================
Another helpful argument for `grep` is the `-v` argument:

```{}
-v, --invert-match
      Invert  the sense of matching, to select non-
      matching lines.
```

basically returning everything that does not match the pattern.

we can examine a sample file from `/data` directory.


Invert grep
========================================================

```
$ cp ~/data/sample_sheet.txt ~/lecture5
$ cd ~/lecture5
$ cat sample_sheet.txt
sample  treatment
AR.1    vehicle
AR.2    vehicle
AR.3    vehicle
AR.DHT.1    DHT
AR.DHT.2    DHT
AR.DHT.3    DHT
```

we can see the metadata file for a hypothetical chipseq experiment.

Invert grep
========================================================

If we don't want to see the _vehicle_ samples we can use the `-v` argument like this:
```{}
$ grep -v vehicle sample_sheet.txt
sample  treatment
AR.DHT.1    DHT
AR.DHT.2    DHT
AR.DHT.3    DHT
```
It returns all the lines in the file **except** the lines that contain _vehicle_.

Redirecting Output To A File
========================================================
## The command to write something to a file is `>`

We can use this to, instead of writing the output of a command to a terminal where we see it, write to a file to store that output.
For example, we can store all the sequences that contain `NNNNNNNNNN` to a new file called `bad_reads.fq`

```{}
$ grep -B 1 -A 2 --no-group-separator NNNNNNNNNN Mov10_oe_1.subset.fq > bad_reads.fq
```

A new file is created called `bad_reads.fq`, and nothing is displayed on the screen.

We can confirm the file was created by entering:
```{}
$ ls -l
```

Redirecting Output and Appending To A File
========================================================
## The redirection command for appending something to an existing file is `>>`

If we were to execute the previous command on a new file `Mov10_oe_2.subset.fq`


```bad
$ grep -B 1 -A 2 --no-group-separator NNNNNNNNNN Mov10_oe_2.subset.fq > bad_reads.fq
```

It would overwrite the results, replacing them with new reads.

However, if we would like to append the bad reads from the new file to the end of the `bad_reads.fq` file, we would use the `>>` command instead.

```{}
$ grep -B 1 -A 2 --no-group-separator NNNNNNNNNN Mov10_oe_2.subset.fq >> bad_reads.fq
$ ls -l bad_reads.fq
```

Redirecting Output To Another Command
========================================================
IMO one of the coolest parts about the *nix environment.

## The redirection command for using the output of a command as input for a different command is `|` (pipe).

<center>
![](keyboard_pipe.png)
</center>

Redirecting Output To Another Command
========================================================

We can pipe the output of a command to `less` for example to be able to scroll through it at our leisure rather than watching it whiz by.

```
$ grep -B 1 -A 2 NNNNNNNNNN Mov10_oe_1.subset.fq | less
```

Or if we only want to see the first few lines of output, we could pipe it into `head`

```
$ grep -B 1 -A 2 NNNNNNNNNN Mov10_oe_1.subset.fq | head -n 8
```

Redirecting Output To Another Command
========================================================

Or maybe we only want to know _how many_ lines from a file match our pattern. We could use the command `wc`.

`wc` stands for ***w***ord ***c***ount, and can count the number of words, lines, and characters in a piece of text given to it.

the argument `-l` indicates to count only the lines.

see `man wc` like we did with grep before to see more options.

```
$ grep NNNNNNNNNN Mov10_oe_1.subset.fq | wc -l
```

Searching other files
========================================================
Lets look at a gene annotation file like the GTF we saw before

<center>
![](color_gtf.png)
</center>

Searching GTF
========================================================
We will download a sample gtf file from gencode:

```{}
$ cd ~
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
$ gunzip gencode.v39.annotation.gtf.gz
$ mv ~/gencode.v39.annotation.gtf ~/lecture5
$ cd ~/lecture5
$ less gencode.v39.annotation.gtf
```
```{}
##description: evidence-based annotation of the human genome (GRCh38), version 39 (Ensembl 105)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2021-09-02
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 3; exon_id "ENSE00002312635.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
```

Some considerations when searching with grep
========================================================
Searching for the word "gene" in a gtf file.

```
$ grep "gene" gencode.v39.annotation.gtf | head
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 3; exon_id "ENSE00002312635.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	transcript	12010	13670	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
chr1	HAVANA	exon	12010	12057	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; exon_number 1; exon_id "ENSE00001948541.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
chr1	HAVANA	exon	12179	12227	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; exon_number 2; exon_id "ENSE00001671638.2"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
chr1	HAVANA	exon	12613	12697	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; exon_number 3; exon_id "ENSE00001758273.2"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
chr1	HAVANA	exon	12975	13052	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; exon_number 4; exon_id "ENSE00001799933.2"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
```

It occurs everywhere - even when it's not the type of interval we are looking for

Some considerations when searching with grep
========================================================
We can search for something that looks just like the "word" gene with the `-w` flag

```
$ grep -w "gene" gencode.v39.annotation.gtf | head
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
chr1	HAVANA	gene	14404	29570	.	-	.	gene_id "ENSG00000227232.5"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; level 2; hgnc_id "HGNC:38034"; havana_gene "OTTHUMG00000000958.1";
chr1	ENSEMBL	gene	17369	17436	.	-	.	gene_id "ENSG00000278267.1"; gene_type "miRNA"; gene_name "MIR6859-1"; level 3; hgnc_id "HGNC:50039";
chr1	HAVANA	gene	29554	31109	.	+	.	gene_id "ENSG00000243485.5"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; level 2; hgnc_id "HGNC:52482"; tag "ncRNA_host"; havana_gene "OTTHUMG00000000959.2";
chr1	ENSEMBL	gene	30366	30503	.	+	.	gene_id "ENSG00000284332.1"; gene_type "miRNA"; gene_name "MIR1302-2"; level 3; hgnc_id "HGNC:35294";
chr1	HAVANA	gene	34554	36081	.	-	.	gene_id "ENSG00000237613.2"; gene_type "lncRNA"; gene_name "FAM138A"; level 2; hgnc_id "HGNC:32334"; havana_gene "OTTHUMG00000000960.1";
chr1	HAVANA	gene	52473	53312	.	+	.	gene_id "ENSG00000268020.3"; gene_type "unprocessed_pseudogene"; gene_name "OR4G4P"; level 2; hgnc_id "HGNC:14822"; havana_gene "OTTHUMG00000185779.1";
chr1	HAVANA	gene	57598	64116	.	+	.	gene_id "ENSG00000240361.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "OR4G11P"; level 2; hgnc_id "HGNC:31276"; havana_gene "OTTHUMG00000001095.3";
chr1	HAVANA	gene	65419	71585	.	+	.	gene_id "ENSG00000186092.7"; gene_type "protein_coding"; gene_name "OR4F5"; level 2; hgnc_id "HGNC:14825"; havana_gene "OTTHUMG00000001094.4";
chr1	HAVANA	gene	89295	133723	.	-	.	gene_id "ENSG00000238009.6"; gene_type "lncRNA"; gene_name "ENSG00000238009"; level 2; tag "overlapping_locus"; havana_gene "OTTHUMG00000001096.2";
```
limiting the types of lines we get back.

GREP -w 
========================================================
```
       -w, --word-regexp
              Select only those lines containing matches
              that form whole words.  The test  is  that
              the  matching  substring must either be at
              the beginning of the line, or preceded  by
              a    non-word    constituent    character.
              Similarly, it must be either at the end of
              the   line   or  followed  by  a  non-word
              constituent  character.   Word-constituent
              characters  are  letters,  digits, and the
              underscore.  This option has no effect  if
              -x is also specified.
```
Some considerations when searching with grep
========================================================
Same goes for when one is inverting a search.
```
$ grep -v "gene" gencode.v39.annotation.gtf
##description: evidence-based annotation of the human genome (GRCh38), version 39 (Ensembl 105)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2021-09-02
```

Some considerations when searching with grep
========================================================
Almost every line as "gene" somewhere in the line.
```
$ grep -v -w "gene" gencode.v39.annotation.gtf | head
##description: evidence-based annotation of the human genome (GRCh38), version 39 (Ensembl 105)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2021-09-02
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 3; exon_id "ENSE00002312635.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	transcript	12010	13670	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:37102"; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000002844.2";
```


Some considerations when searching with grep
========================================================
One must be careful that included in the results are not just areas where you have a substring.

`-w` can be helpful here too
```
$ grep -w gene gencode.v39.annotation.gtf | grep "ATP5MF"
chr6	HAVANA	gene	108907615	108907873	.	-	.	gene_id "ENSG00000228834.1"; gene_type "processed_pseudogene"; gene_name "ATP5MFP2"; level 1; hgnc_id "HGNC:21540"; tag "pseudo_consens"; havana_gene "OTTHUMG00000015333.2";
chr7	HAVANA	gene	99419749	99466197	.	-	.	gene_id "ENSG00000248919.7"; gene_type "protein_coding"; gene_name "ATP5MF-PTCD1"; level 1; hgnc_id "HGNC:38844"; tag "overlapping_locus"; havana_gene "OTTHUMG00000160779.2";
chr7	HAVANA	gene	99448475	99466186	.	-	.	gene_id "ENSG00000241468.8"; gene_type "protein_coding"; gene_name "ATP5MF"; level 2; hgnc_id "HGNC:848"; tag "overlapping_locus"; havana_gene "OTTHUMG00000154609.7";
chr9	HAVANA	gene	77040416	77040673	.	+	.	gene_id "ENSG00000232851.3"; gene_type "processed_pseudogene"; gene_name "ATP5MFP3"; level 1; hgnc_id "HGNC:21286"; tag "pseudo_consens"; havana_gene "OTTHUMG00000020052.2";
chr12	HAVANA	gene	6270168	6270425	.	-	.	gene_id "ENSG00000256103.2"; gene_type "processed_pseudogene"; gene_name "ATP5MFP5"; level 1; hgnc_id "HGNC:33611"; tag "pseudo_consens"; havana_gene "OTTHUMG00000168354.1";
chr15	HAVANA	gene	66414933	66415198	.	-	.	gene_id "ENSG00000261102.2"; gene_type "processed_pseudogene"; gene_name "ATP5MFP6"; level 1; hgnc_id "HGNC:33612"; tag "pseudo_consens"; havana_gene "OTTHUMG00000172827.2";
chr17	HAVANA	gene	64491523	64491789	.	+	.	gene_id "ENSG00000256826.1"; gene_type "processed_pseudogene"; gene_name "ATP5MFP4"; level 1; hgnc_id "HGNC:32451"; tag "pseudo_consens"; havana_gene "OTTHUMG00000132312.1";
chr21	HAVANA	gene	36388878	36389112	.	-	.	gene_id "ENSG00000224421.1"; gene_type "processed_pseudogene"; gene_name "ATP5MFP1"; level 1; hgnc_id "HGNC:849"; tag "pseudo_consens"; havana_gene "OTTHUMG00000086613.1";
```

Some considerations when searching with grep
========================================================
One must be careful that included in the results are not just areas where you have a substring.

`-w` can be helpful here too ... but not always.
```
$ grep -w gene gencode.v39.annotation.gtf | grep -w "ATP5MF"
chr7	HAVANA	gene	99419749	99466197	.	-	.	gene_id "ENSG00000248919.7"; gene_type "protein_coding"; gene_name "ATP5MF-PTCD1"; level 1; hgnc_id "HGNC:38844"; tag "overlapping_locus"; havana_gene "OTTHUMG00000160779.2";
chr7	HAVANA	gene	99448475	99466186	.	-	.	gene_id "ENSG00000241468.8"; gene_type "protein_coding"; gene_name "ATP5MF"; level 2; hgnc_id "HGNC:848"; tag "overlapping_locus"; havana_gene "OTTHUMG00000154609.7";
```

Some considerations when searching with grep
========================================================
We can try to include the `"` marks in our search to isolate just the string `"ATP5MF"` rather than `ATP5MF`
```
$ grep -w gene gencode.v39.annotation.gtf | grep \"ATP5MF\"
chr7	HAVANA	gene	99448475	99466186	.	-	.	gene_id "ENSG00000241468.8"; gene_type "protein_coding"; gene_name "ATP5MF"; level 2; hgnc_id "HGNC:848"; tag "overlapping_locus"; havana_gene "OTTHUMG00000154609.7";
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

Bonus: Cheating at wordle
========================================================

<center>
![](wordle.png)
</center>

```
$ cat /usr/share/dict/words | egrep '^m...h$' | egrep -v 'r|a|i|s|e|u|l|c'
month
```

Bonus: Cheating at wordle
========================================================
## Doesn't always work though

<center>
![](wordle2.png)
</center>
___
Starting at **CANTO**
```
cat /usr/share/dict/words | egrep '^.a...$' | egrep -v 'r|i|s|e|n|o' | grep c | grep t
batch
caput
catch
catty
hatch
latch
match
patch
tacky
watch
yacht
```


Two More Commands
========================================================

## `cut` is a command that extracts columns from files
## `sort` is a command used to sort the contents of a file in a particular order.

these commands both make use of columns of data.

cut To Extract Columns
========================================================

By default, `cut` expects columns in files to be separated by tabs. Like the GTF file above, and many other other common bioinformatics formats.
However, with the argument `-d` we can specify a different delimiter. Some files like `.csv` files can have thier columns separated by a comma (`.csv` stands for comma separated values) and so we can cut columns from them with `cut -d","`

Here we can examine the chromosome and start coordinates from the GTF file.

The `-f` argument allows us to chose the ***f***ield or column we want to see

```{}
$ cut -f1,4 gencode.v39.annotation.gtf | head
##description: evidence-based annotation of the human genome (GRCh38), version 39 (Ensembl 105)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2021-09-02
chr1	11869
chr1	11869
chr1	11869
chr1	12613
chr1	13221
```

cut To Extract Columns
========================================================

By default, `cut` expects columns in files to be separated by tabs. Like the GTF file above, and many other other common bioinformatics formats.
However, with the argument `-d` we can specify a different delimiter. Some files like `.csv` files can have thier columns separated by a comma (`.csv` stands for comma separated values) and so we can cut columns from them with `cut -d","`

Here we can examine the chromosome and start coordinates from the GTF file.

The `-f` argument allows us to chose the ***f***ield or column we want to see

```{}
$ cut -f1,4 gencode.v39.annotation.gtf | grep -v # | head -4
chr1	11869
chr1	11869
chr1	11869
chr1	12613
```

sort To Sort based on Columns
========================================================

By default, `sort` sorts a text file by alpha/numeric order, across all columns.
It is able to, after sorting, showing only unique columns with the argument `-u`.

```{}
$ cut -f1,4 gencode.v39.annotation.gtf | wc -l
3241007
$ cut -f1,4 gencode.v39.annotation.gtf | sort -u | wc -l
624652
```

sort To Sort based on Columns
========================================================

`sort` is also able to sort files based on specific columns with the `-k` 
argument, and therein able to sort by specific methods.

`general-numeric -g, human-numeric -h, month -M, numeric -n, random -R, version -V`

Here we sort column 4 numerically.
```
sort -k4,4n gencode.v39.annotation.gtf | head -5
```
Here we sort column 4 reverse numerically.
```
sort -k4,4nr gencode.v39.annotation.gtf | head -5
```

sort Chromosome names
========================================================
Lets make a file of chromosome names to see how some of the sorts work:

Also just to place things on multiple lines, here we demonstrate that the `\` 
can allow you continue on the next line without executing your command. This can be helpful for readability purposes.

```
$ cut -f1 ~/lecture5/gencode.v39.annotation.gtf | \
> grep -v "#" | \
> sort -u > ~/lecture5/chromosome_names.txt
$ cd ~/lecture5
```


sort Chromosome names
========================================================

```{}
$ shuf chromosome_names.txt | sort -n
chr1
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr2
chr20
chr21
chr22
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chrX
chrY
```

-----------

```{}
shuf chromosome_names.txt | sort -V
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
```
