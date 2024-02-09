# Assignment 4: cut, sort, & advanced grep

##### Due February 16th, 2024  

---

##### <i>Directions  
Show all work. Write complete and coherent answers. Correct answers alone will yield partial credit. Each question is worth 4 points unless otherwise indicated (total 40 pts).  </i>

<b><i>Put single lines of shell code inside single quotes, e.g.:</b></i>  

`du -skh / > ~/hd_usage.txt`

<b><i>Put all commands and terminal output inside triple backquotes to show a bash shell interaction:  </b></i>

```
$ echo "this is terminal output"
this is terminal output
$
```

---

## using cut and sort  

##### 1. Delimiters. In the data folder we downloaded in lecture 5, there's a .gtf file containing the annotations of chromosome 1 from HG19. How many tab-delimited fields are there in this file? (hint: In Rstudio go to Edit>Settings. Choose "Code", and then the "Display" tab. Check the box next to "Show whitespace characters".)  

a. The GTF file is too large to load in Rstudio. Use `head` and redirect `>` to make a snippet of the file called `snippet.gtf` with 20 lines.  
b. Open `snippet.gtf` in Rstudio. What types of whitespace characters can you find?  
c. How many tab-delimited fields are there?  

##### 2. Gene Names.  

a. What field of `snippet.gtf` are gene names in?  
b. Write a pipe between two `cut` commands to extract the subfield containing gene names. You may need to consider using different whitespace characters as delimiters.  
c. Use a similar strategy combined with `sort` and `wc` to determine the number of unique gene names on chromosome 1 with the original .gtf file you used to create the snippet version.  

##### 3. Annotations.  

a. What are the unique annotations present in field 3 of the chromosome 1 file?  
b. How many of each annotation are present in the file?  
c. What are the unique annotations in field 5 of the last field (9) of the chromosome 1 gtf file?  

---

## grep  

##### 4. How many distinct WNT- prefix genes are on chromosome 1?  

##### 5. How many distinct transcripts are there associated with WNT genes from chromosome 1?  

##### 6. How many distinct genes are encoded on the '+' strand of chromosome 1?

##### 7. How many WNT genes are on the '+' strand?

---

## advanced grep

Review the advanced grep slides (slide 29-32 of lecture 5).  

```
^: matches pattern at start of string
$: matches pattern at end of string
.: matches any character except new lines
[]: matches any of enclose characters
[^]: matches any characters *except* ones enclosed (note: is different from ^)
\: "escapes" meta-characters, allows literal matching
```

##### 8. Cheat at Wordle. First guess is "METER", M, T, the second E, and R are green. Formulate an appropriate grep expression (you may chain more than one together with a pipe if needed) to identify candidate 5-letter words for a second guess.  

---

##### Advanced grep: Grouping.  

You can also group extended regular expressions together using parentheses, _e.g._:  

`egrep -w '(^....p$)|(^m...r$)'`  

This strategy is not required but may be helpful to solve the following problems.  

##### 9. Cheat at Wordle (challenge). The first guess is "ARISE". R, S, and E are yellow. How many valid guesses are there on your computer's (or VM's) `/usr/share/dict/words` file? Careful: Proper nouns are not allowed in Wordle.  

##### 10. How many of those guesses begin with the letter T?  

---

### EXTRA CREDIT

##### A. Using any of the methods you've learned, how many variants in the file `gwas.bed` start with the letters exactly matching "rs1234"? (2 pts)  

##### B. Solve any Wordle puzzle using grep. Show your guesses, the responses, and the commands you use to solve it. (2 pts)  
