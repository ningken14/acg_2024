# Assignment 2: More unix commands and a little git

```put all commands and terminal output in triple backquotes```

## Git

1. Either download or create a copy of this file in the repository we created in class. "Add" the file by checking the appropriate box in the git tab of rstudio. Then commit and push to remote. Make a commit message that records the fact that this file has not been worked on (by you) yet.

2. Open the newly created homework file in your homework repo within rstudio. Complete the rest of this assignment. When you're finally ready to finally hand it in, add, commit (with appropriate message) and push. If you need to make an edit before the deadline, you can repeat the commit cycle making a comment in the commit message.

## Man pages

Nearly every command in GNU has a manual, documentation that tells you all the functions and capabilities of each command that we've been learning about. Naturally, there is a program called `man` to access the manual pages (more commonly "man pages") for each command.

3. Let's use `man` to figure out what some other commands do. Run `man` on the `who` command. What do you think the SYNOPSIS section is for?

4. Which option of `who` allows you to determine when your system was last rebooted? When was your system last booted? Show the shell interaction inside triple quotes as we did in the previous assignment.

5. Using `man`, see if you can determine what the `cut` program does? (no answer required)

6. In a web browser, navigate to the following page on Github, which has lots of example datasets people use for machine learning: https://github.com/EpistasisLab/pmlb/tree/master/datasets
Download the dataset called `breast_cancer` (be exact: there are several breast cancer datasets!). You should have a spread sheet in "tab separated values" format (".tsv"). First take a peek at the data using `less`. Use `head` to display the first line only, paste the output of the terminal interaction below.

7. How many fields are there? Which number field corresponds to tumor sizes? Can you write a command to extract this column? Show the command below, but not the output.

8. What does the `wc` command do?

9. How many lines are there in the dataset? How many words? (show work)

10. How many characters are there in the column 'tumor-size'? (hint: use `cut` to isolate the tumor size column, and a pipe `|` to direct the output of this command to an appropriate call of `wc`) (show work)

11. The `sort` command is very useful. Use the `man` page to determine how you would find all the unique values that "tumor-size" can have in this dataset. Using a strategy similar to the previous problem, what are all the unique values in the tumor-size column? (hint: string together 3 commands with 2 pipes) (show work)