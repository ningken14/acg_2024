## Homework 1

---

Open this file in Rstudio if you haven't already! Save frequently as you work.

In this assignment we will practice basic skills from week 2 to populate your virtual linux machine with some typical structure. Read carefully as the instructions have to be followed to the letter to get the right answers. 

You may edit this document directly in Rstudio. When the assignment is finished, email Dr. Hazelett with an edited copy of this file as an attachment making sure to do the last problem which gives your homework a uniquely traceable name. Next week we will begin submitting all homework on github.

To display a terminal interaction in this markdown when you hand it in, simply paste the entire contents between triple single backquotes as follows:

```
hazelettd@54b2bfd85cd6:~$ passwd
Changing password for hazelettd.
Current password: 
New password: 
Retype new password: 
passwd: password updated successfully
```

1. Navigate to your user's home directory. How do you know you're there? What is the shortest command you can type to get you there?

2. In the home directory, create a subdirectory (folder) called 'Desktop'. Now create other folders called 'Projects', 'Music', 'Downtown', and 'Pictures'.

3. The `Downtown` folder was mistakenly named. What are two alternative methods (commands) you could use to correct your mistake and rename the folder? Rename the folder as `Downloads`.

4. Copy the contents of /data into a new folder of the same name on your home directory, _in the 'Desktop' folder_. Hint: A wildcard '*' could be used to copy all the folder's contents at once.

5. Use the 'ls' command to determine the size of the fastq files (only!) in your new ~/Desktop/data folder. If you're not sure which options to use to include size, use the 'man' tool. Can you confirm that the fastq files are the same size as the originals you copied from the root directory? (Show your terminal interactions)

6. The keyword '..' denotes "one folder up in the filesystem tree". Navigate to the Desktop folder you created. Type the following command at your terminal prompt:

`mv data ..`

What is the result of this command?

7. From within the Desktop directory (you may have to navigate there if you changed your present working directory to answer the previous question) move the data folder and its contents back into 'Desktop' using a pathname composed with '..' and using '.' as the target location. '.' is shorthand for "the present working directory".

8. Navigate to your home directory. Type the following command _exactly_:

`rmdir Desktop/data`

From the error message, can you explain why the command didn't work?

9. Now try the following command:

`rm -rf Desktop/data`

What do the '-r' and '-f' flags do in this instance? (Hint: consult the 'man' tool)

10. Use 'man' to figure out which options to re-copy *recursively* the data directory from the root /data and all its contents _in a single command_ back into the home directory ('~/'). Use ls to confirm that it worked, and then remove the directory and its contents once more (we need to conserve space for now and don't want to keep 2 copies after this assignment as some of these files are large).

Bonus: As we discovered in lecture, /data files are off limits to editing. Are you able to edit the contents of files in the copied directory? Why or why not?:

11. Using a terminal window, rename this file _carefully_ by appending your full last name and an underscore character in front of the existing file name. Paste the full shell interaction below (as described at the top in directions). Here is an example of what the new filename should look like:

`Hazelett_assignment_01.Rmd`
