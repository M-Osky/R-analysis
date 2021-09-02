# R analysis
Some short R scripts I wrote to analyse my data.
This is just a sample of some of the work I do with R either for our own research or for the students.

ALL SCRIPTS HAVE EITHER SOME LINES ON THE TOP DESCRIVING THEM OR HELP INFORMATION (--help) IF THEY ARE DONE TO WORK FROM COMMAND LINE

# I will be uploading and editing this continuously, some scripts can be run in desktop R, some other are set up to run in our unix server (using command line).

# Will add descriptions shortly here



Parametric tests for Normal data:

t-test if grouping variables with only two factors.

One way ANOVA for all grouping variables with more than 2 factors
  - If ANOVA is significant: pair-wise pos-hoc tukey test
  - If ANOVA is non-significant: pair-wise Duncan test

Plots will be saved automatically. If there is 4 factors bars between the boxes will show if there is significant differences between groups. For more factors letters in the top of the boxes will mark which ones are not significantly different.

Works with multiple files with multiple grouping variables, but all of them must have the same data structure.

I have new versions of similar analysis for non-parametric and for multivariant analysis (the previous versions are now at R-analysis/old/) will upload when properly commented.

	id	extrainfo	pop	sex	var1	var2	var3	var4	var5	var6
	Lem042	LemuriaNA	LM	fem	0004	0008	0015	0016	0023	0042
	Atl016	Atlantis	AT	fem	0003	0014	0015	0092	0062	0036
	Atl028	AtlantidaSE	AT	mal	0013	0021	0034	0055	0089	0144
  
