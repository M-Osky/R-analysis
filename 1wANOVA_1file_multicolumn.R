#R sucks!

#script to do simple ANOVA with one file containing of all the variables at once.
#Doesn't matter the number of variables



getwd()
#Change working directory
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/AdultsPKPM/Analysis/variances/ANOVA_Ttest")

getwd()			#Check


library(beepr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(ggsignif)



#INPUT FILE:
# The first row can be a title o comment preceded by "#" 
# Second row contains the name/tags of the columns
# First columns: individual or non-usable data. In the example there is two ("indiv" and "label_info") but there can be as many as you want
# further below you will need to set this as "otherstuff", in the example -> otherstuff = 2
# Next: Columns with grouping variables (factor), In the example there is three (sex, pop and group, from column 3 to 5)
# further below you'll need to set how many are in "groupvar", (3 in the example)
# Finally: Columns with quantitative variables, as many as desired, they will all taken into account if they are at the right of the other columns with other information.
# The script will do the tests and plots for every variable

#Something like this:
#
#		#Datafile of phenotypes 1
#		indiv,	label_info,	 sex, pop, group, measure1, m2, m3, m4, m5 ... mn
#		001,	PMf2018_001, fem, PM,  pmf,			04, 08, 15, 16, 23 ... 42
#		002,	PKm2016_002, mas, PK,  pkm,			13, 12, 24, 42, 15 ... 18
#		etc...
#
# Don't mind the tabulation, the  tabs are there only for the data to mach the headers
#Load the packages, set the input file, check the input file and run the script
##The script will do an ANOVA for all grouping variables with more than 2 factors, and a t-test if only two and between each pair of factors.
# the results and plots will be automatically saved, the plots display the box (median and quartiles), a line for the total average (mean) and the outliers.
# the p-value of ANOVA is displayed if more than two factors
# If there is significant differences (t-test) between one of the groups (factors) and the rest, "*" will appear in the top of the box, otherwise "ns"
#Run all the script from "#RUN THE ANALYSIS" till "COSTUMIZED PLOT"
#If you want a beautiful plot for the ANOVA showing the t-test significance of each pair of comparisons you will need to set it manually at the end.
#This is because the scatter of the dots and position of the bar comparing each pair of plots changes according with the magnitude of the measures



rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets


measurespod = read.table("Podadults2017residuals.csv", sep = ',', header =TRUE)	#read the data into a table, read headers
str(measurespod)

#numer of individual or non-relevant data (columns to the left to ignore)
otherstuff = 2
#number of grouping variables. Columns between the ones to ignore and the ones with the variables to analyze
#every cathegory in those columns will be plotted in a different box
groupvar = 3









################################## Below this don't change anythink unless you know what you are doing

coltags <- names(measurespod)
maxcol <- ncol (measurespod)								#number of columns
lastgrup = otherstuff + groupvar							#column with the last grouping variable
firstvar = lastgrup + 1										# column with the first variable
firstgrup = otherstuff+1
numvar = maxcol - lastgrup
varinfo <- measurespod[c(firstvar:maxcol)]									#take a subset including the data of the measured variables
groupinfo <- measurespod[c(firstgrup:lastgrup)]									#take a subset with the grouping variables to test
groupnames <- names (groupinfo)
colgroup = c(firstgrup:lastgrup)
colvariab = c(firstvar:maxcol)
varnames <- names(varinfo)

#CHECK BEFORE PROCEEDING
cat ("Check that this is alright:\n\t", maxcol, "columns in TOTAL\n\t", otherstuff, "with individual/non-relevant data\n\t", groupvar, " columns with grouping variables:", groupnames, "\n\t", numvar, "numeric variables:", varnames, "\n")

#If the text printed makes sense, run it





#################################### RUN THE ANALYSIS ############################################


for (i in colgroup) {
	print (coltags[i])
	cat("\n")
	cat ("\n\n\t-------\n\t", coltags[i], "\n\t-------\n\n", file="anovas_out.txt", sep=" ", append=TRUE)
	for (j in colvariab) {
		print (coltags[j])
		#Store the data to use
		tempdata <- measurespod[c(i,j)]
		tvariable <- tempdata[[2]]
		tgroup <- tempdata[[1]]
		groupnames <- levels(tgroup)
		numgroups <- length(groupnames)
		
		
		minimum <- min(tvariable)
		maximum <- max(tvariable)
		rangeval = maximum - minimum
		dotsize = (rangeval/100)*2
		dotsize
		
		
		cat ("Analysing", numgroups, "groups:", groupnames, "\n")
		nametest <- paste(coltags[i], "_x_", coltags[j], sep="")
		nameplot <- paste(coltags[i], "_x_", coltags[j], ".png", sep="")
		cat("\n",nametest,"\n")
		
		if (numgroups>2)	{
		
		#do the ANalysis Of VAriance
		aovpod = aov(tvariable~tgroup, data=tempdata)
		print(summary (aovpod))
		cat("\n\n")
		
		
		#save results to file
		out <- capture.output(summary(aovpod))
		cat("\n", "\t", "-- ", coltags[j], " --", "\n\n", file="anovas_out.txt", sep="", append=TRUE)
		cat (nametest, "\n", file="anovas_out.txt", sep=" ", append=TRUE)
		cat(out, file="anovas_out.txt", sep="\n", append=TRUE)
		cat(" \n------", " "," ", file="anovas_out.txt",sep="\n", append=TRUE)
		out = ""														#Empty for the next iteration
		
		
		#PLOTTING
		if (numgroups==3)	{
			palette (c("red1", "goldenrod1", "purple3", "black", "green3", "blue", "magenta", "gray"))
			mypalette <- c("red1", "goldenrod1", "purple3", "black", "green3", "blue", "magenta", "gray")
		} else if (numgroups>3) {	
			palette (c("darkorange3", "darkgreen", "orange", "chartreuse3", "red1", "goldenrod1", "purple3", "black", "gray", "blue", "magenta"))
			mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3", "red1", "goldenrod1", "purple3", "black", "gray", "blue", "magenta")
		}
		
		scale_fill_manual(values=mypalette)
		scale_colour_manual(values=mypalette)
		
		
		png(filename=nameplot, units="in", width=5, height=5, res=300)						#Save to file
		print(		ggboxplot (measurespod, x=coltags[i], y=coltags[j], color=coltags[i], legend="none", size=1, outlier.size=3.5, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize , alpha=0.5)) +
		#geom_dotplot(aes(y=measurespod[j], x =measurespod[i]), binaxis = "y", binwidth = 0.014, stackdir = "center") +
		geom_hline (yintercept=mean (tvariable), linetype = 2) +
		stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
		scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)+
		
		#Report ANOVA
		stat_compare_means (method="anova", label.y.npc="bottom", label.x.npc="centre")+
		stat_compare_means (label="p.signif", method="t.test", ref.group=".all.", label.y.npc="top")		
		)			#Compare each group with the rest and report significance
		dev.off()
		
		
		#PAIRWISE COMPARISON
		bypairs <- pairwise.t.test(tvariable, tgroup, p.adjust.method = "bonferroni")
		bypairs
		
		#save
		out <- capture.output(bypairs)
		cat(out, file="anovas_out.txt", sep="\n", append=TRUE)
		cat("\n--------------------", " "," ", file="anovas_out.txt",sep="\n", append=TRUE)
		out = ""														#Empty for the next iteration
		
		
		
		
		} else if (numgroups == 2)	{
		
		#t-test
		ttest = t.test(tvariable~tgroup, data=tempdata)
		print(ttest)
		cat("\n\n")
		
		#save results to file
		out <- capture.output(ttest)
		cat("\n", "\t", "-- ", coltags[j], " --", "\n\n", file="anovas_out.txt", sep="", append=TRUE)
		cat (nametest, "\n", file="anovas_out.txt", sep=" ", append=TRUE)
		cat(out, file="anovas_out.txt", sep="\n", append=TRUE)
		cat(" \n--------------------", " "," ", file="anovas_out.txt",sep="\n", append=TRUE)
		out = ""														#Empty for the next iteration
		
		
		#PLOTTING
		palette (c("black", "red1", "goldenrod1", "purple3", "green3", "blue", "magenta", "gray"))
		mypalette <- c("black", "red1", "goldenrod1", "purple3", "green3", "blue", "magenta", "gray")
		
		scale_fill_manual(values=mypalette)
		scale_colour_manual(values=mypalette)
		
		png(filename=nameplot, units="in", width=5, height=5, res=300)						#Save to file
		print(ggboxplot (tempdata, x=coltags[i], y=coltags[j], size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
		geom_hline (yintercept=mean (tvariable), linetype = 2)+
		stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
		#Report t-test p-value
		stat_compare_means (method="t.test", label.y.npc="bottom", label.x.npc="centre")+
		#stat_compare_means (label="p.signif", method="wilcox.test", ref.group=".all.", label.y.npc="top")+			#Compare each group with the rest and report significance
		geom_signif(comparisons = list(c(groupnames[1], groupnames[2])), test="t.test", map_signif_level=TRUE, tip_length= 0)+
		scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette))
		dev.off()
		
		}
		
	}
	cat("\n\n\n\n")
	cat ("\n***********************************************************************\n\n\n", file="anovas_out.txt", sep=" ", append=TRUE)

}


beep(3)





########################################################END

















#CUSTOMIZED PLOTS
# run manually line by line, adjust values to your data


palette (c("darkorange3", "darkgreen", "orange", "chartreuse3", "goldenrod1", "green3", "purple3", "black", "gray", "blue", "magenta"))
mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3", "goldenrod1", "green3", "purple3", "black", "gray", "blue", "magenta")
scale_fill_manual(values=mypalette)
scale_colour_manual(values=mypalette)

str(measurespod)

i=5		#position of column with grouping variable you want to plot
j=12		#position of column with quantitative Normal variable you want to analyze


tempdata <- measurespod[c(i,j)]
tvariable <- tempdata[[2]]
tgroup <- tempdata[[1]]
tempdata <- measurespod[c(i,j)]
tvariable <- tempdata[[2]]
tgroup <- tempdata[[1]]
groupnames <- levels(tgroup)
Indiv <- measurespod[[1]]
nameplot <- paste("Groovy_", coltags[i], "_x_", coltags[j], ".png", sep="")
group <- factor(measurespod[[i]])


minimum <- min(tvariable)
		maximum <- max(tvariable)
		rangeval = maximum - minimum
		dotsize = (rangeval/100)*2
		dotsize


png(filename=nameplot, units="in", width=5, height=5, res=300)						#Save to file
print(

# you'll need to manually adjust the position of the comparative bars (geom_signif(y_position)) to be placed inmediatly above (the first two) and below the datapoints
# also you'll need to adjust the size of the dots to the magnitude or the measures (ggboxplot(add.params = list(binwidth)))
ggboxplot (tempdata, x=coltags[i], y=coltags[j], size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
geom_hline (yintercept=mean (tvariable), linetype = 2) +
stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
stat_compare_means (method="anova", label.y.npc="centre", label.x.npc="right") +
geom_signif(comparisons = list(c(groupnames[1], groupnames[3])), y_position=0.298, test="t.test", map_signif_level=TRUE, tip_length= 0.02) +
geom_signif(comparisons = list(c(groupnames[2], groupnames[4])), y_position=0.345, test="t.test", map_signif_level=TRUE, tip_length= 0.02) +
geom_signif(comparisons = list(c(groupnames[1], groupnames[2]), c(groupnames[3], groupnames[4])), y_position=-0.298, test="t.test", map_signif_level=TRUE, tip_length=0, vjust=2) +
geom_signif(comparisons = list(c(groupnames[1], groupnames[4])), y_position=-0.325, test="t.test", map_signif_level=TRUE, tip_length= -0.02, vjust=3)+
scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)

)
dev.off()
	
		
		
		
		
		
		
		
		
		
		

########		END			##########

