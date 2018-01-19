#R script to do simple ANOVA of all the possible combination between 
# 3 different grouping variables and any number of numerical variables


#Change working directorysetwd
setwd("D:/Dropbox/MOSKY/CURRO/DATA/Measurements/Analysis/ANOVA")
getwd()			#Check

library(beepr)

rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces

#INPUT FILE:
# The first row can be a title o comment preceded by "#" 
# Secon row contains the name/tags of the columns
# The first 5 columns from the input file MUST BE meta data about the samples. 
# All the other columns (doesn't matter how many) MUST BE variables to analyze"
# Typically first and second col will be the individual tag, and any other info
# columns 3 to 5 MUST BE grouping data, I used "sex", "population" and "group", but doesn't matter
#Something like this:
#
#		#Datafile of phenotypes 1
#		indiv,	label_info,	 sex, pop, group, measure1, m2, m3, m4, m5 ... mn
#		001,	PMf2018_001, fem, PM,  pmf,			04, 08, 15, 16, 23 ... 42
#		002,	PKm2016_002, mas, PK,  pkm,			13, 12, 24, 42, 15 ... 18
#		etc...
#
# Don't mind the tabulation, the  tabs are there only for the data to mach the headers


measurespod = read.table("phenodata1.csv", sep = ',', header =TRUE)	#read the data into a table, read headers
coltags <- names(measurespod)
maxcol <- ncol (measurespod)								#number of columns
cat (maxcol, "columns: ", coltags, "\n")
colgroup = c(3:5)
colvariab = c(6:maxcol)





for (i in colgroup) {
	print (coltags[i])
	cat("\n")
	cat ("\n\n\t-------\n\t", coltags[i], "\n\t-------\n\n", file="anovas_out.txt", sep=" ", append=TRUE)
	for (j in colvariab) {
		
		tempdata <- measurespod[c(i,j)]
		tempdata
		tvariable <- tempdata[[2]]
		tgroup <- tempdata[[1]]

		aovpod = aov(tvariable~tgroup, data=tempdata)					#do the analysis of variance
		nametest <- paste(coltags[i], " ~ ", coltags[j], ".png", sep="")
		cat("\n",nametest,"\n")
		print(summary (aovpod))
		cat("\n\n")

		out <- capture.output(summary(aovpod))							#save results to file
		cat (nametest, "\n", file="anovas_out.txt", sep=" ", append=TRUE)
		cat(out, file="anovas_out.txt", sep="\n", append=TRUE)
		cat(" \n--------------------", " "," ", file="anovas_out.txt",sep="\n", append=TRUE)
		out = ""														#Empty for the next iteration

		#png(filename=nametest)											#save plot
		#boxplot(tvariable~tgroup,data=tempdata)						#graphical summary
		#dev.off()
		
	}
	cat("\n\n\n\n")
	cat ("\n***********************************************************************\n\n\n", file="anovas_out.txt", sep=" ", append=TRUE)
}

beep(3)
