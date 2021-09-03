# R sucks!

#	Change log: 2019-02-26 - M'Ósky
#Implemented to work with multiple files at once
#Implemented Duncan and Tukey tests
#Added comparison bars when 4 groups, and display of letters when multiple groups to show which ones are significantly differents
#There are new restrictions in the input file: don't use "-" in any categorical variable you want to use in the analysis.


##########  USAGE:
# Load the packages, set the input file, check the input file and run the script
# The script will do an ANOVA for all grouping variables with more than 2 factors, and a t-test if only two and between each pair of factors.
# the results and plots will be automatically saved, the plots display the box (median and quartiles), a line for the total average (mean) and the outliers.
# the p-value of ANOVA is displayed if more than two factors
# If there is significant differences (t-test) between one of the groups (factors) and the rest, "*" will appear in the top of the box, otherwise "ns"
# A t-test is also performed for each pair of factors in each grouping variable. If there are 4 factors it will be ploted showing the significance in bars.
# If ANOVA was significant, a post-hoc Tukey test will be applied pair-wise to check from which groups de differences come from.
# If ANOVA was not significant then a Duncan test is used pair-wise.
# If any Tukey or Duncan pair-wise comparison is significant, a box plot showing it will be ploted
# Plots of every test will be saved showing the significance. All Test results will be saved for each input file.
# Also a Multiple-Way ANOVA can be set with two of the grouping variables
##################






###########   INPUT FILE FORMAT: comma separated csv
#
# IMPORTANT: Your categorical variables (used for grouping) can not include "-" in any category/factor
#
# The first row can be a title o comment preceded by "#". Optional
# The Second row contains the name/tags of the columns. If you don't you need to change the option "headers" to "FALSE" whern importing the input file.
#
# First columns: individual tags or non-usable data. In the example there is two ("indiv" and "label_info", from column 1 to 2) but there can be as many as you want.
#   further below you will need to set this as "otherstuff". In this example -> otherstuff = 2
# Next columns: Categorical/grouping explanatory variables (factor). In the example there is three (sex, pop and group, from column 3 to 5) but there can be as many as you want.
#   further below you'll need to set how many are in "groupvar". In this example -> groupvar = 3
# Rest of columns: Quantitative response variables, all of them will be red if they are at the right of the previously defined columns. Add as may as you want.
#   It does not matter if there is a different number of quantitative response variables to analyse in the different files
#   The script will do the tests and plots for every variable found

#Example:
#
#   	#Datafile of phenotypes 1
#   	indiv,  label_info,  sex, pop, group, measure1, m2, m3, m4, m5,   mn
#   	001,   PMf2018_001,  fem,  PM,   pmf,       04, 08, 15, 16, 23,   42
#   	002,   PKm2016_002,  mal,  PK,   pkm,       13, 12, 24, 42, 15,   18,
#   	etc...
#
# Don't mind the space separation in the example, it is there only so the data to mach the headers, the point here is the sorting of columns according to type of variables and the comma separated format.
# Once you are sure everything is alright select everything bellow "RUN THE ANALYSIS" and run it. It may take a while.
##################################################






	#PACKAGES NEEDED
#install.packages("beepr")
library(beepr)
#install.packages("ggplot2")
#install.packages("ggfortify")
library(ggplot2)
library(ggfortify)
#install.packages("ggpubr")
#install.packages("ggsignif")
library(ggpubr)
library(ggsignif)

#install.packages("DescTools")
library(DescTools)
#install.packages("rcompanion")
library(rcompanion)
#install.packages("svglite")
library(svglite)









# CLEAN SCREEN AND MEMORY

rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets









# SETTINGS:


# 1) INPUT FILE(S)


# DIRECTORY

getwd() #check
#Change working directory if needed
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/AdultsPKPM/CalculationStudents/LUKA/AOV")

# Save a list with all the csv files in the working directory
#filelist <- list.files(path = ".", pattern = ".csv")		  #Change "." for the path of a directory with your data files "folder/subdirectory/here"
filelist <- list.files(path = ".", pattern = "*.csv")		  #pattern = ".csv" to read all the files with the same format
# You can also make pattern more specific to filter out some files: pattern= "_final.csv"
# Or directly write the name of a single input file, so it only reads that one: pattern="my_input_file.csv"
filelist   			#Check



# DATA

# Define how many variables of which kind there are in the columns in ALL YOUR INPUT FILES
# Check above in the "Input file format" the order the columns must have according the their type

# OPTIONAL: Open one of the input files to check the variables, rows and columns:
checkfile = read.table(filelist[1], sep = ',', header =TRUE)   			#read the data into a table, read headers
str(checkfile)











# 2) ANALYSIS SETTINGS


# VARIABLES

#Set number of variables. Those two parameters need to be the same for all files.
otherstuff = 1   		#Number of individual or non-relevant variables. Number of columns to the rigth.
groupvar = 3   			#Number of categorical groupin variables. Number of columns before the numerical variables to analyse.
# If this is not the same for all your input files, add dummy variables to the ones with less variables or analyse the one by one.
# Names of the variables dont need to be the same, just their type.



# P-VALUE ANALYSIS

#set which value is the critical limit, any p-value bellow it will be considered significant (before applying any FDR correction.
pcritsignif=0.05000



# PLOT TITLE

#Write a text to add the beginning of the title from the plots with the raw data. Your kind of data is a good choice
titlehead = "Residual distribution of"
#titlehead = 0   			#set to 0 if you do not want anything added to the beggining of the title. Worry not, the 0 will not appear in the title

#Write a text to add at the end of the title from the plots with the raw data. Your species is/are a good choice
titletail = "in 2018 'P. siculus' adults"
#titletail = 0   			#set to 0 if you do not want anything added to the end of the title. Worry not, the 0 will not appear in the title

# It will look something like: "titlehead" "name of variable analysed" per "name of the grouping variable used" "titletail"
cat("Example of a generic plot title:\n", titlehead, names(checkfile[otherstuff+groupvar+1]), "per", names(checkfile[otherstuff+groupvar]), titletail, "\n", sep=' ')











# For the sake of automatic personalization this includes lots of loops inside loops, so it will take a while
#############################################################################################################################
##################################   			RUN THE SCRIPT!!			   ##############################################
#############################################################################################################################



numfiles = length(filelist)		  #how many files there are
myrange = c(1:numfiles)			  #range from 1 to nmber of files

#k=1
#i=4
#j=6

for (k in myrange) {
	csvfile <- filelist[[k]]		  #save a different file name in each loop
	namefile <- gsub(".csv", "", csvfile)		  #delete the extension (.csv) from file
	outfile =paste("param_tests-", namefile, ".txt", sep='')		  #generate a file name including the old file name
	
	
	measurespod = read.table(csvfile, sep = ',', header =TRUE)		  #open a file

	coltags <- names(measurespod)
	maxcol <- ncol (measurespod)								#number of columns
	lastgrup = otherstuff + groupvar							#column with the last grouping variable
	firstvar = lastgrup + 1										# column with the first variable
	firstgrup = otherstuff+1
	numvar = maxcol - lastgrup
	varinfo <- measurespod[c(firstvar:maxcol)]									#take a subset including the data of the measured variables
	groupinfo <- measurespod[c(firstgrup:lastgrup)]									#take a subset with the grouping variables to test
	groupsfound <- names (groupinfo)
	colgroup = c(firstgrup:lastgrup)
	colvariab = c(firstvar:maxcol)
	varnames <- names(varinfo)

	#Print info
	cat ("Analysing file", csvfile, ":\n\t", maxcol, "columns in TOTAL\n\t", otherstuff, "with individual/non-relevant data\n\t", groupvar, " columns with grouping variables:", groupsfound, "\n\t", numvar, "numeric variables:", varnames, "\n")


	#################################### RUN THE ANALYSIS ############################################


	for (i in colgroup) {
		print (coltags[i])
		cat("\n")
		cat ("\n\n\t-------\n\t", coltags[i], "\n\t-------\n\n", file=outfile, sep=" ", append=TRUE)
		for (j in colvariab) {
			print (coltags[j])
			#Store the data to use
			tempdata <- measurespod[c(i,j)]
			tvariable <- tempdata[[2]]
			tgroup <- tempdata[[1]]
			groupnames <- levels(tgroup)
			numgroups <- length(groupnames)
			
			# Calculate the size of the dots according to the range of values and the position of the tags
			minimum <- min(tvariable)
			maximum <- max(tvariable)
			rangeval = maximum - minimum
			dotsize = (rangeval/100)*2
			unit <- rangeval/100
			Qart <- quantile(tvariable)
			Quartiles <- as.data.frame(Qart)
			
			q1 <- Quartiles[1,] - unit
			q5 <- Quartiles[5,] + unit
			
			
			nametest <- paste(coltags[i], "_x_", coltags[j], sep="")
			cat("\n",nametest,"\n")
			
			if (numgroups>2)	{
			
				#ANALYSIS OF VARIANCE
				aovpod = aov(tvariable~tgroup, data=tempdata)
				print(summary(aovpod))
				cat("\n\n")
				summaraov <- unlist(summary(aovpod))
				pval = summaraov["Pr(>F)1"]
				
				#save results to file
				out <- capture.output(summary(aovpod))
				cat("\n", "\t", "-- ", coltags[j], " --", "\n\n", file=outfile, sep="", append=TRUE)
				cat (nametest, "\n\n\tANOVA\n", file=outfile, sep=" ", append=TRUE)
				cat(out, file=outfile, sep="\n", append=TRUE)
				cat(" \n------", " "," ", file=outfile,sep="\n", append=TRUE)
				out = ""														#Empty for the next iteration
				
				
				
				
				#PLOTTING
				if (numgroups==3)	{
					palette (c("red1", "goldenrod1", "purple3", "black", "green3", "blue", "magenta", "gray"))
					mypalette <- c("red1", "goldenrod1", "purple3", "black", "green3", "blue", "magenta", "gray")
				} else if (numgroups==4) {	
					palette (c("darkorange3", "darkgreen", "orange", "chartreuse3", "red1", "goldenrod1", "purple3", "black"))
					mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3", "red1", "goldenrod1", "purple3", "black")
				} else if (numgroups<8) {	
					palette (c("red1", "goldenrod1", "purple3", "olivedrab4", "black", "firebrick4", "slateblue4"))
					mypalette <- c("red1", "goldenrod1", "purple3", "olivedrab4", "black", "firebrick4", "slateblue4")
				} else if (numgroups>=8) {	
					palette (c("black", "slateblue4", "purple2", "magenta", "red1", "firebrick4", "darkorange", "goldenrod1", "olivedrab4", "chartreuse3", "turquoise4", "blue", "maroon4", "darkorange3", "darkgreen", "gray"))
					mypalette <- c("black", "slateblue4", "purple2", "magenta", "red1", "firebrick4", "darkorange", "goldenrod1", "olivedrab4", "chartreuse3", "turquoise4", "blue", "maroon4", "darkorange3", "darkgreen", "gray")
				}
				
				scale_fill_manual(values=mypalette)
				scale_colour_manual(values=mypalette)
				
				
				
				
				
				#PAIRWISE COMPARISON
				
				#PostHoc
				if (pval < pcritsignif)	{
					pairtest="Tukey"
					cat(pairtest, " multiple comparison of means", "\n", sep='')
					bypairs1 <- TukeyHSD(aovpod)
					poshoc <-capture.output( TukeyHSD(aovpod))
					cat(poshoc, file=outfile, sep="\n", append=TRUE)
					cat(" \n------", " "," ", file=outfile,sep="\n", append=TRUE)
					bypairs1$tgroup
					bypairs1 <- as.data.frame((bypairs1$tgroup))
					DTpairsOut <- data.frame(Comparison=factor(row.names(bypairs1)), differences=as.double(bypairs1$diff), P.adj=as.double(bypairs1$"p adj"))   			#Save the information in a dataframe
				} else if (pval >= pcritsignif) {
					pairtest = "Duncan"
					cat(pairtest, "'s new multiple range test", "\n", sep='')
					bypairs1 <- PostHocTest(aovpod, method = "duncan")
					poshoc <-capture.output(PostHocTest(aovpod, method="duncan"))
					cat(poshoc, file=outfile, sep="\n", append=TRUE)
					cat(" \n------", " "," ", file=outfile,sep="\n", append=TRUE)
					bypairs1$tgroup
					bypairs1 <- as.data.frame((bypairs1$tgroup))
					DTpairsOut <- data.frame(Comparison=factor(row.names(bypairs1)), differences=as.double(bypairs1$diff), P.adj=as.double(bypairs1$pval))   			#Save the information in a dataframe
				}
				
				#t.test#
				bypairs2 <- pairwise.t.test(tvariable, tgroup, p.adjust.method = "bonferroni")
				#save
				out2 <- capture.output(bypairs2)
				
				cat("\n--------------------", " "," ", file=outfile,sep="\n", append=TRUE)
				cat(out2, file=outfile, sep="\n", append=TRUE)
				cat("\n--------------------\n\n\n", " "," ", file=outfile,sep="\n", append=TRUE)
				out2 = ""														#Empty for the next iteration
				
				
				#define significantly different groups with letters
				
				if (any(DTpairsOut$P.adj < pcritsignif)) {
					completters <- cldList(P.adj ~ Comparison, data = DTpairsOut, threshold = pcritsignif)
					
					#Plot
					Labels <- as.vector(completters$Letter)    #extract only the letters
				
					X = 1:numgroups
					Y = q5+(unit*2)
					
					nameplot3 <- paste(namefile, "-", coltags[i], "_x_", coltags[j], "-PosHoc", pairtest, ".png", sep="")
					png(filename=nameplot3, units="in", width=5.5, height=5, res=500)						#Save to file
					print(
						ggboxplot (tempdata, x=coltags[i], y=coltags[j], title=paste(pairtest, "pair-wise differences of", coltags[j], "per", coltags[i]), size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
						theme(plot.title = element_text(hjust = 0.5, size=12))+
						geom_hline (yintercept=mean (tvariable), linetype = 2) +
						stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
						scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)+
						stat_compare_means (method="anova", label.y.npc="bottom", label.x.npc="centre")+
						annotate("text", x = X, y = Y, label = Labels)
					)
					dev.off()
					
					
					nameplot35 <- paste(namefile, "-", coltags[i], "_x_", coltags[j], "-PostHoc", pairtest, ".svg", sep="")
					image3=ggboxplot (tempdata, x=coltags[i], y=coltags[j], title=paste(pairtest, "pair-wise differences of", coltags[j], "per", coltags[i]), size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
						theme(plot.title = element_text(hjust = 0.5, size=12))+
						geom_hline (yintercept=mean (tvariable), linetype = 2) +
						stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
						scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)+
						stat_compare_means (method="anova", label.y.npc="bottom", label.x.npc="centre")+
						annotate("text", x = X, y = Y, label = Labels)
					ggsave(file=nameplot35, plot=image3, width=5.5, height=5, dpi=350)
					
				} else {
					nameplot3 <- paste(namefile, "-", coltags[i], "_x_", coltags[j], "-ANOVA.png", sep="")
					
					rawtitle = paste(coltags[j], "per", coltags[i])
					
					if (titlehead != 0) { rawtitle = paste(titlehead, rawtitle) }
					if (titletail != 0) { rawtitle = paste(rawtitle, titletail) }
					
					png(filename=nameplot3, units="in", width=6, height=5, res=500)						#Save to file
					print(
						ggboxplot (tempdata, x=coltags[i], y=coltags[j], title=rawtitle, size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
						theme(plot.title = element_text(hjust = 0.5, size=12))+
						geom_hline (yintercept=mean (tvariable), linetype = 2) +
						stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
						stat_compare_means (method="anova", label.y.npc="bottom", label.x.npc="centre")+
						scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)
					)
					dev.off()
					
					
					nameplot35 <- paste(namefile, "-", coltags[i], "_x_", coltags[j], "-Clean.svg", sep="")
					image3=ggboxplot (tempdata, x=coltags[i], y=coltags[j], title=rawtitle, size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
						theme(plot.title = element_text(hjust = 0.5, size=12))+
						geom_hline (yintercept=mean (tvariable), linetype = 2) +
						stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
						stat_compare_means (method="anova", label.y.npc="bottom", label.x.npc="centre")+
						scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)
					ggsave(file=nameplot35, plot=image3, width=5.5, height=5, dpi=350)
					
					
					
				}
				
				
				
				if (numgroups==4)	{
					  
					nameplot2 <- paste(namefile, "-", coltags[i], "_x_", coltags[j], "-TtestBars.png", sep="")
					
					png(filename=nameplot2, units="in", width=5.5, height=5, res=500)						#Save to file
					print(
						
						ggboxplot (tempdata, x=coltags[i], y=coltags[j], title=paste("t-test pair-wise differences of", coltags[j], "per", coltags[i]), size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
						theme(plot.title = element_text(hjust = 0.5, size=12))+
						geom_hline (yintercept=mean (tvariable), linetype = 2) +
						stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
						geom_signif(comparisons = list(c(groupnames[1], groupnames[3])), y_position=q5+(unit*2), test="t.test", map_signif_level=TRUE, tip_length= 0.02, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[2], groupnames[4])), y_position=q5+(unit*8.5), test="t.test", map_signif_level=TRUE, tip_length= 0.02, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[1], groupnames[2])), y_position=q1-unit, test="t.test", map_signif_level=TRUE, tip_length=0, vjust=2, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[3], groupnames[4])), y_position=q1-unit, test="t.test", map_signif_level=TRUE, tip_length=0, vjust=2, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[2], groupnames[3])), y_position=q1-(unit*2), test="t.test", map_signif_level=TRUE, tip_length=0, vjust=2, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[1], groupnames[4])), y_position=q1-(unit*8), test="t.test", map_signif_level=TRUE, tip_length= -0.02, vjust=3, textsize=3)+
						scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)
					)
					dev.off()
					
					nameplot25 <- paste(namefile, "-", coltags[i], "_x_", coltags[j], "-TtestBars.svg", sep="")   	#save a name for your plot with your variable names
					
					image=ggboxplot (tempdata, x=coltags[i], y=coltags[j], title=paste("t-test pair-wise differences of", coltags[j], "per", coltags[i]), size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
						theme(plot.title = element_text(hjust = 0.5, size=12)) +
						geom_hline (yintercept=mean (tvariable), linetype = 2) +
						stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
						geom_signif(comparisons = list(c(groupnames[1], groupnames[3])), y_position=q5+(unit*2), test="t.test", map_signif_level=TRUE, tip_length= 0.02, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[2], groupnames[4])), y_position=q5+(unit*8.5), test="t.test", map_signif_level=TRUE, tip_length= 0.02, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[1], groupnames[2])), y_position=q1-unit, test="t.test", map_signif_level=TRUE, tip_length=0, vjust=2, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[3], groupnames[4])), y_position=q1-unit, test="t.test", map_signif_level=TRUE, tip_length=0, vjust=2, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[2], groupnames[3])), y_position=q1-(unit*2), test="t.test", map_signif_level=TRUE, tip_length=0, vjust=2, textsize=3) +
						geom_signif(comparisons = list(c(groupnames[1], groupnames[4])), y_position=q1-(unit*8), test="t.test", map_signif_level=TRUE, tip_length= -0.02, vjust=3, textsize=3)+
						scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)
					
					ggsave(file=nameplot25, plot=image, width=5.5, height=5, dpi=350)
				} 
			
			
			
			
			} else if (numgroups == 2)	{
				
				#T test
				ttestout = t.test(tvariable~tgroup, data=tempdata)
				print(ttestout)
				cat("\n\n")
				
				#save results to file
				out <- capture.output(ttestout)
				cat("\n", "\t", "-- ", coltags[j], " --", "\n\n", file=outfile, sep="", append=TRUE)
				cat (nametest, "\n", file=outfile, sep=" ", append=TRUE)
				cat(out, file=outfile, sep="\n", append=TRUE)
				cat(" \n--------------------", " "," ", file=outfile,sep="\n", append=TRUE)
				out = ""														#Empty for the next iteration
				
				
				#PLOTTING
				palette (c("black", "red1", "goldenrod1", "purple3", "green3", "blue", "magenta", "gray"))
				mypalette <- c("black", "red1", "goldenrod1", "purple3", "green3", "blue", "magenta", "gray")
				
				scale_fill_manual(values=mypalette)
				scale_colour_manual(values=mypalette)
				
				nameplot1 <- paste(namefile, "-", coltags[i], "_x_", coltags[j], "-Ttest.png", sep="")
				png(filename=nameplot1, units="in", width=5, height=5, res=300)						#Save to file
				print(
				ggboxplot (tempdata, x=coltags[i], y=coltags[j], size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
					geom_hline (yintercept=mean (tvariable), linetype = 2)+
					stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
					#Report t-test p-value, choose either a bar with "*" or "ns" or the p-value 
					stat_compare_means (method="t.test", label.y.npc="bottom", label.x.npc="centre")+
					#geom_signif(comparisons = list(c(groupnames[1], groupnames[2])), test="t.test", map_signif_level=TRUE, tip_length= 0)+
					scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)
					)
				
				dev.off()
				
				
				nameplot05 = paste(namefile, "-", coltags[i], "_x_", coltags[j], "-Ttest.svg", sep="")
				
				image2=ggboxplot (tempdata, x=coltags[i], y=coltags[j], size=1, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5), color=coltags[i], legend="none", outlier.size=3.5) +
					geom_hline (yintercept=mean (tvariable), linetype = 2)+
					stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
					#Report t-test p-value
					stat_compare_means (method="t.test", label.y.npc="bottom", label.x.npc="centre")+
					#geom_signif(comparisons = list(c(groupnames[1], groupnames[2])), test="t.test", map_signif_level=TRUE, tip_length= 0)+
					scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)
				
				ggsave(file=nameplot05, plot=image2, width=5, height=5, dpi=350)
				
				
				cat ("\n\n")
			
			} 
		}
		cat("\n\n")
		cat ("\n***********************************************************************\n\n\n", file=outfile, sep=" ", append=TRUE)
		
	}


	cat ("File ", csvfile, " done!\n\n#########################################\n\n\n\n", sep='')

}

warnings()


beep(4)


cat ("\n\nDONE!!\nAll files analysed!\n\n", sep='')






