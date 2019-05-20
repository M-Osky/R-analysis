#R sucks! <-- One day you'll be thankfull that you specify the language you are using in each script
# M'Ósky 2019

# USAGE: This script should automatically analyze all the csv files in the wd
# It will output box plots and a file with the Kruskal-Wallis test for each grouping variable defined and a pari-wise wilcoxon test

#I just wrote it today and I'm going to upload it, so it may need some debugging
#Input file as in my other scripts, columns with individual tags or sueless data, then columns with grouping variables, then columns with the variables to analyse
#observations in rows

getwd()
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/F1crosses/Petra")





#LOAD LIBRARIES
library("ggplot2")
library("ggfortify")
library("ggpubr")
library("ggsignif")
library("beepr")

rm(list = ls())		#Remove all objects



# SETTINGS:

# 1) INPUT FILE

# Save a list with all the csv files in the working directory
filelist <- list.files(path = ".", pattern = ".csv")		  #Change "." for the path of a directory with your data files "folder/subdirectory/here"; You can also make pattern more specific to filter out some files: "_final_version.csv"
filelist

# Define how are the columns in your input file
firstgroup = 2		#First column with a grouping variable
firstvar = 7		#First column with a measure (numerical variable)





################################################################ THIS SHOULD RUN SMOOTHLY

numfiles = length(filelist)		  #how many files there are
myrange = c(1:numfiles)			  #range from 1 to nmber of files

lastgroup = firstvar - 1		  #Calculate the column that contains the last grouping variable by resting one to hte column with the first numerical variable

# LOOP THROUGHOUT ALL THE FILES

#k=1

for (k in myrange) {
	csvfile <- filelist[[k]]		  #save a different file name in each loop
	# OPEN INPUT FILE
	cat("Analysing ", csvfile, "\n")		  #Print on screen
	
	f1measures = read.table(csvfile, sep = ',', header =TRUE)		  #open a file
	str(f1measures)		  #General structure and organization of dataset
	
	maxcol <- ncol (f1measures)		  #Extract the number of columns of your file
	
	coltags=names(f1measures)		  #extract the names of the columns (headers)
	groupcol=c(firstgroup:lastgroup)		  #save a range with the position of the columns with the groups
	variablecol=c(firstvar:maxcol)		  #save a range with the position of the columns with the variables
	
	# Results file
	namefile <- gsub(".csv", "", csvfile)		  #delete the extension (.csv) from file
	outfile =paste("nonparam_test-", namefile, ".txt", sep='')		  #generate a file name including the old file name
	titlefile=paste("Analyzed data from: ", csvfile, sep="")
	cat("DIFFERENCES BETWEEN GROUPS", "--------------------------", "", titlefile, "####################", "\n", sep = '\n', file = outfile)		 #create a output file

#i=3
#j=14

	# LOOP THROUGH ALL THE GROUPS AND VARIABLES IN THE FILE
	
	for (i in groupcol) {
		
		tgroup=f1measures[[i]]		  #extract one group per loop
		tgroup
		
		for (j in variablecol) {
			
			
			
			tvariable=f1measures[[j]]		  #extract one variable per loop
			tvariable
			
			groupnames=levels(tgroup)		  #extract group names (factors)
			numgroups=length(groupnames)		  #extract number of factors (groups)

			# Calculate the size of the dots according to the range of values
			minimum <- min(tvariable)
			maximum <- max(tvariable)
			rangeval = maximum - minimum
			dotsize = (rangeval/100)*2
			dotsize
			
			
			#SAVE A LIST OF COLORS AS A PALETTE FOR PLOTS
			if (numgroups==2)	{
				palette (c("red1", "black", "green3", "orange", "magenta", "gray", "blue"))
				mypalette <- c("red1", "black", "green3", "orange", "magenta", "gray", "blue")
			} else if (numgroups==3)	{
				palette (c("red1", "goldenrod1", "purple3", "black", "green3", "orange", "magenta", "blue"))
				mypalette <- c("red1", "goldenrod1", "purple3", "black", "green3", "orange", "magenta", "blue")
			} else if (numgroups==4) {	
				palette (c("darkorange3", "darkgreen", "orange", "chartreuse3", "black", "red1", "goldenrod1", "purple3", "blue", "magenta"))
				mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3", "black", "red1", "goldenrod1", "purple3", "blue", "magenta")
			} else if (numgroups>4) {	
				palette (c("darkorange3", "darkgreen", "purple3", "goldenrod1", "red1", "black", "orange", "chartreuse3"))
				mypalette <- c("darkorange3", "darkgreen", "purple3", "goldenrod1", "red1", "black", "orange", "chartreuse3")
			}
			scale_fill_manual(values=mypalette)
			scale_colour_manual(values=mypalette)

			if (numgroups==2) {
				
				# Perfor Wilcoxon test to compare both groups
				wmwtest = wilcox.test(tvariable~tgroup, data=f1measures)
				wilcoxout=capture.output(wmwtest)		  #Save output
				cat(coltags[[i]], " x ", coltags[[j]], "\n", file=outfile, sep = ' ', append = TRUE)		  #generate title for the analysis
				cat(wilcoxout, file=outfile, sep = '\n', append = TRUE)		  #print to the file the output
				cat("", "--------------------------------", "\n", sep="\n", file=outfile, append=TRUE)


				# Plotting
				
				plotname=paste(namefile, "_", coltags[[i]],"-", coltags[[j]], ".png", sep = "")
				png(filename=plotname, units="in", width=5, height=5, res=300)      #Save to file
				print(ggboxplot (f1measures, x=coltags[i], y=coltags[j], color=coltags[i], legend="none", size=1, outlier.size=3.5, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5)) +
					geom_hline (yintercept=mean (tvariable), linetype = 2) +
					stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
					scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)+
			
					#Report Wilcoxon p-value
					stat_compare_means (method="wilcox.test", label.y.npc="bottom", label.x.npc="centre"))
				dev.off()
			}
			
			else if (numgroups > 2) {
			
				# Perform Kruskal-Wallis analysis
				kwfirst=kruskal.test(tvariable~tgroup,data=f1measures)
				kwout=capture.output(kwfirst)		  #save output
				cat(coltags[[i]], " x ", coltags[[j]], "\n", file=outfile, sep = ' ', append = TRUE)		  #Print analysis name and variable names to a file
				cat(kwout, file=outfile, sep = '\n', append = TRUE)		  ##save Kruskal Walllis output to a file
				cat("", "--------------------------------", "\n", sep="\n", file=outfile, append=TRUE)


				# Pair-wise comparisons with Wilcoxon test
				pairtest = pairwise.wilcox.test(tvariable, tgroup, p.adjust.method = "bonferroni")		  #perform wilcoxon test for all pairs of groups, adjust the p-value for multiple comparisons with Bonferroni correction
				wilcox_out=capture.output(pairtest)		  #save output
				cat(coltags[[i]], " x ", coltags[[j]],"   Wilcoxon test", "\n", file=outfile, sep = "", append=TRUE)
				cat(wilcox_out, "", "", "######################################################", "\n", "\n", file=outfile, sep = '\n', append=TRUE)


				# Plotting
				
				plotname=paste(namefile, "_", coltags[[i]],"-", coltags[[j]], ".png", sep = "")
				png(filename=plotname, units="in", width=7, height=5, res=300)      #Save to file
				print(ggboxplot (f1measures, x=coltags[i], y=coltags[j], color=coltags[i], legend="none", size=1, outlier.size=3.5, add="dotplot", add.params = list(color = coltags[i], binwidth=dotsize, alpha=0.5)) +
					geom_hline (yintercept=mean (tvariable), linetype = 2) +
					stat_summary(fun.y=mean, geom="point", shape=95, size=6) +
					scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette)+
					
					# Report Kruskal-Wallis p-value
					stat_compare_means (method="kruskal.test", label.y.npc="bottom", label.x.npc="centre")+
					# Report significance of the pair-wise comparisons
					stat_compare_means (label="p.signif", method="wilcox.test", ref.group=".all.", label.y.npc="top"))
					
				dev.off()
			}
		}
	}
	
	cat ("File ", csvfile, " done!\n", sep='')

}




beep(4)


cat ("\n\nALL FILES ANALYSED!\n\n", sep='')




















