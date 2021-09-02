#Change working directorysetwd
#Plot variables, test for normality while aplying multiple transformations

# =====================================================================
# For compatibility with Rscript.exe (only if not loading libraries correctly): 
# =====================================================================
if(length(.libPaths()) == 1){
    # We're in Rscript.exe
    possible_lib_paths <- file.path(Sys.getenv(c('USERPROFILE','R_USER')),
                                    "R","win-library",
                                    paste(R.version$major,
                                             substr(R.version$minor,1,1),
                                             sep='.'))
    indx <- which(file.exists(possible_lib_paths))
    if(length(indx)){
       .libPaths(possible_lib_paths[indx[1]])
    }
    # CLEAN UP
    rm(indx,possible_lib_paths)
} 



#install.packages("ggplot2")
library(ggplot2)
#install.packages("ggfortify")
library(ggfortify)
#install.packages("ggpubr")
library(ggpubr)

#install.packages("ggpmisc")
library(ggpmisc)
#install.packages("npsurv")		# not available anymore
#library(npsurv)		# not available anymore


#install.packages("fitdistrplus")
library(fitdistrplus)
#install.packages("broom")
library(broom)
#install.packages("psych")
library(psych)


#trying this packages for standardisation and normalisations
#install.packages("caret")
library("caret")
#install.packages("heatmaply")
library("heatmaply")





# The input file is:
#		tag, indiv, sex,pop,group,measure1,m2,m3,m4,m5 ... n
#		Am001   001, M,  A,  MA,		04,08,15,16,23 ... 42
#		Bm002   002, F,  B,  FB,		02,71,82,81,82 ... 84
#...etc
# Two first columns must be individual or non-usable data
# Columns with grouping variables (3 in the example, from 3 to 5)
# Columns with your variables
# (In the example tabs and spaces are there only for the data to match the position of the header)



### TRY TO LOOP THROUGH ALL YOUR FILES AND VARIABLES

rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets

#Check working directory
getwd()

setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/nat_pop_Iva/pheno-env/jun21_RDAnew/residuals_norm/norm/bysex")

##		FILE LIST, INPUT PATH AND PATTERN
filelist <- list.files(path = ".", pattern = ".csv$")		  #Change "." for the path of a directory with your data files "folder/subdirectory/here"; You can also make pattern more specific to filter out some files: "_final_version.csv"
filelist		  #check that all the fileas are alright

#open one file to check
wrk <- read.delim(filelist[[1]], sep = ',', header =TRUE)		  #comma separated file with column headers
str(wrk)		  #check structure of the file
head(wrk)		  #check file


#Individual tags if any must be in the left most column
#Number of Columns with Stuff (to ignore) -including individual tags and factor/grouping variables, before the first numerical variable to analyze-
otherstuff=5

# CHECK A FILE
coltags <- names(wrk)
maxcol <- ncol (wrk)								#number of columns
lastgrup = otherstuff 							#column with the last grouping variable
firstvar = lastgrup + 1										# column with the first variable
allvar = maxcol - lastgrup
varnames = coltags[c(firstvar:maxcol)]
varnum = length(varnames)


cat ("\n\nFile checked has", varnum, "variables:", varnames, "\n will be tested for Normality\n\n")







################## RUN EVERYTHING BELLOW IF IT SEEMS ALRIGHT ########################

#warnings()

rawtime = Sys.time()
hour<- sub(".* ", "", rawtime)
onlydate<-sub(" .*", "", rawtime)
printime <- gsub(":", "-", hour)

numfiles = length(filelist)		  # extract how many files there are
myrange = c(1:numfiles)			  #range from 1 to number of files

file1 = paste("Normality_", onlydate, "_", printime, "_out.txt", sep="")
cat ( "Shapiro test of Normality p-values\n", "File_name\tvariable\ttransformation\tvalue_added\tShapiroW pval\tNormality?", sep="\n", file=file1, append=FALSE)

#k=1

#generate an empty data frame for the summary table
summarytable = data.frame(matrix("NA", ncol = 3, nrow = 1))
names(summarytable) <- c("File", "Variable", "Normal transformations")
rown=1

for (k in myrange) {
	csvfile <- filelist[[k]]		  #choose one of the file names in each loop
	# OPEN INPUT FILE
	cat("Processing ", csvfile, "\n")		  #Print on screen
	measurespod <- read.delim(filelist[[k]], sep = ',', header =TRUE)
	
	refrow = c()
	coltags <- names(measurespod)
	#coltags
	maxcol <- ncol (measurespod)								#number of columns
	lastgrup = otherstuff							#column with the last grouping variable
	firstvar = lastgrup + 1										# column with the first variable
	allvar = maxcol - lastgrup
	varnames = coltags[c(firstvar:maxcol)]
	#varnames
	varrange = c(firstvar:maxcol)
	#varrange
	varnum = length(varnames)
	
	
	#extract the variables to test

	#i=15
	for (i in varrange) {
		
		normlist="NONE"
		varname <- coltags[i]
		cat(varname, sep="\n")
		
		
		
		#extract the variable
		variable <- as.vector(measurespod[[i]])
		#variable
		
		#if there is negative values calculate how much you should add to all of them for some transformations
		minim <- min(variable)
		minim
		
		if(minim < -0.1) {
			toadd=0.1
			negmin=abs(minim)+toadd		   #minimum negative value, this will added to all the values, so there is no negative values
			#negmin
		} else if (minim < 0) {
			toadd=0.01
			negmin=abs(minim)+toadd		   #minimum negative value, this will added to all the values, so there is no negative values
			#negmin
		} else { negmin = 0 }
		
		#negmin
		#toadd
		
		
		
		#variable
		
		
		##### SHAPIRO-WILK
		
		#raw
		#perform test
		testout <- shapiro.test (((variable)))
		#declare if something was added to the variable
		added=0
		
		#save if it is significantly different from Normal
		pval = round(testout$p.value, digits=6)
		if (pval < 0.05) {
			isnormal="NO"
		} else {
			isnormal="YES!"
			normlist="no-transf"
		}
		
		#isnormal
		#declare which transformation  was used with the data
		transformed = "no-transf"
		Nans = "no"
		#print table
		cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		
		plotname = paste(csvfile, "_", varname, "-distrb-plot.png", sep="")
		
		png(plotname, units="in", width=10, height=7, res=300)
		descdist(variable, boot=1000)
		dev.off()
		
		
		#log
		transformed = "nat log"
		added=negmin
		#check if there are missing
		checktrans = any(is.na(log(variable+negmin)))
		if(checktrans==FALSE) {
			testout <- shapiro.test ((log(variable+negmin)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		#log10
		transformed = "log 10"
		added=negmin
		#check if there are missing
		checktrans = any(is.na(log10(variable+negmin)))
		if(checktrans==FALSE) {
			testout <- shapiro.test ((log10(variable+negmin)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else { 
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		#squareroot
		added=negmin
		transformed = "sqr root"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^(1/2))))
		if(checktrans==FALSE) {
			testout <- shapiro.test (((variable+negmin)^(1/2)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		#cubic root
		added=negmin
		transformed = "cubic root"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^(1/3))))
		if(checktrans==FALSE) {
			testout <- shapiro.test (((variable+negmin)^(1/3)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		#4th root
		added=negmin
		transformed = "4th root"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^(1/4))))
		if(checktrans==FALSE) {
			checktrans = any(is.na(((variable+negmin)^(1/4))))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		#fifth root
		added=negmin
		transformed = "5th root"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^(1/5))))
		if(checktrans==FALSE) {
			testout <- shapiro.test (((variable+negmin)^(1/5)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		##exponential
		added=0
		transformed = "exp"
		#check if there are missing
		checktrans = any(is.na(exp(variable)))
		checknums = length(unique(exp(variable)))
		if (checknums < 2) { checktrans = TRUE }
		if(checktrans==FALSE) {
			testout <- shapiro.test ((exp(variable)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		#logistic
		added=negmin
		transformed = "logistic"
		#check if there are missing
		checktrans = any(is.na(logistic(variable+negmin)))
		checknums = length(unique(logistic(variable)))
		if (checknums < 2) { checktrans = TRUE }
		
		if(checktrans==FALSE) {
			testout <- shapiro.test ((logistic(variable)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#power of two
		added=negmin
		transformed = "pwr2"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^2)))
		if(checktrans==FALSE) {
			testout <- shapiro.test (((variable)^2))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		
		#power of three
		added=negmin
		transformed = "pwr3"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^3)))
		if(checktrans==FALSE) {
			testout <- shapiro.test (((variable)^3))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		
		##power of four
		added=negmin
		transformed = "pwr4"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^4)))
		if(checktrans==FALSE) {
			testout <- shapiro.test (((variable)^4))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		
		#power of five
		added=negmin
		transformed = "pwr5"
		#check if there are missing
		checktrans = any(is.na(((variable+negmin)^5)))
		if(checktrans==FALSE) {
			testout <- shapiro.test (((variable)^5))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#gamma
		added=0
		transformed = "gamma"
		#check if there are missing
		checktrans = any(is.na(gamma(variable)))
		checknums = length(unique(gamma(variable)))
		if (checknums < 2) { checktrans = TRUE }
		if(checktrans==FALSE) {
			testout <- shapiro.test ((gamma(variable)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#sin
		added=0
		transformed = "sin"
		#check if there are missing
		checktrans = any(is.na(sin(variable)))
		if(checktrans==FALSE) {
			testout <- shapiro.test ((sin(variable)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#asin
		added=0
		transformed = "asin"
		#check if there are missing
		checktrans = any(is.na(asin(variable)))
		if(checktrans==FALSE) {
			testout <- shapiro.test ((asin(variable)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#cos
		added=0
		transformed = "cos"
		#check if there are missing
		checktrans = any(is.na(cos(variable)))
		if(checktrans==FALSE) {
			testout <- shapiro.test ((cos(variable)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#tan
		added=0
		transformed = "tan"
		#check if there are missing
		checktrans = any(is.na(tan(variable)))
		if(checktrans==FALSE) {
			testout <- shapiro.test ((tan(variable)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#scaled sin
		added=0
		transformed = "std sin"
		#check if there are missing
		checktrans = any(is.na(sin(scale((variable), center=TRUE, scale=TRUE))))
		if(checktrans==FALSE) {
			testout <- shapiro.test (sin(scale((variable), center=TRUE, scale=TRUE)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		#scaled cos
		added=0
		transformed = "std asin"
		#check if there are missing
		checktrans = any(is.na(asin(scale((variable), center=TRUE, scale=TRUE))))
		if(checktrans==FALSE) {
			testout <- shapiro.test (asin(scale((variable), center=TRUE, scale=TRUE)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		#scaled cos
		added=0
		transformed = "std cos"
		#check if there are missing
		checktrans = any(is.na(cos(scale((variable), center=TRUE, scale=TRUE))))
		if(checktrans==FALSE) {
			testout <- shapiro.test (cos(scale((variable), center=TRUE, scale=TRUE)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		
		
		
		
		#scaled tan
		added=0
		transformed = "std tan"
		#check if there are missing
		checktrans = any(is.na(tan(scale((variable), center=TRUE, scale=TRUE))))
		if(checktrans==FALSE) {
			testout <- shapiro.test (tan(scale((variable), center=TRUE, scale=TRUE)))
			checkshapir = any(is.na(testout$p.value))
			checkshapir
			pval = round(testout$p.value, digits=6)
			if (checkshapir == TRUE) {
				cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE)
			} else if (pval < 0.05) {
				isnormal="NO"
			} else {
				isnormal="YES!"
				if (normlist=="NONE") { normlist = transformed } else { normlist = c(normlist, transformed) }
			}
			#isnormal
			cat(paste(csvfile, varname, transformed, added, pval, isnormal, sep="\t"), sep="\n", file=file1, append=TRUE)
		} else { cat(paste(csvfile, varname, transformed, added, "#N/A!", "#N/A!", sep="\t"), sep="\n", file=file1, append=TRUE) }
		
		
		printlist = paste(normlist, collapse=" / ")
		addvector=c(csvfile, varname, printlist)
		summarytable[rown,] <- addvector
		
		rown=rown+1
	}
	
	
	
	cat ("\nAll variables in  ", csvfile, "  analyzed.\n\n", sep ="")
	
}

summaryfile = paste("SummaryTable", onlydate, "_", printime, "_Normality.txt", sep="")
write.table(summarytable, file=summaryfile, sep="\t", quote=FALSE, row.names=FALSE)


cat ("\n\nALL FILES ANALYZED!\n\n\n")






