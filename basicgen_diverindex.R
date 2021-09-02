#!/usr/bin/Rscript

# DESCRIPTION
# First R script attempt for evoleco
# Use basicgen.R --help for more information
# you need a genepop input file. As genepop doesn't keep the population names,
# if you want proper population names you will need also a Structure input file
# and the population code should be part of the sample ID name in the Structure file
# Both files should have the same name but different extensions: ".str" for Structure and ".gen" for Genepop

#Check the packages needed below and make sure that they are installed, then you can run the script just calling it: basicgen.R -i inputfilename







############################################################################
####      Script, you know... don't touch anything unless... bla bla bla
###########################################################################

#######################   		Parsing command line arguments
cat("\n\n#################################################\n\nStarting: Reading command line arguments...\n\n")
library(optparse)


#> set interface
option_list <- list(

  make_option(c("-f", "--flag"),    action="store_true", default=FALSE,   help="Flag [default: %default]"),
  make_option(c("-i", "--input"),   type="character",    default=FALSE,   help="Genepop input file name without the extension. This name must be the same for the Structure file (if any)"),
  make_option(c("-b", "--bootspval"), type="integer",      default=999, help="MCMC iterations for HWE test sigificance and div stats CIs [default: %default]"),
  make_option(c("-r", "--bootsci"), type="integer",      default=100,     help="bootstrap replicates to generate confidence intervals for Fst [default: %default]"),
  make_option(c("-n", "--namepops"),  type="character",    default="yes", help="Do you want to extract the population names from the sample names of a Structure file? Mandatory to calculate Fst per population! (yes/no) [default: %default]"),
  make_option(c("-p", "--poplength"), type="integer",      default=2,     help="How many characters long is the population ID at the beginning of each sample name? [default: %default]")
)

parser <- OptionParser(usage="%prog [options]\n\nUse this script to calculate basic diversity indexes and population genetics indexes from a Genepop file.\n\nIf you get an \'X11 display\' error it means you need to set a virtual framebuffer X11 server.\nCheck \"Initiate and declare a virtual graphical environment\" in the beginning of the script.\n", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)

opt <- args$options
file <- args$args

#< set interface
#######################   		LIBRARIES

cat("\n\n>Setting virtual display environment and loading packages.\n")

###### Initiate and declare a virtual graphical environment #########
# The line below is very important, if there is not display environment for X11 you'll need a virtual framebuffer X11 server: Xvfb
	#Uncomment the line below the first time you run it, and don't worry if you get a "Fatal server error: Server is already active for display 0" 
	#Unix is very melodramatic but that just means the command was not needed. Line below for virtual framebuffer X11 server:

#system("Xvfb :0 -ac -screen 0 1960x2000x24 &")   		#<-- This line
Sys.setenv("DISPLAY"=":0")


###### R Packages ##########
library(BiocManager)
#install.packages("devtools")   	#multipurpose package with plenty tools and for handling other packages
library("devtools")
cat("\n")
library("adegenet")
##install.packages("dartR")
#install_github("green-striped-gecko/dartR")
##install.packages("dartR", repos="https://mirrors.nic.cz/R/")
library(parallel)
library(dartR)
#install.packages("StAMPP")
library(StAMPP)
library(diveRsity)
library ("hierfstat")   	#package for genetic distances





inputname <- opt$input
proper_popnames <- opt$namepops
poplength <- opt$poplength
bootrep<- opt$bootspval
bootci<- opt$bootsci




######################################   		Setting it from the script if needed
#
#inputname = "psic567_1g2g"   		#withouth the extension 
#inputname = "psic567_1g2g"   		#withouth the extension 

#proper_popnames = "yes"   		#You want to extratc the population codes from samples IDs in a different Structure file? (yes/no)
#proper_popnames = "no"   		#You want to extratc the population codes from samples IDs in a different Structure file? (yes/no)
#proper_popnames = "yes"   		#You want to extratc the population codes from samples IDs in a different Structure file? (yes/no)


# If you want the proper population names in the tables, you will need to set also:
#poplength = 2 #how many characters of the sample name belong to population name
#poplength = 2 #how many characters of the sample name belong to population name



######################################



# Open the input file to add it to Adegenet

genepopfile = paste(inputname, ".gen", sep="")   		#write Genepop file name

cat("Done.\n----------------------------\n", sep="")
# rawgenind <- read.genepop(genepopfile, ncode=3L)   		#open it in adegenet
# rawgenind
# popnames <- unique(pop(rawgenind))
# cat("\nDefault population names are:", popnames, "\n", sep=" ")



#replace popnames
if (proper_popnames=="yes") {
	
	strucfile = paste(inputname, ".str", sep="")  		#write Structure file name

	# Open the input file to check how many SNPs, samples, etc.
	cat("\nOpening Structure file (", strucfile, ") to read population names from the first ", poplength, " characters of each sample name\n", sep="")
	dataset = read.table(strucfile, sep = '\t', header =TRUE)

	# Extract the population list per sample from the sample names, the number of loci, etc
	snp_num <- ncol(dataset)-2   		#number of loci
	samplerows <- length(dataset$X)   		#number of samples
	samplenum = samplerows/2
	odd_indexes<-seq(1,samplerows, 2)   			#list of odd indexes, because there is two rows per sample
	pops <- as.data.frame(factor(dataset[[1]]))   		#list of samples in the file
	newnames = substr(as.vector(pops[odd_indexes,1]), 1, poplength)   		# Extract only odd sample names and take only the first characters of its name
	popnames <- as.vector(unique(newnames))   		#list of populations	
	#pop(rawgenind) <- newnames
	cat("\n\nFound ", samplenum, " individuals and ", snp_num, " SNPs\n\n", sep="" )
	cat("\nDone. New population names are:", popnames, "\n\n Saving rest of  the str file information. This may take some time...\n\n", sep=" ")
	rawgenind <- read.structure(strucfile, row.marknames=1, onerowperind=FALSE, n.ind=samplenum, n.loc=snp_num, col.lab=1, col.pop=2, NA.char = "0", ask=FALSE)
	pop(rawgenind) <- newnames
	rawgenind
	
} else if (proper_popnames=="no") {
	cat ("\nUsing first sample per population as popname\n")
} else {
	cat ("\nERROR: Something wrong with the settings, keeping default genepop pop names.\n")
	q(save = "no", status = 1)
}




############## DIVERSITY

cat("\nNow opening ", genepopfile, " and calculating allelic diversity indexes, this will take long... ", sep="")
divstats1 = paste("RAWdivRstats_", inputname, "_", bootrep, "repl.txt", sep="")

cat("\n\nDone!\n\nNow saving.\n\n")

allel_div <- basicStats(infile = genepopfile, outfile = divstats1, fis_ci = TRUE, ar_ci = TRUE, fis_boots = bootrep, ar_boots = bootrep, mc_reps = bootrep, rarefaction = FALSE, ar_alpha = 0.05,fis_alpha = 0.05)
cat("Done! Processing results.\n\n", sep="")
totals <- lapply(allel_div$main_tab, function(x) (x)$overall)

#Let's assign the right names to the populations

if (proper_popnames=="yes") {
	names(totals) <- popnames
}


#lapply generates a list of results but we need a proper dataframe to be able to plot it with the tags and everything
#transform the output of lapply to a proper dataframe
total.df <- as.data.frame(totals)

#problem, it didn't keep the name of each index (just numbers)
#now we need to identify each value with the index it belongs to
#those names are stored for each population in the output
indexnames <- rownames(allel_div$main_tab[[1]])

#add he column names
rownames(total.df) <- indexnames
indexes <- as.data.frame(t(total.df))   		#translocate
indexes

divstats2 = paste("divRstats_", inputname, "_", bootrep, "repl.txt", sep="")
cat("Diversity indexes (diveRsity): ", "", sep="\n",file=divstats2, append=FALSE)
write.table (format (indexes, digits=4),file=divstats2, sep="\t", append=TRUE, quote=FALSE)




#let's try to plot a simple table with per population Fst

if (proper_popnames=="yes") {

	#basicstats <- basic.stats(rawgenind, diploid = TRUE)
	#basicstats
	#hierfstCI <- betas(rawgenind, nboot=bootrep,lim=c(0.025,0.975),diploid=TRUE,betaijT=FALSE)
	
	cat("Now calculating Fst per population\n")
	hierfstCI <- betas(rawgenind, nboot=bootci,lim=c(0.025,0.975),diploid=TRUE)

	#str(hierfstCI)

	#hierfstCI$betaiovl
	#hierfstCI$ci
	#popnames <- names(hierfstCI$betaiovl)
	
	#cat("\npopnames: ", popnames, "\n\n", sep="")
	#popnames
	avgFstperpop <- as.data.frame(hierfstCI$betaiovl)
	#avgFstperpop

	limitsFstperpop <- as.data.frame(hierfstCI$ci)
	
	
	#cat("\n\ntable\n")
	#cat(capture.output(limitsFstperpop), sep="\n")
	#cat("\n\n")
	#cat(capture.output(str(limitsFstperpop)), sep="\n")
	#cat("\n\n")
	
	#print1 = as.vector(names(limitsFstperpop))
	#print2 = as.vector(colnames(limitsFstperpop))
	#print3 = as.vector(rownames(limitsFstperpop))
	
	#cat("\nReplacing data frame names:", print1, "\ncolumn names:", print2, "\nor row names:", print3, "\nwith population names:", popnames, "\n\n\n", sep=" ")
	
	
	#names(limitsFstperpop) <- popnames
	
	#addednames <- as.vector(names(limitsFstperpop))
	#cat("\nAdded names are now:", addednames, "\n\n\n", sep=" ")
	
	transloc <- as.data.frame(t(limitsFstperpop) )
	
	#aftertrans <- as.vector(rownames(transloc))
	
	#cat("\nRownames (after translocate) are:", aftertrans, "\n\n\n", sep=" ")
	
	#transloc
	#names(transloc)
	names(transloc) <- c("lowerCI", "upperCI")
	#colnames(transloc) <- c("lowerCI", "upperCI")
	#rownames(transloc) <- popnames
	#transloc
	avgFstperpop$lowerCI <- transloc$lowerCI
	avgFstperpop$upperCI <- transloc$upperCI
	avgFstperpop$Pops <- popnames
	names(avgFstperpop) <- c("avgFst", "lowerCI", "upperCI", "Population")
	printfsts <- avgFstperpop[, c("Population", "avgFst", "lowerCI", "upperCI")]
	fstcitable = paste("FstsCIperpop_", inputname, "_", bootci, "repl.txt", sep="")
	#write.table (format (avgFstperpop, digits=5),file=fstcitable, sep="\t", append=FALSE, quote=FALSE, col.names=TRUE, row.names=TRUE)
	write.table (format (printfsts, digits=5),file=fstcitable, sep="\t", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)


}





cat("\nAll results saved!\n")

q(save = "no", status = 0)
