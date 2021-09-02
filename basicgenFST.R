#!/usr/bin/Rscript

# DESCRIPTION
# First R script attempt for evoleco
# Use basicgen.R --help for more information
# you need a genepop input file. As genepop doesn't keep the population names,
# if you want proper population names you will need also a Structure input file
# and the population code should be part of the sample ID name in the Structure file
# Both files should have the same name but different extensions: ".str" for Structure and ".gen" for Genepop

#Check the packages needed below and make sure that they are installed, then you can run the script just calling it: basicgen.R -i inputfilename

#updated 03022021 to output CI


#######################   		LIBRARIES


library(optparse)




############################################################################
####      Script, you know... don't touch anything unless... bla bla bla
###########################################################################

#######################   		Parsing command line arguments
cat("\n\n#################################################\n\nStarting: Reading command line arguments...\n\n")


#> set interface
option_list <- list(

  make_option(c("-f", "--flag"),    action="store_true", default=FALSE,   help="Flag [default: %default]"),
  make_option(c("-i", "--input"),   type="character",    default=FALSE,   help="Genepop input file name without the extension. This name must be the same for the Structure file (if any)"),
  make_option(c("-b", "--bootspval"), type="integer",      default=99999, help="bootstrap iterations to generate p-values for Fsts and MCMC for HWE test [default: %default]"),
  make_option(c("-r", "--bootsci"), type="integer",      default=999,     help="bootstrap replicates to generate confidence intervals for Fis and Ar [default: %default]"),
  make_option(c("-n", "--namepops"),  type="character",    default="yes", help="Do you want to extract the population names from the sample names of a Structure file? (yes/no) [default: %default]"),
  make_option(c("-p", "--poplength"), type="integer",      default=2,     help="How many characters long is the population ID at the beginning of each sample name? [default: %default]")
)

parser <- OptionParser(usage="%prog [options]\n\nUse this script to calculate basic diversity indexes and population genetics indexes from a Genepop file.\n\nIf you get an \'X11 display\' error it means you need to set a virtual framebuffer X11 server.\nCheck \"Initiate and declare a virtual graphical environment\" in the beginning of the script.\n", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)

opt <- args$options
file <- args$args

#< set interface
cat("\n\n>Setting virtual display environment and loading packages.\n")

###### Initiate and declare a virtual graphical environment #########
# The line below is very important, if there is not display environment for X11 you'll need a virtual framebuffer X11 server: Xvfb
	#Uncomment the line below the first time you run it, and don't worry if you get a "Fatal server error: Server is already active for display 0" 
	#Unix is very melodramatic but that just means the command was not needed. Line below for virtual framebuffer X11 server:

system("Xvfb :0 -ac -screen 0 1960x2000x24 &")   		#<-- Uncomment this line if needed
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





# Figure format
fig <- paste (opt$output,"svg", sep = ".") 


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
#rawgenind <- read.structure(strucfile, row.marknames=1, onerowperind=FALSE, n.ind=samplenum, n.loc=snp_num, col.lab=1, col.pop=2, NA.char = "0", ask=FALSE)

genepopfile = paste(inputname, ".gen", sep="")   		#write Genepop file name

cat("Done.\n----------------------------\n\nReading input file (", genepopfile, ") this will take a while..\n", sep="")
rawgenind <- read.genepop(genepopfile, ncode=3L)   		#open it in adegenet
rawgenind
popnames <- unique(pop(rawgenind))
cat("\nDefault population names are:", popnames, "\n", sep=" ")



#replace popnames
if (proper_popnames=="yes") {
	
	strucfile = paste(inputname, ".str", sep="")  		#write Structure file name

	# Open the input file to check how many SNPs, samples, etc.
	cat("\nOpening Structure file (", strucfile, ") to read population names from the first ", poplength, " characters of the sample name\n", sep="")
	dataset = read.table(strucfile, sep = '\t', header =TRUE)

	# Extract the population list per sample from the sample names, the number of loci, etc
	snp_num <- ncol(dataset)-2   		#number of loci
	samplenum <- (length(dataset$X) - 1)   		#number of samples

	odd_indexes<-seq(1,samplenum, 2)   			#list of odd indexes, because there is two rows per sample
	pops <- as.data.frame(factor(dataset[[1]]))   		#list of samples in the file
	newnames = substr(as.vector(pops[odd_indexes,1]), 1, poplength)   		# Extract only odd sample names and take only the first characters of its name
	popnames <- unique(newnames)   		#list of populations	
	pop(rawgenind) <- newnames
	cat("\nDone. New population names are:", popnames, "\n", sep=" ")
	
} else if (proper_popnames=="no") {
	cat ("\nUsing first sample per population as popname\n")
} else {
	cat ("\nERROR: Something wrong with the settings, keeping default genepop pop names.\n")
	q(save = "no", status = 1)
}



# Transform Genind format to Genlight format

cat("\nTransforming data format to genlight...", sep="")
genlux <- gi2gl(rawgenind)   		#genind to genlight
cat(" Done!\n", sep="")

cat(" File:\n", sep="")

genlux

#cat("\n\n\n\nSummary:\n", sep="")

#sum(genlux)

cat("\n\n", sep="")

############## PAIR-WISE FSTs

fstfile = paste("pairFst_", inputname, ".txt", sep="")

# Calculate pair-wise Fsts with p-values
cat("\nCalculating pairwise Fsts, this may take a while... ", sep="")
pairFst <- stamppFst(genlux, nboots = bootrep, percent = 95, nclusters = 1)
cat("Done!\nSaving results.\n\n", sep="")
lowest <- min(as.numeric(unlist(pairFst$Pvalues)), na.rm=TRUE)
highest <- max(as.numeric(unlist(pairFst$Pvalues)), na.rm=TRUE)
fsts <- pairFst$Fsts
pvals <- pairFst$Pvalues


simplefile=paste("FSTmatrix_", inputname, ".txt", sep="")

#save
cat("", "Pair-wise Fsts", "", "Fsts:", sep="\n",file=fstfile, append=FALSE)
write.table (format (fsts, digits=5),file=fstfile, sep="\t", append=TRUE, quote=FALSE)
write.table (format (fsts, digits=5),file=simplefile, sep="\t", append=FALSE, quote=FALSE)
cat("\n", "p-values: ", lowest, " - ", highest, "\n", sep="",file=fstfile, append=TRUE)
write.table (format (pvals, digits=5),file=fstfile, sep="\t", append=TRUE, quote=FALSE)

#table with CI and rest of Info
lastcol=  ncol(pairFst$Bootstraps)
firstcol= lastcol - 3
allfstdat <- pairFst$Bootstraps[,c(1,2,firstcol:lastcol)]
fsttable=paste("TableAllFst_", inputname, ".txt", sep="")
write.table (format (allfstdat, digits=5),file=fsttable, sep="\t", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)



############ AMOVA

#cat("\n\nNow computing AMOVA\n")
#gendist <- stamppNeisD(genlux, pop = TRUE)
#stamppAmova(gendist, genlux, perm = 100)


############## DIVERSITY

# cat("\nCalculating allelic diversity indexes, this may take a while... ", sep="")
# divstats1 = paste("RAWdivRstats_", inputname, sep="")

# allel_div <- basicStats(infile = genepopfile, outfile = divstats1, fis_ci = TRUE, ar_ci = TRUE, fis_boots = bootci, ar_boots = bootci,mc_reps = bootrep, rarefaction = FALSE, ar_alpha = 0.05,fis_alpha = 0.05)
# cat("Done! Processing results.\n\n", sep="")
# totals <- lapply(allel_div$main_tab, function(x) (x)$overall)

# #Let's assign the right names to the populations

# if (proper_popnames=="yes") {
	# names(totals) <- popnames
# }


# #lapply generates a list of results but we need a proper dataframe to be able to plot it with the tags and everything
# #transform the output of lapply to a proper dataframe
# total.df <- as.data.frame(totals)

# #problem, it didn't keep the name of each index (just numbers)
# #now we need to identify each value with the index it belongs to
# #those names are stored for each population in the output
# indexnames <- rownames(allel_div$main_tab[[1]])

# #add he column names
# rownames(total.df) <- indexnames
# indexes <- as.data.frame(t(total.df))   		#translocate
# indexes

# divstats2 = paste("divRstats_", inputname, ".txt", sep="")
# cat("Diversity indexes (diveRsity): ", "", sep="\n",file=divstats2, append=FALSE)
# write.table (format (indexes, digits=4),file=divstats2, sep="\t", append=TRUE, quote=FALSE)


# cat("\nAll results saved!\n")

# q(save = "no", status = 0)
