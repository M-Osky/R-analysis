#R sucks!
#R script to filter out loci not in HWe in a given dataset 





#######################################################################
########		changelog 14/09/2020 Fixed some error		###########
#at some point the function to delete loci gl.drop.loc
#checked the loci metadata before deleting loci and if there
#was none, it will fail. I deleted that "metadata check" from the
#source code of the function loci gl.drop and it worked fine
#but then I realised that would not be a good solution because it
#will work only for me.
#then I noticed that the only thing the function was checking was 
# to check if the metadata dimensions agreed with the number of loci
# so just added a dataframe with the loci to the metadata and it worked-
###########################################################################






getwd()

#install.packages("adegenet")   	#main genetic package
library ("adegenet")   	#main genetic packages
#install.packages("pegas")
library("pegas")   	# We will use it mainly for testing HWeq

#install.packages("beepr")
library("beepr")
#install.packages("tidyverse")
#library("tidyverse")

#install.packages("BiocManager")
library("BiocManager")

#BiocManager::install("SNPRelate")
library("SNPRelate")

#BiocManager::install("qvalue")
library("qvalue")

#install.packages("dartR")
library("dartR")




setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/PopGenet/genpopR/sicula596x39k/sicula585x39K_noutfix")




rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets



############################################################## PARAMETERS ########

# classic structure input file, replacing soaces with tabs: sed -i 's/ /\t/g' file.str
# don't need to edit population or sample names. further in the script there is some lines to do it

# Input file name (no extension)
inputname = "sicula585x39k_m05R7r6h6DP420_1p0p_0m0p"



# How many characters of the sample name belong to the population code
poplength=3

#threshold (critical) value below which a pvalue is considered significant
refval = 0.05

# number of iterations to run to calculate the p-value for the HWe test
iterations=99

#proportion of populations in which a locus must be out of HWE in order for it to be deleted
ref_outhw = 0.6

#######################################################################



# open and get the information
strucfile = paste(inputname, ".str", sep="")  		#write Structure file name

cat("\nOpening Structure file (", strucfile, ") to read dataset parameters\n", sep="")
dataset = read.table(strucfile, sep = '\t', header =TRUE)
beep(4)

str(dataset)


#locinames
headers <- as.vector(names(dataset))
head(headers)
locilist <- headers[-1]
locilist <- locilist[-1]
head(locilist)
length(locilist)


# Extract the population list per sample from the sample names, the number of loci, etc
snp_num <- ncol(dataset)-2   		#number of loci
samplerows <- length(dataset$X)   		#number of samples
samplenum <- (length(dataset$X)/2)
odd_indexes<-seq(1,samplerows, 2)   			#list of odd indexes, because there is two rows per sample
samps <- as.data.frame(factor(dataset[[1]]))   		#list of samples in the file


#SAVE SAMPLES
newnames = substr(as.vector(samps[odd_indexes,1]), 1, poplength)   		# Extract only odd sample names and take only the first characters of its name
samplelist = as.data.frame(samps[odd_indexes,1]) # Extract a list of samples without duplicates

popnames <- unique(newnames)   		#list of populations
popnames
tempop <- newnames

####### IF YOU WANT TO CHANGE THE POPULATION NAMES THIS IS THE MOMENT
	
# group F1s from the same cross-type together
fixedpops <- gsub('^W([A-Z]{2})$', 'F\\1', tempop)
fixedpops <- gsub('^X([A-Z]{2})$', 'F\\1', fixedpops)

# ONLY IF ALL WILD POPS: Delete "N"
fixedpops <- gsub('^N([A-Z]{2})$', '\\1', fixedpops)




newnames <- fixedpops
popnames <- unique(newnames)   		#list of populations
popnames




	#######################################
	#####    IMPORT FILE TO ADEGENET  #####
	#######################################




rawgenind <- read.structure(strucfile, row.marknames=1, onerowperind=FALSE, n.ind=samplenum, n.loc=snp_num, col.lab=1, col.pop=2, NA.char = "0", ask=FALSE)
beep(3)
cat("\nDone.\n", sep=" ")
#Compare both imports
rawgenind

#replace nums by popnames
pop(rawgenind)
pop(rawgenind) <- newnames
pop(rawgenind)
maxpop <- length(levels(pop(rawgenind)))


#create groups per sp and year
newnames
info <- newnames
popyear <- gsub('^Z([A-Z]{2}.*)$', 'Psic 17 \\1', info)
popyear <- gsub('^Y([A-Z]{2}.*)$', 'Psic 18 \\1', popyear)
popyear <- gsub('^F([A-Z]{2}.*)$', 'Zoo F1s \\1', popyear)
popyear <- gsub('^N(VC)$', 'Psic 19 \\1', popyear)
popyear <- gsub('^N(DU)$', 'Psic 19 \\1', popyear)
popyear <- gsub('^N(OS)$', 'Psic 19 \\1', popyear)
popyear <- gsub('^N(KL)$', 'Psic 19 \\1', popyear)
popyear <- gsub('^N(RK)$', 'Psic 19 \\1', popyear)
popyear <- gsub('^N(OB)$', 'Psic 19 \\1', popyear)
popyear <- gsub('^N(ST)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(PJ)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(SC)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(BJ)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(KP)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(PK)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(PM)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(PG)$', 'Psic 16 \\1', popyear)
popyear <- gsub('^N(.*)$', 'Pmel 16 \\1', popyear)
popyear


#add info
extrainfo <- as.data.frame(popyear)
# either all
extrainfo["Pops"] <- newnames
strata(rawgenind) <- extrainfo
head(strata(rawgenind))








		########################################
		####   HARDY-WEINBERG EQUILIBRIUM   ####
		########################################


##### PER POPULATIONS

#fixnames loci
head(locilist)
locilist <- gsub('^X([0-9]*).([0-9]*)..$', '\\1:\\2', locilist)

#build a dataframe with locinames on the rows
locinames <- as.data.frame(locilist)
basetable <- as.data.frame(locilist)
head(basetable)

#split datafile per pop
genpopdata <- seppop(rawgenind, ~Pops)
beep(4)



for (k in genpopdata) { 

	popname <- print(levels(k@pop))
	cat("\n\n", popname, "\n", sep="")
	pophw <- hw.test(k, B=iterations)
	beep(1)
	popHWtable <- as.data.frame(pophw)
	#head(popHWtable)
	#popHWtable$Pr.exact
	signifpval <- popHWtable[which (popHWtable$Pr.exact < refval),]
	#signifpval
	
	#length(rownames(signifpval))
	locinum <- length(locinames$locilist)
	signum <- length(rownames(signifpval))
	proportionsig = signum / locinum
	percentsig = round(proportionsig*100, digits=2)
	cat("\n", popname, ": ", signum, " out of ", locinum, " were significantly out of HWe (", percentsig, "% p-values were below ", refval, ")\n\n", sep="")
	


	#subset dataframe
	#str(popHWtable)
	head(popHWtable)
	infoHWtable <- as.data.frame(popHWtable$Pr.exact)
	#head(infoHWtable)
	
	
	#fix locinames (may be not necessary)
	longnames <- rownames(popHWtable)
	#head(longnames)
	locishort <- gsub('^([0-9]*):([0-9]*):.*$', '\\1:\\2', longnames)
	head(locishort)
	infoHWtable$locilist <- locishort
	
	#add population name
	#head(infoHWtable)
	newhead <- c(popname, "locilist")
	names(infoHWtable) <- newhead
	
	#merge
	merged_table <- merge(basetable, infoHWtable, by.x="locilist", by.y="locilist", sort=TRUE)
	head(merged_table)
	
	
	
	basetable <- merged_table
	cat("\n\n", sep="")
}


head(basetable)
beep(4)


#subset dataframe with only p-values
lastpop = maxpop + 1
lociinfo <- basetable[,c(2:lastpop)]
head(lociinfo)

#add a column to basetable to hold the number of pops in which the loci is out of hwe
locinum <- length(locilist)
outhw <- numeric(locinum)
basetable$outhw <- outhw
lastcol = lastpop + 1
locirange = c(1:locinum)


#loop and add the number of populations in which each loci is out of hwe

deletelist = c()


for (i in locirange) {
	#number of pops in whic locus is out of hw
	outnum <- sum(lociinfo[i,] < refval)
	#proportion with total num of pops
	hwi_rate <- (outnum/maxpop)
	#add to the table
	basetable[i, lastcol] <- hwi_rate
	#if proportion is too high, add locus name to delete list
	locusname <- as.character(basetable[i,1])
	#cat("\n", locusname, "\n\n", sep="")
	if (hwi_rate >= ref_outhw) { deletelist <- c(deletelist, locusname) }
}

deletelength <- length(deletelist)
finalnum = locinum - deletelength
#deletelength
cat("\n\nThe \"worst\" locus is out of HWe in ", max(basetable$outhw)*maxpop, " out of ", maxpop, " populations.\nThat is a frequency of ", max(basetable$outhw), "; only loci above ", ref_outhw, " were removed.\nThat is ", deletelength, " loci removed from a total of ", locinum, ".\n", sep="")
deletelist
cat("There will be ", finalnum, " loci left\n\n", sep="")
beep(6)




#transform to genlight

#fix locinames
#names(rawgenind@all.names) <- locinames$locilist
#head(names(rawgenind@all.names))
rawgenind
genlightdata <- gi2gl(rawgenind, parallel = TRUE)
genlightdata
#add loci names
genlightdata@loc.names <- locinames$locilist
genlightdata@other$loc.metrics <- as.data.frame(locinames)


#delete loci
genlightequil <- gl.drop.loc(genlightdata, deletelist, v=2)

length(genlightequil@loc.names)
genlightequil

#save
outname = paste(inputname, "_HWe.str", sep="")

#genlightequil@ind.names

gl2structure(genlightequil, indNames = genlightequil@ind.names, addcolumns = genlightequil@pop, ploidy = 2, exportMarkerNames = TRUE, outfile = outname, outpath = ".", v = 1)
cat("\nDONE!!\n\n")
beep(5)

























##################################################################################################################
##################################################################################################################
##################################################################################################################
########################		NOT NEEDED SO FAR, DID NOT FINISH YET		######################################
##################################################################################################################
##################################################################################################################
##################################################################################################################




##############        GLOBAL HWE        ###############


#Do the test
allHWpval <- hw.test(rawgenind, B=iterations)
beep(4)





#check output
head(allHWpval, n=10L)
str(allHWpval)

#weird structure? How to know how many are significant?
# Usually we will check the p-values of HW like this:
allHWpval$Pr.exact
#unsurprisingly the HWeq test of a gengind object is not a dataframe

is.data.frame(allHWpval)
is.atomic(allHWpval)
# Is a vector, we can not access variables with "$"

#transform
allHWtable <- as.data.frame(allHWpval)
head(allHWtable)
str(allHWtable)
nrow(allHWtable)
is.data.frame(allHWtable)

#Now extract the significant p-values from the dataframe
#we use the which to extract tto make a subset, but now instead of the rows which are "male"/"female" or "PK"/"PM"
#we extract the rows wich have a value below 0.05 to a new subset
sigNOHW <- allHWtable[which (allHWtable$Pr.exact <refval),]
#check the number of rows of the new dataframe to know how many loci
nrow(sigNOHW)

noHWloci <- row.names(sigNOHW)

#   2.1 - Check the dataset we opened at the begining of the class, and the loci to delete
str(dataset, list.len=10)
head(noHWloci)


#fixnames loci
lociweird <- names(dataset)
head(lociweird)
locifixed <- gsub('^X([0-9]*).([0-9]*)..$', '\\1:\\2', lociweird)
head(locifixed)
names(dataset) <- locifixed
head(names(dataset))

shitloci <- as.vector(gsub('^([0-9]*:[0-9]*):.*$', '\\1', noHWloci))
head(shitloci)


#Make a subset of dataset with the columns that do not (!) match (%in%) the names stored
filtered <- dataset[ , !(names(dataset) %in% shitloci)]

dim(filtered)
head(filtered)


#Save as a file

savename = paste(inputname, "_WHe.str", sep="")
write.table(filtered,file=savename,sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

# delete the first two headers from the file (X and X.1)
#


