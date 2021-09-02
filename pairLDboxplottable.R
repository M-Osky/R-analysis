#!/usr/bin/Rscript

#script to output a table in which each LD r2 value is categorized in size groups

#to read options 
library(optparse)

#######################   		Parsing command line arguments
cat("\n\n#################################################\n\nStarting: Reading command line arguments...\n\n")


#> set interface
option_list <- list(

  make_option(c("-f", "--flag"),    action="store_true", default=FALSE,   help="Flag [default: %default]"),
  make_option(c("-r", "--ld"),   type="character",    default=FALSE,   help="Plink LD analysis output.ld; analysis done for all pairs per scaffold. Required"),
  make_option(c("-b", "--bin"), type="integer",      default=30,     help="Number of distance bins in which to clasify the R-squared values [default: %default]")
)

parser <- OptionParser(usage="%prog [options]\n\nUse this script to print summary tables from Plink LD analysis.\n\nIf you get an \'X11 display\' error it means you need to set a virtual framebuffer X11 server.\nCheck \"Initiate and declare a virtual graphical environment\" in the beginning of the script.\n", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)

opt <- args$options
file <- args$args



#######################   		LIBRARIES

cat("\n\n>Setting virtual display environment and loading packages.\n")

###### Initiate and declare a virtual graphical environment #########
# The line below is very important, if there is not display environment for X11 you'll need a virtual framebuffer X11 server: Xvfb
	#Uncomment the line below the first time you run it, and don't worry if you get a "Fatal server error: Server is already active for display 0" 
	#Unix is very melodramatic but that just means the command was not needed. Line below for virtual framebuffer X11 server:

#system("Xvfb :0 -ac -screen 0 1960x2000x24 &")   		#<-- This line
Sys.setenv("DISPLAY"=":0")




############################################################################
####      Script, you know... don't touch anything unless... bla bla bla
###########################################################################

#< set interface
options(scipen = 999)




inputfile <- opt$ld
nbin <- opt$bin








#open and check
inputname <-gsub('^(.*).ld$', '\\1', inputfile)
cat("\n\nReading ", inputname, "\n\n", sep="")
dataset = read.table(inputfile, header=TRUE)


#calculate pairwise bp distances
dataset$DIST = abs(dataset$BP_B - dataset$BP_A)

cat(capture.output(head(dataset)), sep="\n")






#############################################


r = pretty(max(dataset$DIST), n=1)[2]
r

#select only distance and r2
str(dataset)
subdata = dataset[,8:7]
str(subdata)
names(subdata) <-c("binsize", "R2")

#length(subdata$DIST)
#str(subdata)
cat("\n\nNow clasifying all distances in ", nbin, " bins\n\n", sep="")
#define size of the bins
binfreq=r/nbin
cat("\nEach bin every ", binfreq, "bp\n", sep="")
myseq = seq(0, r, by=binfreq)   		#save the seq of limits between bins
#myseq
numbins = length(myseq)-1
#numbins

myrg = c(1:numbins)
#myrg


#b=2  ####debug
maxdata=nrow(subdata)
wholedata=c(1:maxdata)

eachval =1  #debug
for (eachval in wholedata) {
	cat(eachval, " of ", maxdata, "\n", sep="")
	mydist = subdata[eachval, 1]
	mydist
	b=1   #debug
	#calculate the middle value
	for (b in myrg) {
		minor = myseq[b]
		minor
		major = myseq[b+1]
		major
		midval=mean(c(minor, major))
		midval
		label=round(midval/1000)
		label
		
		if (mydist > minor && mydist < major) {
			subdata[eachval, 3] <- label
			
			liminf = round((minor/1000),0)
			limsup = round((major/1000),0)
			interv=paste(liminf,"-",limsup, "kb", sep="")
			interv
			subdata[eachval, 1] <- interv
		}
	}
}

boxeddata <- subdata[,c(1, 3, 2)]
names(boxeddata)<-c("interval", "midvalue", "R2" )

head(boxeddata)
str(boxeddata)



forboxplots=paste(inputname, "allsize-r2.txt", sep="")
write.table(boxeddata, file=forboxplots, quote=FALSE, sep="\t", row.names=FALSE)

cat("\n\n\nAll done!!\n\n\n")
