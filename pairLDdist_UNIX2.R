#!/usr/bin/Rscript

#script to calculate average LD and pb distances
#input file is a R-squared output from plink --r2
#it should end in .ld (or change the extennsion filtering)
#column headers matter, so don't change them
#it is done for analysis calculated per scaffold/chrom so if the flag --inter-chr was used results may not be accurate

#to read options 
library(optparse)

#######################   		Parsing command line arguments
cat("\n\n#################################################\n\nStarting: Reading command line arguments...\n\n")


#> set interface
option_list <- list(

  make_option(c("-f", "--flag"),    action="store_true", default=FALSE,   help="Flag [default: %default]"),
  make_option(c("-r", "--ld"),   type="character",    default=FALSE,   help="Plink LD analysis output.ld; analysis done for all pairs per scaffold. Required"),
  make_option(c("-b", "--bim"), type="character",      default=FALSE,     help="Plink input file.bim with all the loci information. Optional")
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




inputfile <- opt$ld
inputbim <- opt$bim








#open and check
inputname <-gsub('^(.*).ld$', '\\1', inputfile)

cat("\n\nReading ", inputname, "\n\n", sep="")

dataset = read.table(inputfile, header=TRUE)
cat(capture.output(head(dataset)), "\n", sep="\n")



################################################################
########				RUN!							########
################################################################
options(scipen = 999)

scaffoldlist = unique(dataset$CHR_A)
numscaf = length(scaffoldlist)
cat("\n\n", numscaf, " scaffolds detected\n\n", sep="")
cat("\n", scaffoldlist, "\n", sep=" ")
rangescaf = c(1:numscaf)

#calculate pairwise bp distances
dataset$DIST = abs(dataset$BP_B - dataset$BP_A)
#head(dataset)

#create an empty dataframe with same format
bigscafdf = dataset[which(is.na(dataset$CHR_A)), ]
#bigscafdf
#k=1 #debug


#create dataframes for the outputs
#general scaffold data information
graltable = data.frame(scaffold=character(numscaf), min_length=integer(numscaf), numSNPs=integer(numscaf), min_dis=integer(numscaf), max_dis=integer(numscaf), mean_dis=double(numscaf), SD_dis=double(numscaf), min_r2=double(numscaf), max_r2=double(numscaf), mean_r2=double(numscaf), SD_r2=double(numscaf), r2bel02=double(numscaf), r2mid=double(numscaf), r2great08=double(numscaf), num_paircomp=integer(numscaf))
#str(graltable)


#number of size bins per scaffold (in 100000 bp increments)
nbins=10
rbins=4 #number of r2 classes
sizebins=c(1:nbins)
#headers for the distance-r2 table
numrowdist = numscaf*rbins*nbins
#numrowdist = numscaf*nbins
#binsdf = data.frame(scaffold=character(numrowdist), dist_bin=integer(numrowdist), num_pairs=integer(numrowdist), mean_r2=double(numrowdist), sd_r2=double(numrowdist), r2_00_02=double(numrowdist), r2_02_04=double(numrowdist), r2_04_06=double(numrowdist), r2_06_10=double(numrowdist))
binsdf = data.frame(scaffold=character(numrowdist), sizebin_Mb=character(numrowdist), r2_range=character(numrowdist), r2_propor=double(numrowdist), num_pairs=integer(numrowdist))
#str(binsdf)


i=1
for (k in rangescaf) {
	
	#k=2 #debug
	
	#get scaffold name
	scafname = scaffoldlist[k]
	cat("\n", scafname, "\n", sep="")
	graltable[k,1] <- scafname
	#graltable
	
	#select only rows from that scaffold (first column)
	scafdata <- dataset[which (dataset$CHR_A == scafname),]
	
	#save number of comparisons
	numpairwcomp=length(scafdata$SNP_B)
	graltable[k,15] <- numpairwcomp

	#calculate length of scaffold
	alldist=c(scafdata$BP_A, scafdata$BP_B)
	minlength = max(alldist)
	#minlength
	graltable[k,2] <- minlength
	#graltable
	
	
	#number of SNPs
	allnames = c(scafdata$SNP_A, scafdata$SNP_B)
	numsnps = length(unique(allnames))
	graltable[k,3] <- numsnps
	
	#dces
	mindis = min(scafdata$DIST)
	graltable[k,4] <- mindis
	maxdis = max(scafdata$DIST)
	graltable[k,5] <- maxdis
	meandis = mean(scafdata$DIST)
	graltable[k,6] <- meandis
	sddis = sd(scafdata$DIST)
	graltable[k,7] <- sddis
	
	#r-squared
	minr2 = min(scafdata$R2)
	graltable[k,8] <- minr2
	maxr2 = max(scafdata$R2)
	graltable[k,9] <- maxr2
	meanr2 = mean(scafdata$R2)
	graltable[k,10] <- meanr2
	sdr2 = sd(scafdata$R2)
	graltable[k,11] <- sdr2
	
	#r2 percent
	totalpairs=length(scafdata$R2)
	#totalpairs
	lowr2 <- scafdata[which (scafdata$R2 < 0.2),]
	lownum = length(lowr2$R2)
	#lownum
	#percent of R2 < 0.2
	lowperc = (lownum / totalpairs)*100
	graltable[k,12] <- lowperc
	
	inter2 <- scafdata[which (scafdata$R2 >= 0.2),]
	midr2 <- scafdata[which (inter2$R2 <= 0.8),]
	midnum = length(midr2$R2)
	#midnum
	#percent of R2 E [0.2, 0.8]
	midperc = (midnum / totalpairs)*100
	graltable[k,13] <- midperc
	
	higr2 <- scafdata[which (scafdata$R2 > 0.8),]
	hignum = length(higr2$R2)
	#hignum
	#percent of R2 > 0.8
	higperc = (hignum / totalpairs)*100
	graltable[k,14] <- higperc
	#head(graltable)
	
	
	
	
	
	#DECAY CALCULATIONS PER CHROM
	
	#create the headers for the table
	#str(binsdf)
	
	#head(binsdf)
	#r2 per dist
	
	if (length(scafdata[which (scafdata$DIST > 1000000),])>=4) {
	
	bigscafdf <- rbind(bigscafdf, scafdata)   #save all valid scaffold data for later
	
		for (n in sizebins) {
			#add scaffold name to df
			
			#n=2   #debug
			
			#calculate range of distances for this loop
			distmax = n*100000
			#distmax
			below=n-1
			distmin = below*100000
			#distmin
			if (n < 10 ) { above = paste("0.", n, sep="") } else { above =format(n/10, nsmall = 1)}
			binname <- paste("0.", below, "-", above, sep="") 
			#binname
			
			
			
			
			#select pairs with that dist and save some descriptive
			dist00_max <- scafdata[which (scafdata$DIST < distmax),]
			distmin_max <- dist00_max[which (dist00_max$DIST >= distmin),]
			pairnumdist = length(distmin_max$DIST)
			
			#pairnumdist
			#avgr2dist = mean(distmin_max$R2)
			#sdr2dist = sd(distmin_max$R2)
			#binsdf[i,3] <- pairnumdist
			#binsdf[i,4] <- avgr2dist
			#binsdf[i,5] <- sdr2dist
			#head(binsdf)
			
			#select pairs with a r2 range
			#R2 00-0.2
			r2_00_02 <- distmin_max[which (distmin_max$R2 < 0.2),]
			#i=1
			binsdf[i,1] <- scafname
			binsdf[i,2] <- binname
			binsdf[i,3] <- "0.0-0.2"
			numr2_00_02=length(r2_00_02$R2)
			propr2_00_02 = (numr2_00_02 / pairnumdist)
			binsdf[i,4] <- propr2_00_02
			binsdf[i,5] <- numr2_00_02
			i=i+1
			
			#R2 0.2-0.4
			r2_00_04 <- distmin_max[which (distmin_max$R2 < 0.4),]
			r2_02_04 <- r2_00_04[which (r2_00_04$R2 >= 0.2),]
			binsdf[i,1] <- scafname
			binsdf[i,2] <- binname
			binsdf[i,3] <- "0.2-0.4"
			numr2_02_04=length(r2_02_04$R2)
			propr2_02_04 = (numr2_02_04 / pairnumdist)
			binsdf[i,4] <- propr2_02_04
			binsdf[i,5] <- numr2_02_04
			i=i+1
			
			#R2 0.4-0.6
			r2_00_06 <- distmin_max[which (distmin_max$R2 < 0.6),]
			r2_04_06 <- r2_00_06[which (r2_00_06$R2 >= 0.4),]
			binsdf[i,1] <- scafname
			binsdf[i,2] <- binname
			binsdf[i,3] <- "0.4-0.6"
			numr2_04_06=length(r2_04_06$R2)
			propr2_04_06 = (numr2_04_06 / pairnumdist)
			binsdf[i,4] <- propr2_04_06
			binsdf[i,5] <- numr2_04_06
			i=i+1
			
			#R2 0.6-1
			r2_06 <- distmin_max[which (distmin_max$R2 >= 0.6),]
			binsdf[i,1] <- scafname
			binsdf[i,2] <- binname
			binsdf[i,3] <- "0.6-1.0"
			
			numr2_06=length(r2_06$R2)
			propr2_06_10 = (numr2_06 / pairnumdist)
			binsdf[i,4] <- propr2_06_10
			binsdf[i,5] <- numr2_06
			i=i+1
			#do the next set of distance bins
			
		}
	
	}
	#check next scaffold
	
	
}


#save
head(graltable)
gralname=paste(inputname, "_scaffoldsummary.txt", sep="")
write.table(graltable, file=gralname, quote=FALSE, sep="\t", row.names=FALSE)


# #plot size x num SNPs
# graltable$lengthkb <- graltable$min_length/1000

# plot1=paste(inputname, "_size-num.png", sep="")
# png(plot1, units="in", width=10, height=10, res=300)
# plot(x=graltable$lengthkb, y=graltable$numSNPs, main="Number of SNPs per Scaffold length", xlab="min scaffold length (kb)", ylab="Number of SNPs per scaffold", col="darkolivegreen4", pch=19)
# dev.off()


# #plot size x mean R2
# plotsr=paste(inputname, "_size-r2_all.png", sep="")
# png(plotsr, units="in", width=10, height=10, res=300)
# plot(x=graltable$lengthkb, y=graltable$mean_r2, main=expression(paste("Average R"^2," per Scaffold length")), xlab="min scaffold length (kb)", ylab=expression(paste("Mean R"^2, " between loci pairs")), col="darkorange", pch=19)
# dev.off()




# #plot number of snps per scaffold
# graltable$scafnum <-gsub('^[A-Z]*[a-z]*([0-9]*)$', '\\1', graltable$scaffold)
# ymaxlim = pretty(max(graltable$numSNPs), n=1)[2]

# plot2=paste(inputname, "_scafnum-snps.png", sep="")
# png(filename=plot2, units="in", width=25, height=7, res=900)
# print(
# ggplot(data=graltable, aes(x=reorder(scafnum, as.numeric(scafnum)), y=numSNPs)) + geom_bar(stat="identity", fill="darkolivegreen4") + ggtitle ("Number of SNPs per scaffold") + 
 # theme(axis.text.x = element_text(angle = 90, size=5, vjust=0.2, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) +
 # scale_y_continuous(expand = c(0,20), breaks=seq(0, ymaxlim, 100)) + scale_x_discrete(expand = c(0.005,0)) + xlab("Scaffold number") + ylab ("Number of SNPs")
 
# )
# dev.off()


# #plot length per scaffold
# ymaxlim = pretty(max(graltable$lengthkb), n=1)[2]

# plot3=paste(inputname, "_scafnum-length.png", sep="")
# png(filename=plot3, units="in", width=25, height=7, res=900)
# print(
# ggplot(data=graltable, aes(x=reorder(scafnum, as.numeric(scafnum)), y=lengthkb)) + geom_bar(stat="identity", fill="goldenrod3") + ggtitle ("Scaffold length (kb)") + xlab("Scaffold number") + ylab ("") +
 # theme(axis.text.x = element_text(angle = 90, size=5, vjust=0.2, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) +
 # scale_y_continuous(expand = c(0,1000), breaks=seq(0, ymaxlim, 10000)) + scale_x_discrete(expand = c(0.005,0))
 
# )
# dev.off()




#save
head(binsdf)

scafbinsname=paste(inputname, "_perscaffold-bins.txt", sep="")
write.table(binsdf, file=scafbinsname, quote=FALSE, sep="\t", row.names=FALSE)

#plot bins

# #filter scaffolds
# head(binsdf)
# badscaffold <- unique(binsdf[which (binsdf$num_pairs < 5),]$scaffold)
# allscafbins = unique(binsdf$scaffold)

# scaflist = setdiff(allscafbins, badscaffold)

# for (scafname in scaflist) {
	# plotfile = paste(inputname, "_", scafname, "r2bins.png", sep="")
	# binsonesc <- binsdf[which (binsdf$scaffold == scafname),]
	# mypaletteheat=c("tomato3", "sienna2", "tan1", "khaki")
	
	# png(filename=plotfile, units="in", width=10, height=7, res=300)
	# print(
	# ggplot(data=binsonesc, aes(x=sizebin_Mb, y=r2_propor)) + geom_col(aes(fill=r2_range), width=0.95, position = position_fill(reverse = TRUE)) +
	 # scale_fill_manual(expression(paste("Ranges of r"^2)), values=mypaletteheat)+ scale_colour_manual(values=mypaletteheat) + ggtitle(scafname) + xlab("Distance bins (Mb)") + ylab (expression(paste("proportion of r"^2," classes"))) +
	 # theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) + 
	 # scale_y_continuous(expand = c(0,0))
	
	# )
	# dev.off()

# }

#now, all
cat("\n\nAll done, now calculating proportion of SNPs per R2 bin overall Scaffolds\n\n")


#bigscafdf

binstot = binsdf[which(is.na(binsdf$scaffold)), ]

i=1
scafname="All scaffold with distances above above 1Mb"
for (n in sizebins) {
	
	#calculate range of distances for this loop
	distmax = n*100000
	#distmax
	below=n-1
	distmin = below*100000
	#distmin
	if (n < 10 ) { above = paste("0.", n, sep="") } else { above =format(n/10, nsmall = 1)}
	binname <- paste("0.", below, "-", above, sep="") 
	#binname
	
	
	
	#select pairs with that dist and save some descriptive
	dist00_max <- bigscafdf[which (bigscafdf$DIST < distmax),]
	distmin_max <- dist00_max[which (dist00_max$DIST >= distmin),]
	pairnumdist = length(distmin_max$DIST)
	
	
	
	#select pairs with a r2 range
	#R2 00-0.2
	r2_00_02 <- distmin_max[which (distmin_max$R2 < 0.2),]
	#i=1
	binstot[i,1] <- scafname
	binstot[i,2] <- binname
	binstot[i,3] <- "0.0-0.2"
	numr2_00_02=length(r2_00_02$R2)
	propr2_00_02 = (numr2_00_02 / pairnumdist)
	binstot[i,4] <- propr2_00_02
	binstot[i,5] <- numr2_00_02
	i=i+1
	
	#R2 0.2-0.4
	r2_00_04 <- distmin_max[which (distmin_max$R2 < 0.4),]
	r2_02_04 <- r2_00_04[which (r2_00_04$R2 >= 0.2),]
	binstot[i,1] <- scafname
	binstot[i,2] <- binname
	binstot[i,3] <- "0.2-0.4"
	numr2_02_04=length(r2_02_04$R2)
	propr2_02_04 = (numr2_02_04 / pairnumdist)
	binstot[i,4] <- propr2_02_04
	binstot[i,5] <- numr2_02_04
	i=i+1
	
	#R2 0.4-0.6
	r2_00_06 <- distmin_max[which (distmin_max$R2 < 0.6),]
	r2_04_06 <- r2_00_06[which (r2_00_06$R2 >= 0.4),]
	binstot[i,1] <- scafname
	binstot[i,2] <- binname
	binstot[i,3] <- "0.4-0.6"
	numr2_04_06=length(r2_04_06$R2)
	propr2_04_06 = (numr2_04_06 / pairnumdist)
	binstot[i,4] <- propr2_04_06
	binstot[i,5] <- numr2_04_06
	i=i+1
	
	#R2 0.6-1
	r2_06 <- distmin_max[which (distmin_max$R2 >= 0.6),]
	binstot[i,1] <- scafname
	binstot[i,2] <- binname
	binstot[i,3] <- "0.6-1.0"
	
	numr2_06=length(r2_06$R2)
	propr2_06_10 = (numr2_06 / pairnumdist)
	binstot[i,4] <- propr2_06_10
	binstot[i,5] <- numr2_06
	i=i+1
	#do the next set of distance bins
	
}



wholebinsname=paste(inputname, "all-bins.txt", sep="")
write.table(binstot, file=wholebinsname, quote=FALSE, sep="\t", row.names=FALSE)
#binstot


# ###plot
# plotfile = paste(inputname, "_All_bigScaff_r2bins.png", sep="")
# mypaletteheat=c("tomato3", "sienna2", "tan1", "khaki")

# png(filename=plotfile, units="in", width=10, height=7, res=300)
# print(
# ggplot(data=binstot, aes(x=sizebin_Mb, y=r2_propor)) + geom_col(aes(fill=r2_range), width=0.95, position = position_fill(reverse = TRUE)) +
 # scale_fill_manual(expression(paste("Ranges of r"^2)), values=mypaletteheat)+ scale_colour_manual(values=mypaletteheat) + ggtitle(scafname) + xlab("Distance bins (Mb)") + ylab (expression(paste("proportion of r"^2," classes"))) +
 # theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) + 
 # scale_y_continuous(expand = c(0,0))
# )
# dev.off()






cat("\n\nAll done, now calculating average R2 (no num per bin) overall Scaffolds\n\n")

sizegroups=c(100000, 1000000, 10000000, 90000000)



#r=90000000   #debug
#r=10000000   #debug
#r=1000000   #debug
#r=100000   #debug

for (r in sizegroups) {
	cat("calculating max dist: ", r, "\n", sep="")
	#select only pairs with distance below certain number
	subdata <- dataset[which (dataset$DIST <= r),]
	#length(subdata$DIST)
	#str(subdata)
	
	r2decayname=paste(inputname, "avgR2dist", r, "overall.txt", sep="")
	cat("avgr2", "bins\n", file=r2decayname, sep="\t")
	
	
	#define size of the bins
	binfreq=r/10
	#binfreq
	myseq = seq(0, r, by=binfreq)   		#save the seq of limits between bins
	numbins = length(myseq)-1
	#numbins
	
	#create a dataframe for the mean R2 results per bin
	tablebins=data.frame(meanR2=double(numbins))
	#tablebins
	myrg = c(1:numbins)
	
	#b=1   #debug
	#myseq
	#myrg
	#calculate the middle value
	for (b in myrg) {
		minor = myseq[b]
		major = myseq[b+1]
		label=round(mean(c(minor, major))/1000)
		if (b == 1) {
			rowlabels = label
		} else {
			rowlabels=c(rowlabels, label)
		}
		
	}
	
	#add middle values from each bin as rownames
	rownames(tablebins) <- myrg
	tablebins$bins <- rowlabels
	names(tablebins) <- c("avgr2", "bins")
	#tablebins
	
	#generate a table with the mean r2 per bin
	for (g in myrg) {
		#g=2 ##debug##
		minsize = myseq[g]
		#minsize
		maxsize = myseq[g+1]
		#maxsize
		binmax <- subdata[which (subdata$DIST <= maxsize),]
		length(binmax$R2)
		bindist <- binmax[which (binmax$DIST >= minsize),]
		length(bindist$R2)
		binmean = mean(bindist$R2)
		#binmean
		tablebins[g,1] <- binmean
		
	}
	
	
	
	
	tablebins
	
	write.table(tablebins, file=r2decayname, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
	
	# #summarize the distances size used for the plot
	# if(r>1000000) { 
		# sumsize = r/1000000
		# titsize = paste("0-", sumsize, " Mb", sep="")
	# } else if (r >= 1000) {
		# sumsize = r/1000
		# titsize = paste("0-", sumsize, " kb", sep="")
	# } else { 
		# titsize = paste("0-", r, sep="")
	# }
	
	
	# #plot
	# mytitle=paste("Overall LD decay for intermarquer distances: ", titsize, sep="")
	
	# #y limits
	# ymax = max(tablebins$avgr2, na.rm = TRUE)+0.1
	# if (ymax > 1) { ymax = 1 }
	# yrange=c(0, ymax)
	
	# plotname=paste(inputname, "_LDdecay_", titsize, ".png", sep="")
	
	# png(plotname, units="in", width=12, height=7, res=300)
	# #plot(x=tablebins$bins, y=tablebins$avgr2, main=mytitle, xlab="Maker distance (kb)", ylab="Average r2", ylim=yrange, xaxt="n", col="darkolivegreen4", pch=19)
	# #plot(x=tablebins$bins, type="b", y=tablebins$avgr2, main=mytitle, xlab="Maker distance (kb)", ylab=expression("Average r"^2), ylim=yrange, xaxt="n", col="darkolivegreen4", pch=19)
	# plot(x=tablebins$bins, type="b", y=tablebins$avgr2, main=mytitle, xlab="Maker distance (kb)", ylab="", ylim=yrange, xaxt="n", col="darkolivegreen4", pch=19)
	# title(ylab=expression("Average r"^2), line=2.2, cex.lab=1)
	# axis(side=1, at=rowlabels)
	# dev.off()
	
 }


#repeat with median R2

for (r in sizegroups) {
	cat("max dist: ", r, "\n", sep="")
	r2decayname=paste(inputname, "medianR2dist", r, "overall.txt", sep="")
	cat("medianr2", "bins\n", file=r2decayname, sep="\t")
	#select only pairs with distance below certain number
	subdata <- dataset[which (dataset$DIST <= r),]
	#length(subdata$DIST)
	#str(subdata)
	
	#define size of the bins
	binfreq=r/10
	#binfreq
	myseq = seq(0, r, by=binfreq)   		#save the seq of limits between bins
	numbins = length(myseq)-1
	#numbins
	
	#create a dataframe for the mean R2 results per bin
	tablebins=data.frame(meanR2=double(numbins))
	#tablebins
	myrg = c(1:numbins)
	
	#b=1   #debug
	#myseq
	#myrg
	#calculate the middle value
	for (b in myrg) {
		minor = myseq[b]
		major = myseq[b+1]
		label=round(mean(c(minor, major))/1000)
		if (b == 1) {
			rowlabels = label
		} else {
			rowlabels=c(rowlabels, label)
		}
		
	}
	
	#add middle values from each bin as rownames
	rownames(tablebins) <- myrg
	tablebins$bins <- rowlabels
	names(tablebins) <- c("medianval", "bins")
	#tablebins
	
	#generate a table with the mean r2 per bin
	for (g in myrg) {
		#g=10 ##debug##
		minsize = myseq[g]
		#minsize
		maxsize = myseq[g+1]
		#maxsize
		binmax <- subdata[which (subdata$DIST <= maxsize),]
		length(binmax$R2)
		bindist <- binmax[which (binmax$DIST >= minsize),]
		length(bindist$R2)
		binmedian = median(bindist$R2)
		#binmedian
		tablebins[g,1] <- binmedian
		
	}
	
	
	#tablebins
	write.table(tablebins, file=r2decayname, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
	
	# #summarize the distances size used for the plot
	# if(r>1000000) { 
		# sumsize = r/1000000
		# titsize = paste("0-", sumsize, " Mb", sep="")
	# } else if (r >= 1000) {
		# sumsize = r/1000
		# titsize = paste("0-", sumsize, " kb", sep="")
	# } else { 
		# titsize = paste("0-", r, sep="")
	# }
	
	
	# #plot
	# mytitle=paste("Overall LD decay for intermarquer distances: ", titsize, sep="")
	
	# #y limits
	# ymax = max(tablebins$medianval, na.rm = TRUE)+0.1
	# if (ymax > 1) { ymax = 1 }
	# yrange=c(0, ymax)
	
	# plotname2=paste(inputname, "_LDdecay_", titsize, "Median.png", sep="")
	
	# png(plotname2, units="in", width=12, height=7, res=300)
	# plot(x=tablebins$bins, type="b", y=tablebins$medianval, main=mytitle, xlab="Maker distance (kb)", ylab="", ylim=yrange, xaxt="n", col="darkolivegreen4", pch=19)
	# title(ylab=expression("Median r"^2), line=2.2, cex.lab=1)
	# axis(side=1, at=rowlabels)
	# dev.off()
	
}










############################################################
############################################################
############################################################





#now the same but for all boxplots

r = pretty(max(dataset$DIST), n=1)[2]
r
#select only pairs with distance below certain number
str(dataset)
subdata = dataset[,8:7]
str(subdata)
names(subdata) <-c("binsize", "R2")

#length(subdata$DIST)
#str(subdata)

#define size of the bins
binfreq=r/20
binfreq
myseq = seq(0, r, by=binfreq)   		#save the seq of limits between bins
myseq
numbins = length(myseq)-1
numbins

myrg = c(1:numbins)
myrg


#b=2  ####debug
maxdata = nrow(subdata)
wholedata=c(1:maxdata)
cat("\n\nClassifying all comparisons in sizebins\n\n", sep="")

#eachval =1  #debug
for (eachval in wholedata) {
	cat(eachval " of ", maxdata, "\n", sep="")
	mydist = subdata[eachval, 1]
	mydist
	#b=1   #debug
	#calculate the middle value
	for (b in myrg) {
		minor = myseq[b]
		minor
		major = myseq[b+1]
		major
		label=round(mean(c(minor, major))/1000)
		label
		
		if (mydist > minor && mydist < major) {
			subdata[eachval, 1] <- label
		}
	}
}

head(subdata)
str(subdata)

##### needs work

forboxplots=paste(inputname, "allsize-r2.txt", sep="")
write.table(subdata, file=forboxplots, quote=FALSE, sep="\t", row.names=FALSE)




# #summarize the distances size used for the plot
# if(r>1000000) { 
	# sumsize = r/1000000
	# titsize = paste("0-", sumsize, " Mb", sep="")
# } else if (r >= 1000) {
	# sumsize = r/1000
	# titsize = paste("0-", sumsize, " kb", sep="")
# } else { 
	# titsize = paste("0-", r, sep="")
# }


# #plot
# mytitle=paste("Overall LD decay for intermarquer distances: ", titsize, sep="")

# #y limits
# ymax = max(tablebins$avgr2, na.rm = TRUE)+0.1
# if (ymax > 1) { ymax = 1 }
# yrange=c(0, ymax)

# plotname2=paste(inputname, "_LDdecay_", titsize, "Median.png", sep="")

# png(plotname2, units="in", width=12, height=7, res=300)
# #plot(x=tablebins$bins, y=tablebins$avgr2, main=mytitle, xlab="Maker distance (kb)", ylab="Average r2", ylim=yrange, xaxt="n", col="darkolivegreen4", pch=19)
# #plot(x=tablebins$bins, type="b", y=tablebins$avgr2, main=mytitle, xlab="Maker distance (kb)", ylab=expression("Average r"^2), ylim=yrange, xaxt="n", col="darkolivegreen4", pch=19)
# plot(x=tablebins$bins, type="b", y=tablebins$avgr2, main=mytitle, xlab="Maker distance (kb)", ylab="", ylim=yrange, xaxt="n", col="darkolivegreen4", pch=19)
# title(ylab=expression("Median r"^2), line=2.2, cex.lab=1)
# axis(side=1, at=rowlabels)
# dev.off()









################################################################################################
################################################################################################

# now check the bim file with the general loci info
bimfile <-gsub('^(.*).bim$', '\\1', inputbim)
cat("\n\nreading bim file: ", bimfile, "\n\n", sep="")



if (inputbim != FALSE) { 

	#open and check
	bimdata = read.table(inputbim, header=FALSE)
	#add headers
	bimdata$totpos <- "NA"
	names(bimdata) <- c("Scaffold", "ID", "morgans", "POS", "ref", "alt", "totpos")
	str(bimdata)
	head(bimdata)
	prevchrom = bimdata[1,1]

	uniquescaf = length(unique(bimdata$Scaffold))

	#create dataframes for the outputs
	#general scaffold data information
	bimscaftable = data.frame(scaffold=character(uniquescaf), min_length=integer(uniquescaf), numSNPs=integer(uniquescaf))
	#str(graltable)
	
	#now a table to plot snps per each 1Mb per scaffold
	bimbins = data.frame(scaffold=character(), numscaf=integer(), Mbsize=integer(), snps=integer(), sorted=integer())
	inistepb=1000000			#bins of 1 Mb
	snpbind=1				#add snps in each loop
	sortbins=1				#bin order
	Mbsize=1				#bin value
	
	stepb=inistepb
	
	#calculate position + scaffold
	numSNPs = length(bimdata$ID)
	rangelines=c(1:numSNPs)
	flag=0
	refchrom=bimdata[1,1]
	scaffcount=1
	snpcount=0
	#snp=1   #debug
	for (snp in rangelines) {
		#take position and chrom name
		size = bimdata[snp,4]
		#size
		chrom = bimdata[snp,1]
		#chrom
		snpcount=snpcount+1
		
		
		if(flag==0 && refchrom == chrom) {
			#if it's the first chromosome just save the position 
			#save length and chrom name to compare
			refchrom = chrom
			refsize = size
			bimdata[snp, 7] <- refsize
		} else if (flag==0 && refchrom != chrom && snpcount > 1) {
			#when new chromosome is encountered
			#save previous chrom data
			bimscaftable[scaffcount, 1] <- refchrom  		#save last chromosome name
			bimscaftable[scaffcount, 2] <- refsize  		#save last SNP position
			snp_pevious=snpcount-1
			bimscaftable[scaffcount, 3] <- snp_pevious  		#save number of snps reached in the previous chrom
			
			refchrom = chrom
			#position in second is equal to length of first crhomosome + position
			newsize = refsize + size
			bimdata[snp, 7] <- newsize
			flag=1
			scaffcount=scaffcount+1
			snpcount=1
		} else if (flag==1 && refchrom == chrom) {
			#for the rest of chromosomes
			refchrom = chrom
			#keep refference size from previous chromosome + position in new
			newsize = refsize + size
			lastsize=size
			bimdata[snp, 7] <- newsize
		} else if (flag==1 && refchrom != chrom) {
			#when a new chrom is found
			#save previous chrom data
			bimscaftable[scaffcount, 1] <- refchrom  		#save last chromosome name
			bimscaftable[scaffcount, 2] <- lastsize  		#save last SNP position
			snp_pevious=snpcount-1
			bimscaftable[scaffcount, 3] <- snp_pevious  		#save number of snps reached in the previous chrom
			
			
			refchrom = chrom
			#take as new reference size the one from last chromosome
			refsize=newsize
			lastsize=size
			#position in second is equal to length of previous crhomosome + position
			newsize = refsize + size
			bimdata[snp, 7] <- newsize
			scaffcount=scaffcount+1
			snpcount=1
		} 
		
		
		
		
		#now a table to plot snps per each 1Mb per scaffold
		#bimbins = data.frame(scaffold=character(), numscaf=integer(), Mbsize=integer(), snps=integer(), sorted=integer())
		#toadd = data.frame(scaffold="scaffold3", numscaf=3, binnum=2, snps=3)
		#newdf <- rbind(bimbins, toadd)
		#stepb=1000000			#bins of 1 Mb
		#snpbind=1				#add snps in each loop
		#sortbins=1				#bin order
		
		
		snppos = bimdata[snp,4]
		chromosome = bimdata[snp,1]
		
		if (snppos < stepb && chromosome == prevchrom) {
			prevchrom = chromosome
			
			if (sortbins==numSNPs) {
				cat("\nAll done, saving data!\n\n")
				sizemb=stepb/inistepb
				numerscaf <-gsub('^[A-Z]*[a-z]*([0-9]*)$', '\\1', prevchrom)
				toadd = data.frame(scaffold="prevchrom", numscaf=numerscaf, Mbsize=sizemb, snps=snpbind, sorted=sortbins)
				bimbins <- rbind(bimbins, toadd)
				sortbins=sortbins+1
				stepb=inistepb
				snpbind=1
				prevchrom = chromosome
			}
			
			snpbind = snpbind + 1
			
		} else if (snppos > stepb && chromosome == prevchrom) {
			#change bin, but first save data from the previous bin
			sizemb=stepb/inistepb
			numerscaf <-gsub('^[A-Z]*[a-z]*([0-9]*)$', '\\1', prevchrom)
			toadd = data.frame(scaffold="prevchrom", numscaf=numerscaf, Mbsize=sizemb, snps=snpbind, sorted=sortbins)
			bimbins <- rbind(bimbins, toadd)
			sortbins=sortbins+1
			stepb=stepb+inistepb
			snpbind=1
			prevchrom = chromosome
		} else if (chromosome != prevchrom) {
			sizemb=stepb/inistepb
			numerscaf <-gsub('^[A-Z]*[a-z]*([0-9]*)$', '\\1', prevchrom)
			toadd = data.frame(scaffold="prevchrom", numscaf=numerscaf, Mbsize=sizemb, snps=snpbind, sorted=sortbins)
			bimbins <- rbind(bimbins, toadd)
			sortbins=sortbins+1
			stepb=inistepb
			snpbind=1
			prevchrom = chromosome
			
		}
	
	}
	#save last chrom data
	bimscaftable[scaffcount, 1] <- refchrom  		#save last chromosome name
	bimscaftable[scaffcount, 2] <- lastsize  		#save last SNP position
	bimscaftable[scaffcount, 3] <- snpcount  		#save number of snps reached in the previous chrom


	str(bimscaftable)
	head(bimscaftable)
	tail(bimscaftable)


	#save
	bimname=paste(bimfile, "_summarybim.txt", sep="")
	write.table(bimscaftable, file=bimname, quote=FALSE, sep="\t", row.names=FALSE)

	bimbinname=paste(bimfile, "_sortedSNPperBin.txt", sep="")
	write.table(bimbins, file=bimbinname, quote=FALSE, sep="\t", row.names=FALSE)
	


	# #plot size x num SNPs
	# bimscaftable$lengthkb <- bimscaftable$min_length/1000

	# plot1b=paste(bimfile, "_size-num_bim.png", sep="")
	# png(plot1b, units="in", width=10, height=10, res=300)
	# plot(x=graltable$lengthkb, y=graltable$numSNPs, main="Number of SNPs per Scaffold length (raw bim file)", xlab="min scaffold length (kb)", ylab="Number of SNPs per scaffold", col="firebrick4", pch=19)
	# dev.off()


	# #plot with number of snps per scaffold
	# bimscaftable$scafnum <-gsub('^[A-Z]*[a-z]*([0-9]*)$', '\\1', bimscaftable$scaffold)
	# ymaxlim = pretty(max(bimscaftable$numSNPs), n=1)[2]

	# plot2b=paste(bimfile, "_scafnum-snps_bim.png", sep="")
	# png(filename=plot2b, units="in", width=25, height=7, res=900)
	# print(
	# ggplot(data=bimscaftable, aes(x=reorder(scafnum, as.numeric(scafnum)), y=numSNPs)) + geom_bar(stat="identity", fill="firebrick4") + ggtitle ("Number of SNPs per scaffold (bim raw file)") + 
	 # theme(axis.text.x = element_text(angle = 90, size=5, vjust=0.2, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) +
	 # scale_y_continuous(expand = c(0,20), breaks=seq(0, ymaxlim, 100)) + scale_x_discrete(expand = c(0.005,0)) + xlab("Scaffold number") + ylab ("Number of SNPs")
	 
	# )
	# dev.off()


	# #plot with length per scaffold
	# ymaxlim = pretty(max(bimscaftable$lengthkb), n=1)[2]

	# plot3b=paste(bimfile, "_scafnum-length_bim.png", sep="")
	# png(filename=plot3b, units="in", width=25, height=7, res=900)
	# print(
	# ggplot(data=graltable, aes(x=reorder(scafnum, as.numeric(scafnum)), y=lengthkb)) + geom_bar(stat="identity", fill="purple4") + ggtitle ("Length (kb) from all scaffolds in bim file") + xlab("Scaffold number") + ylab ("") +
	 # theme(axis.text.x = element_text(angle = 90, size=5, vjust=0.2, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) +
	 # scale_y_continuous(expand = c(0,1000), breaks=seq(0, ymaxlim, 10000)) + scale_x_discrete(expand = c(0.005,0))
	 
	# )
	# dev.off()
	
	
	
	


	# Plot number of SNPs per sorted bin




















}

cat("\n\n\nAll done!\n\n\n\n")
warnings()
