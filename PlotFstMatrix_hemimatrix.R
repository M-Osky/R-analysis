#R sucks!

#plot pairwiseFst calculated in other packages

#packages
library ("adegenet")   	#main genetic package
library ("adegraphics")   	#main genetic package
library(devtools)		#handle github packages
#install_github("alexploner/Heatplus")
library(RColorBrewer)	#colours
library("Heatplus")		#heatmaps
library("ape") #phylogenetic trees






rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets



getwd() #check
#Change working directory if needed
setwd("")


#input file must be an FSTmatrix from Stammp package
#hemimatrix with NA on the top-right and populations in rows and headers
inputname="FSTmatrix_sicula585x39k_m05R7r6h6DP420_1p0p_0m0p"

info="All P. sicula natural populations (LD and HWe filtered before subsampling)"
#info="whatever"


#open and check
inputfile=paste(inputname, ".txt", sep="")
dataset = read.table(inputfile, header=TRUE, row.names=1)
str(dataset)
head(dataset)



# DO YOU NEED THE POPULATIONS IN A SPECIFIC ORDER?

#declare a vector with the order
	#neworder=c("PK", "PM", "KP")
#import a file with the order (one population per line, as in distruct)
	sortedpops="sortedpops"
	sorted = read.table(sortedpops, header=FALSE)
	str(sorted)
	neworder <- as.vector(sorted$V1)
# either if you are alright with the input file order, set "neworder" to "no"
	neworder="no"
#

neworder




########################
########### RUN! #######
########################


#convert negatives to 0
dataset[dataset<0] <- 0

#covert hemimatrix to symmetric matrix
prematrix <- as.matrix(as.dist(dataset))
prematrix


# personalized order
# do not worry with the warning about length > 1
if(neworder != "no") {
	unsorteddf <- as.data.frame(prematrix)
	columnsdf <- unsorteddf[neworder]
	dataset <- columnsdf[order(match(rownames(columnsdf), neworder)), , drop = FALSE]
	datamatrix = as.matrix(dataset)
} else {
	datamatrix=prematrix
}

#datamatrix






#invert order of columns respect rows
revmatrix=datamatrix[,order(ncol(datamatrix):1)]

#calculate the range of values
minval=min(revmatrix)
maxval=max(revmatrix)
rangevals = c(minval,maxval)
#rangevals


#calculate the breaks for the colours in the plot
mybreaks = niceBreaks(rangevals, 7)
numbreaks = length(mybreaks)

#colours
heatpalette = brewer.pal(n = (numbreaks-1), name = 'YlOrRd')


#plots
plotnamedendro=paste(inputname, "_HeatDendro.png", sep="")
png(filename=plotnamedendro, units="in", width=14, height=11, res=300)
plot(regHeatmap(revmatrix, scale="none", dendrogram=NULL, breaks=mybreaks, col=heatpalette))
dev.off()

#par("mar")
#par(mar=c(1,1,1,1))

plotnameheat=paste(inputname, "_HeatMap.png", sep="")
png(filename=plotnameheat, units="in", width=14, height=10, res=300)
heatmap_2(revmatrix, scale="none", do.dendro=c(FALSE, FALSE), breaks=mybreaks, col=heatpalette, Rowv=NA, Colv=NA, legend=2)
dev.off()


#To appreciate even better the genetic distances between each pair build a tree


#use the command nj() to create a neighbour joining tree from: FstsMatrix
psic.tree <- nj(datamatrix)	   #NETWORK JOINING TREE
#plot

mypalette=rainbow(nrow(datamatrix))

plotnametree1=paste(inputname, "_FstTREEunrEDIT_Ivaedited.png", sep="")
plotnametree2=paste(inputname, "_FstTREEfan.png", sep="")
plotnametree3=paste(inputname, "_FstTREEunr.svg", sep="")
plotnametree4=paste(inputname, "_FstTREEfan.svg", sep="")
titled=paste(info, " pair-wise Fst's tree-plot", sep="")

############################ EDITED PLOT PAPER #######################

ivapalette=c("deepskyblue2", "chocolate1", "lightblue1", "firebrick3", "gray66", "lavender", "darkslategray4", "darksalmon", "gold", "green4", "royalblue3", "darkorchid3", "wheat1", "hotpink4")
png(filename=plotnametree1, units="in", width=8, height=5, res=300)
plot(psic.tree, type="unr", tip.col=ivapalette, font=2, cex=0.9)
add.scale.bar(0.2,0, cex=0.7)
dev.off()

####################################################################
ivapalette=c("deepskyblue2", "chocolate1", "lightblue1", "firebrick3", "gray66", "lavender", "darkslategray4", "darksalmon", "gold", "green4", "royalblue3", "darkorchid3", "wheat1", "hotpink4")
ivaedit=c("deepskyblue2", "chocolate1", "skyblue", "firebrick3", "gray66", "thistle3", "darkslategray4", "darksalmon", "goldenrod", "green4", "royalblue3", "darkorchid3", "khaki3", "hotpink4")
png(filename=plotnametree1, units="in", width=8, height=5, res=300)
plot(psic.tree, type="unr", tip.col=ivaedit, font=2, cex=0.9)
add.scale.bar(0.2,0, cex=0.7)
dev.off()

####################################################################

png(filename=plotnametree1, units="in", width=8, height=5, res=600)
plot(psic.tree, type="unr", tip.col=mypalette, font=2, cex=0.9)
add.scale.bar(0.2,0, cex=0.5)
dev.off()

png(filename=plotnametree2, units="in", width=15, height=10, res=300)
plot(psic.tree, type="fan", tip.col=mypalette, font=2, main=titled, cex=0.7, cex.main= 1)
add.scale.bar(0.5,0, cex=0.5)
dev.off()
#add a bar to indicate the relationship between the length and the Fst value

#now svg
svg(filename=plotnametree3, width=15, height=10)
plot(psic.tree, type="unr", tip.col=mypalette, font=2, main=titled, cex=0.6, cex.main= 1)
add.scale.bar()
dev.off()

svg(filename=plotnametree4, width=15, height=10)
plot(psic.tree, type="fan", tip.col=mypalette, font=2, main=titled, cex=0.6, cex.main= 1)
add.scale.bar()
dev.off()



