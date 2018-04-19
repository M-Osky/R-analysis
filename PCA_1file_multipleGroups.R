#R sucks!

#PCA
#Script to perform a PRINCIPAL COMPONENT ANALYSIS with different methods and see the scores and other information.
#It works with multiple grouping variables at the same time (does the analysis for each one)
#You just need to identify the grouping variables (usually cathegorical as population, sex, etc) in the data and run it
#It uses Iva's script to perform an ANOVa from the GLM scores of the first PC 

#load the packages, set the input file, check that it identifies correctly the variables, and run the rest of the script
#plots and results will be saved in the working directory

#Check working directory
getwd()
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/F1crosses/Analysis/Multivariant")
getwd()


#					Packages					#

# =====================================================================
# For compatibility with Rscript.exe: 
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

#source("https://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")
#a

library(pcaMethods)
library(beepr) # to tell us when the analysis finished (turn on volume)
library(ggplot2) # install+load installr
library(ggfortify) # install+load installr
library (dplyr)
library(ellipse)
#library(psych) # install+load installr



########################################



# The input file is:
#		tag, indiv, sex, pop, group, measure1, m2, m3, m4, m5 ... mn
#		A001, 001,	M,	  A,	MA,		04,		08,15,16,  23 ... 42
#		B002, 002,	F,	  B,	FB,		02,		71,82,81,  82 ... 84
#...etc
# First: columns with individual or non-usable data. In the example there is two ("tag" and "indiv") but there can be as many as desired, further below you will need to set this as "otherstuff" 
# Next: Columns with grouping variables, In the example there is three (sex, pop and group, from 3 to 5), further below you'll need to set how many in "groupvar"
# Finally: Columns with quantitative variables, as many as desired, they will all taken into account if they are at the right of the other columns with other information.
# (In the example tabs and spaces are there only for the data to match the position of the header)






rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets


############		SET UP THE INFORMATION ABOUT THE INPUT DATA FILE		###############

#input file
measurespod = read.table("Loglength2017f1_2nout.csv", sep = ',', header =TRUE)	#read the data into a table, read headers
str(measurespod)				#chek the data



################################# OPTIONAL: Analise only a subset#############
#str(measurespod)	#check sataset   										###
# # Choose on of the subsets and perform the analysis  						###
# # Subset 1  																###
#Mmother <- measurespod[which (measurespod$Mod =="M"),]   					###
#myvars <- names(Mmother) %in% c("Mod", "Cross")  							###
#subseted <- Mmother[!myvars]												###
# # Subset 2  																###
#Mfather <- measurespod[which (measurespod$Fad =="M"),]   					###
#myvars <- names(Mfather) %in% c("Fad", "Cross")   							###
#subseted <- Mfather[!myvars]												###
# # Subset 3  																###
#Kmother <- measurespod[which (measurespod$Mod =="K"),]  					###
#myvars <- names(Kmother) %in% c("Mod", "Cross")   							###
#subseted <- Kmother[!myvars]												###
# # Subset 4  																###
#Kfather <- measurespod[which (measurespod$Fad =="K"),]   					###
#myvars <- names(Kfather) %in% c("Fad", "Cross")   							###
#subseted <- Kfather[!myvars]												###
																			###
# #For any subset chosen: check new dataset and rename it     				###
#str(subseted)																###
#measurespod <- subseted  													###
#str(measurespod)															###
##############################################################################


##				IMPORTANT			##
#Set up the number of initial columns with individual or non-relevant/usable data
otherstuff=2
#Set up the number of columns with grouping variables
groupvar = 3



########################	BELOW THIS DON'T CHANGE ANYTHING UNLESS YOU KNOW WHAT YOU ARE DOING!	###########


coltags <- names(measurespod)
maxcol <- ncol (measurespod)								#number of columns

lastgrup = otherstuff + groupvar							#column with the last grouping variable
firstvar = lastgrup + 1										# column with the first variable
firstgrup = otherstuff+1
numvar = maxcol - lastgrup
varinfo <- measurespod[c(firstvar:maxcol)]									#take a subset including the data of the measured variables

groupinfo <- measurespod[c(firstgrup:lastgrup)]									#take a subset with the grouping variables to test
groupnames <- names (groupinfo)
colgroup = c(firstgrup:lastgrup)
colvariab = c(firstvar:maxcol)
varnames <- names(varinfo)

#CHECK BEFORE PROCEEDING
cat ("Check that this is alright:\n\t", maxcol, "columns in TOTAL\n\t", otherstuff, "with individual/non-relevant data\n\t", groupvar, " columns with grouping variables:", groupnames, "\n\t", numvar, "variables:", varnames, "\n")


#### IF EVERYTHING WAS ALRIGHT RUN THE REST OF THE SCRIPT
#				START PCA


for (i in colgroup) {

	if (groupvar == 1) {
		colgroup = 1
		i = firstgrup
	}
	i

	groups <- factor(measurespod[[i]])
	group <- (measurespod[[i]])

	#extract the levels/group names
	tags <- levels(measurespod[[i]]) 
	tags

	grouping <- names(measurespod[i])
	grouping
	podinfo = measurespod[c(i, firstvar:maxcol)]
	str(podinfo)
	
	numgroups <- length(tags)

	
	#colors
	if (numgroups<3)	{
		palette (c("black", "red1", "goldenrod1", "purple3", "green3", "orange", "magenta", "gray", "blue"))
		mypalette <- c("black", "red1", "goldenrod1", "purple3", "green3", "orange", "magenta", "gray", "blue")
		myshapes=c(22, 21, 23:25)
	} else if (numgroups==3)	{
		palette (c("red1", "goldenrod1", "purple3", "black", "green3", "orange", "magenta", "blue"))
		mypalette <- c("red1", "goldenrod1", "purple3", "black", "green3", "orange", "magenta", "blue")
		myshapes=c(22, 21, 25:23)
	} else if (numgroups>3) {	
		palette (c("darkorange3", "darkgreen", "orange", "chartreuse3", "black", "red1", "goldenrod1", "purple3", "blue", "magenta"))
		mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3", "black", "red1", "goldenrod1", "purple3", "blue", "magenta")
		myshapes=c(22,22,21,21,23:25)
	}

	

	#See variables by pairs
	titlepca <- paste(grouping, "PCA", "_pair_distrib", ".png", sep="")
	png(	filename=titlepca, units="in", width=5, height=5, res=300)
	pairs(varinfo, col = groups, upper.panel = NULL, pch = 19, cex = 0.2)
	legend("topright", bty = "n", legend = tags, pch = 19, pt.cex=1, text.width=0.05 ,col = mypalette, xpd = T, cex = 0.7, y.intersp = 0.8)
	dev.off() # clear the format from the previous plot


	#Do the default R PCA
	phenopodPCA <- prcomp(scale(varinfo))

	#decide how many PC to take
	pcomp <- ("PCA's PCs")
	out <- capture.output(summary(phenopodPCA))							#save results to file
	cat("", "----", grouping, "----", "\n", pcomp, "", out, file="podarcisPCAs_out.txt", sep="\n", append=TRUE)
	cat("\n--", "", file="podarcisPCAs_out.txt",sep="\n", append=TRUE)
	out = ""														#Empty

	titlepca <- paste(grouping, "PCA", "_PCs1and2", ".png", sep="")
	png(filename=titlepca, units="in", width=5, height=5, res=300)
	screeplot(phenopodPCA, type="lines")
	dev.off() # clear the format from the previous plot

	#To see the loadings, that is how each variable affects each component
	
	
	loadings <- phenopodPCA$rotation
	cat("", "Loadings", "", sep="\n",file="podarcisPCAs_out.txt", append=TRUE)
	write.table (format (loadings, digits=4),file="podarcisPCAs_out.txt", sep="\t", append=TRUE)
	
	
	
	#See the PCA then
	titlepca <- paste(grouping, "PCA", "_basic", ".png", sep="")
	png(filename=titlepca, units="in", width=5, height=5, res=300)
	plot(phenopodPCA$x[,1:2], col = groups, pch = 19)
	title(main=grouping)
	#legend ("bottomright", legend = tags, pch=19, col = c("black", "red1", "goldenrod1", "purple3", "green3", "blue", "magenta", "gray"))
	dev.off() # clear the format from the previous plot
	
	
	####################################################################
	
	####significance aov####
	
	#Extract scores
	pcascores<-phenopodPCA$x
	
	#extract the information of individual tags and grouping variables from original file
	allinfo <- measurespod[c(1:lastgrup)]

	#save scores and info
	write.table(pcascores,file="pca_scores.csv",sep=",", append=FALSE)
	write.table(allinfo,file="tags.csv",sep=",", append=FALSE)
	
	#join information from both
	importscores = read.table("pca_scores.csv", sep = ',', header =FALSE, skip=1)
	importgroups = read.table("tags.csv", sep = ',', header =FALSE, skip=1)
	allimport <- bind_cols(importgroups, importscores)
	
	#Depurate
	numPCs <- ncol(pcascores)
	firstscore <- (lastgrup+3)
	lastscore <- (numPCs+lastgrup+1)
	lastinfocol <- (lastgrup+1)
	
	pca_scores <- allimport[c(2:lastinfocol, firstscore:lastscore)]
	str(pca_scores)
	
	#select PC and grouping
	pc1 <- allimport[[firstscore]]
	group <- pca_scores[[i]]
	pc1
	#GLM
	mod1<-glm(pc1~group, data=pca_scores)
	mod1
	aov(mod1)->aov_pc1
	aov_pc1
	summary(aov_pc1)->aov_pod_pca_status
	
	#save
	out <- capture.output(aov_pod_pca_status)
	cat("", "--------", "", "", "GLM of the first PC scores by group", "", out, file="podarcisPCAs_out.txt", "\n", sep="\n", append=TRUE)
	out=""
	
	#done
	
	
	
	
	
	
	
	###############################################################
	#fancier
	palettecolor=mypalette
	scale_fill_manual(values=mypalette)
	scale_colour_manual(values=mypalette)
	
	titlepca <- paste(grouping, "_PCA", "_autoplot", ".png", sep="")
	prop.pca = phenopodPCA$sdev^2/sum(phenopodPCA$sdev^2)
	Xlabel <- paste("PC1 (", round(100*(prop.pca[1]), digits=2), "%)", sep="")
	Ylabel <- paste("PC2 (", round(100*(prop.pca[2]), digits=2), "%)", sep="")

	
	png(filename=titlepca, units="in", width=5, height=5, res=300)						#Save to file
	print(autoplot(phenopodPCA, data=podinfo) +
	labs(x = Xlabel, y = Ylabel) + coord_fixed(ratio = 1) +
	#theme(axis.line = element_line(colour = "gray70"), panel.background = element_blank(), legend.key=element_blank(), legend.background = element_rect(linetype="solid", colour ="gray70")) +
	theme_bw()+
	geom_point(aes(shape = group, color=group, fill = group), size = 2.5) +
	scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette) +
	scale_shape_manual(values=myshapes))
	dev.off()
	

	

	
	
	#Now lets try a different method
	podphenPCAmethods <- pca(podinfo[,-1], scale = "uv", center = T, nPcs = 2, method = "svd")


	#Let's check the variance explained by each component
	str(podphenPCAmethods) # slots are marked with @
	podphenPCAmethods@R2


	titlepca2 <- paste(grouping, "_loadingsPCs", ".png", sep="")
	png(filename=titlepca2, units="in", width=5, height=5, res=300)
	slplot(podphenPCAmethods, scoresLoadings = c(F,T), scol = groups, pch=19)
	#legend ("top", bty = "n",legend = tags, pch=1, col = c("black", "red1", "goldenrod1", "purple3", "green3", "blue", "magenta", "gray"))
	dev.off() # clear the format from the previous plot


	#To identify the samples
	pca2 <- prcomp(varinfo, center=TRUE, scale=TRUE)
	prop.pca = pca2$sdev^2/sum(pca2$sdev^2)
	dataset = data.frame(groups = podinfo[,grouping], pca2 = pca2$x)
	titlepca2 <- paste("PCA_ggplot2", ".png", sep="")

	Indiv <- measurespod[[1]]
	png(filename=titlepca2, units="in", width=10, height=10, res=500)
	print(ggplot(dataset, aes(x=pca2.PC1, y=pca2.PC2, label=Indiv)) +
	geom_point(aes(pca2.PC1, pca2.PC2, colour = groups), size = 0.5) +
	geom_text(aes(label=Indiv), size=1, hjust=0.5, vjust=0, nudge_y=0.09) +
	labs(x = paste("PC1 (", round(100*(prop.pca[1]), digits=2), "%)", sep=""),
	y = paste("PC2 (", round(100*(prop.pca[2]), digits=2), "%)", sep="")) +
	scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette) +
	coord_fixed(2))
	dev.off()

	
	cat("", "----------------------------------", "", file="podarcisPCAs_out.txt", "\n", sep="\n", append=TRUE)

	
	
}


beep(3)




