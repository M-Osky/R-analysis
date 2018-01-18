#ANOVAs

#Change working directorysetwd
setwd("D:/Dropbox/MOSKY/CURRO/podarcis/Meassurements/Analysis/Multivariate")
getwd()			#Check

source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")


#library(ade4) # multivariate analysis
#library(grid) # has the viewport function
#library(sinkr)

library(beepr) # to tell us when the analysis finished (turn on volume)
library(pcaMethods)
library(ggplot2)
library(ggfortify)




measurespod = read.table("phenodataout.csv", sep = ',', header =TRUE)	#read the data into a table, read headers

# The input file is:
#		tag, sex,pop,group,measure1,m2,m3,m4,m5 ... m15
#		001, M,  A,  MA,		04,08,15,16,23 ... 42
#		002, F,  B,  FB,		02,71,82,81,82 ... 84
# tabs and spaces are there only for the data to match the position of the header

measurespod															#chek
podinfo <- measurespod[c(4, 6:16)]									#take a subset (all)

podinfo																#check

varnames<- colnames (podinfo)
varnames

# The first column corresponds to the classes
groups <- factor(podinfo[[1]])
groups

#extract the levels/group names
tags <- levels(podinfo[[1]]) 
tags

#colors
palette (c("black", "red1", "goldenrod1", "purple4", "green3", "blue", "magenta", "gray"))

pairs(podinfo[,-1], col = groups, upper.panel = NULL, pch = 19, cex = 0.5)
legend("topright", bty = "n", legend = tags, pch = 19, col = c("black", "red1", "goldenrod1", "purple4"), xpd = T, cex = 1.5, y.intersp = 0.8)

dev.off() # clear the format from the previous plot


#Do the default R PCA
phenopodPCA <- prcomp(scale(podinfo[,-1]))

#decide how many PC to take
pcomp <- ("Principal Components Variance")
out <- capture.output(summary(phenopodPCA))							#save results to file
cat(pcomp, out, file="dataPCAs.txt", sep="\n", append=TRUE)
cat(" \n----", " ", file="dataPCAs.txt",sep="\n", append=TRUE)
out = ""														#Empty

screeplot(phenopodPCA, type="lines")

dev.off() # clear the format from the previous plot

#To see the loadings, that is how each variable affects each component
PCloadings <- ("Principal Components Loads for each Variable")
out1 <- capture.output(phenopodPCA$rotation[,1])							#save results to file
out2 <- capture.output(phenopodPCA$rotation[,2])							#save results to file
cat(PCloadings, out1, out2, file="dataPCAs.txt", sep="\n", append=TRUE)
cat(" \n--------------------", " "," ", file="dataPCAs.txt",sep="\n", append=TRUE)
out1 = ""														#Empty
out2 = ""														#Empty



#See the PCA then
plot(phenopodPCA$x[,1:2], col = groups, pch = 19, labels=TRUE, label.size=3)
legend ("bottomright", legend = tags, pch=19, col = c("black", "red1", "goldenrod1", "purple4"))

dev.off() # clear the format from the previous plot

#fancier
autoplot(phenopodPCA, data=podinfo, colour='group', label=TRUE)

dev.off() # clear the format from the previous plot




#Now lets try a different method
podphenPCAmethods <- pca(podinfo[,-1], scale = "uv", center = T, nPcs = 2, method = "svd")


#Let's check the variance explained by each component
str(podphenPCAmethods) # slots are marked with @
podphenPCAmethods@R2

slplot(podphenPCAmethods, scoresLoadings = c(T,T), scol = groups, pch=19 )
legend ("top", bty = "n",legend = tags, pch=1, col = c("black", "red1", "goldenrod1", "purple4"))

dev.off() # clear the format from the previous plot



beep(3)



###################################################################################
###################################################################################
###################################################################################
###################################################################################



#Let's try a different approach, I'm still working on this, so it may be not completely right.
#LDA
#Linear Discriminant Analysis: find the linear combinations of the original variables that gives the best possible separation between the groups

library("MASS")
podphenLDA = lda(podinfo$group ~ podinfo$HH + podinfo$HL + podinfo$HW + podinfo$SL + podinfo$LJL + podinfo$LJO + podinfo$SVL + podinfo$LTH + podinfo$BLOT + podinfo$WBL + podinfo$BITE)
# In this example "measure1" is "HH", "measure2" is "HL", etc. 

podphenLDA					#check
podphenLDA$scaling[,1]		#check the loadings of the first discriminant function

#calculate the value of each discriminant function, by substituting the variables’ values into the linear combination for the discriminant function
stdedvariables <- as.data.frame(scale(podinfo[,-1])) # standardise the variables
valuesLDA <- predict(podphenLDA, stdedvariables)


calcWithinGroupsVariance <- function(variable,groupvariable) {
	# find out how many values the group variable can take
	groupvariable2 <- as.factor(groupvariable[[1]])
	levels <- levels(groupvariable2)
	numlevels <- length(levels)
	# get the mean and standard deviation for each group:
	numtotal <- 0
	denomtotal <- 0
	for (i in 1:numlevels) {
		leveli <- levels[i]
		levelidata <- variable[groupvariable==leveli,]
		levelilength <- length(levelidata)
		# get the standard deviation for group i:
		sdi <- sd(levelidata)
		numi <- (levelilength - 1)*(sdi * sdi)
		denomi <- levelilength
		numtotal <- numtotal + numi
		denomtotal <- denomtotal + denomi
	}
	# calculate the within-groups variance
	Vw <- numtotal / (denomtotal - numlevels)
	return(Vw)
}


calcBetweenGroupsVariance <- function(variable,groupvariable) {
	# find out how many values the group variable can take
	groupvariable2 <- as.factor(groupvariable[[1]])
	levels <- levels(groupvariable2)
	numlevels <- length(levels)
	# calculate the overall grand mean:
	grandmean <- mean(variable)
	# get the mean and standard deviation for each group:
	numtotal <- 0
	denomtotal <- 0
	for (i in 1:numlevels) {
		leveli <- levels[i]
		levelidata <- variable[groupvariable==leveli,]
		levelilength <- length(levelidata)
		# get the mean and standard deviation for group i:
		meani <- mean(levelidata)
		sdi <- sd(levelidata)
		numi <- levelilength * ((meani - grandmean)^2)
		denomi <- levelilength
		numtotal <- numtotal + numi
		denomtotal <- denomtotal + denomi
	}
	# calculate the between-groups variance
	Vb <- numtotal / (numlevels - 1)
	Vb <- Vb[[1]]
	return(Vb)
}


calcSeparations <- function(variables,groupvariable) {
	# find out how many variables we have
	variables <- as.data.frame(variables)
	numvariables <- length(variables)
	# find the variable names
	variablenames <- colnames(variables)
	# calculate the separation for each variable
	for (i in 1:numvariables) {
		variablei <- variables[i]
		variablename <- variablenames[i]
		Vw <- calcWithinGroupsVariance(variablei, groupvariable)
		Vb <- calcBetweenGroupsVariance(variablei, groupvariable)
		sep <- Vb/Vw
		print(paste("variable",variablename,"Vw=",Vw,"Vb=",Vb,"separation=",sep))
	}
}

podinfo[1]

#calculate the separations achieved by the two linear discriminant functions
calcSeparations(valuesLDA$x,podinfo[1])


















beep(3)

