#R sucks!

#LDA script, works at once with multiple grouping variables
#LINEAR DISCRIMINANT ANALYSIS it is in fact a regression analysis but it looks for the combination of variables that best explains the variance.
# the main difference from PCA is that it is not blind this takes your groups into account when calculating the equivalent to the PCs
# Good thing, it allows prediction

#This script is still under construction. It should work, but is not complete

#Change working directory if needed
getwd()			#Check
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/AdultsPKPM/Analysis/PCA_LDA")
getwd()


#					Packages					#
# RUN this if you experience any problem with the installed packages location
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

# You'll need to install this packages and load them in this order

library("MASS")
library(beepr) # to tell us when the analysis finished (turn on volume)
library(ggplot2)
library(scales)
library(ggpubr)
library(ggfortify)
library(gridExtra)
library(mvtnorm)
library(Momocs)
#install.packages("Momocs")

######################################



# The input file is:
#		tag, indiv, sex, pop, group, measure1, m2, m3, m4, m5 ... m16
#		A001, 001,	M,	  A,	MA,		04,		08,15,16,  23 ... 42
#		B002, 002,	F,	  B,	FB,		02,		71,82,81,  82 ... 84
#...etc
# First: columns with individual or non-usable data. In the example there is two ("tag" and "indiv") but there can be as many as desired, state how many below in "otherstuff" 
# Next: Columns with grouping variables, In the example there is three (sex, pop and group, from 3 to 5), state how many in "groupvar"
# Finally: Columns with quantitative variables, as many as desired, they will all taken into account if they are at the right of the other columns with other information.
# (In the example tabs and spaces are there only for the data to match the position of the header)

#I do NOT recommend to use LDA to analyse grouping variables with only two factors/groups.


############		SET UP THE INFORMATION ABOUT THE INPUT DATA FILE		###############


rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets

#input file
measurespod = read.table("Podadults2017residuals.csv", sep = ',', header =TRUE)	#read the data into a table, read headers
str(measurespod)															#chek



#Set up the number of initial columns with individual or non-relevant/usable data
otherstuff=2
#Set up the number of columns with grouping variables
groupvar = 3



############	BELOW THIS DON'T CHANGE ANYTHING UNLESS YOU KNOW WHAT YOU ARE DOING!	###########




coltags <- names(measurespod)
maxcol <- ncol (measurespod)								#number of columns
#maxcol=11
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



#Before runing the PCA check this:
cat ("Check that this is alright:\n\t", maxcol, "columns in TOTAL\n\t", otherstuff, "with individual/non-relevant data\n\t", groupvar, " columns with grouping variables:", groupnames, "\n\t", numvar, "variables:", varnames, "\n")







#################################################START ANALYSIS

analtype <- ("Linear Discriminant Analysis")
cat(" ", analtype, "----------------------------", "", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)
cat(numvar, " variables: ", file="podarcisLDAs_out.txt", sep="", append=TRUE)
cat(varnames, "\n", file="podarcisLDAs_out.txt", sep=" ", append=TRUE)


options(digits=5)


for (i in colgroup) {

	

	groups <- factor(measurespod[[i]])
	groups
	
	grouping <- names(measurespod[i])
	grouping

	podinfo = measurespod[c(i, firstvar:maxcol)]
	str(podinfo)

	tgroup <- measurespod[[i]]
	groupnames <- levels(tgroup)
	groupnames

	numgroups <- length(groupnames)
	numgroups

	tvar <- measurespod[c(firstvar:maxcol)]
	nvar <- ncol(tvar)
	nvar

	# BASIC LDA
	#Linear Discriminant Analysis: find the linear combinations of the original variables (the 13 chemical concentrations here) that gives the best possible separation between the groups


	cat("\n", "\n", "\n", grouping, "\n", "----", "\n", file="podarcisLDAs_out.txt", sep="", append=TRUE)
	cat("Levels:", groupnames, "\n", file="podarcisLDAs_out.txt", sep=" ", append=TRUE)
	phenopodLDA = lda(tgroup ~ ., data=tvar)
	phenopodLDA$counts		# number of indiv in each group
	phenopodLDA$means		# group specific means for each covariate
	#If wrong, we can specify them in the lda function with "prior="
	priortxt <- paste ("Prior probabilities of class membership (for each group/", grouping, ")", sep="") 
	cat("", "", priortxt, "", phenopodLDA$prior, "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)

	svdtxt <- paste ("Singular Values of Stadard Deviation for each LD")
	moretxt <-paste ("> Number of LDs = number of groups - 1 = ", numgroups, "-1 = ", numgroups-1, sep="")
	cat(moretxt, "", svdtxt, "", phenopodLDA$svd, "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE) #Singular values that gives the ratio of the between- and within-group standard deviations on the linear discriminant variables
	
	btw_gvar = phenopodLDA$svd^2/sum(phenopodLDA$svd^2)		# Use the singular values to compute the amount of the between-group variance that is explained by each linear discriminant
	# Check how many between group variance explains each linear discriminant
	btwgvar_txt <- ("Betwen-group variance explained by each LD")
	cat(btwgvar_txt, "", btw_gvar, "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)
	scaltxt <- ("Linear combination coefficients (scaling) for each LD")		# linear combination coefficients (scaling) for each linear discriminant (number of groups - 1 linear discriminants)
	

	scaling <- phenopodLDA$scaling
	scaling
	varrg <- (1:nvar)
	cat (scaltxt, "", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)
	for (asc in varrg)	{
		bsc=asc+numvar
		csc=bsc+numvar
		if (numgroups==2)	{ cat(varnames[asc], round(scaling[asc], digits=4), "\n", file="podarcisLDAs_out.txt", sep="\t", append=TRUE) }
		else if (numgroups==3)	{ cat(varnames[asc], round(scaling[asc], digits=4), round(scaling[bsc], digits=4), "\n", file="podarcisLDAs_out.txt", sep="\t", append=TRUE) }
		else if (numgroups>3)	{ cat(varnames[asc], round(scaling[asc], digits=4), round(scaling[bsc], digits=4), round(scaling[csc], digits=4), "\n", file="podarcisLDAs_out.txt", sep="\t", append=TRUE) }
		}
	cat("\n", "\n", "__________________________", "\n", "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)


	# LDA with Leave-one-out cross validation

	podlda = lda(tgroup ~ ., data=tvar, CV=TRUE)	#We assign the same prior probability to each group and set the CrossValidation
	head(podlda$class)
	head(podlda$posterior, numgroups-1)		#posterior probabilities for the classes
	
	cat("\n", " ---- Leave-One-Out Cross Validation analysis ---- ", "\n", file="podarcisLDAs_out.txt", sep="", append=TRUE)
	#cat(podlda$posterior, "\n\n", file="podarcisLDAs_out.txt", sep="", append=TRUE)
	posteriortxt <- ("posterior probabilities for each sample and group")
	posteriorprob <- podlda$posterior
	write.table (posteriorprob, file="podarcisLDAs_out.txt",sep="\t", append=TRUE)
	cat("\n", "\n", "__________________________", "\n", "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)
	
	
	
	# Plot the results of the LDA
	

	if (numgroups==2)	{
		palette (c("purple3", "goldenrod1", "red1", "black", "green3", "orange", "magenta", "gray", "blue"))
		mypalette <- c("purple3", "goldenrod1", "red1", "black", "green3", "orange", "magenta", "gray", "blue")
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
	
	str(phenopodLDA)
	scale_fill_manual(values=mypalette)
	scale_colour_manual(values=mypalette)
	scale_shape_manual(values=myshapes)

	if (numgroups==2)	{
		titlelda1 <- paste(grouping, "LDA_LDscores_pergroup_freepriorP", ".png", sep="")
		png(filename=titlelda1)						#Save to file
		plot(phenopodLDA, type="both", col="purple3")
		dev.off()
	} else if (numgroups>2)	{
		lda <- phenopodLDA
		plda <- predict(object = lda, newdata = podinfo)
		dataset = data.frame(groups = podinfo[,grouping], lda = plda$x)
		prop.lda = phenopodLDA$svd^2/sum(phenopodLDA$svd^2)
		
		titlelda2 <- paste(grouping, "LDA", "_freepriorP", ".png", sep="")
		Xlabel <- paste("LD1 (", round(100*(prop.lda[1]), digits=2), "%)", sep="")
		Ylabel <- paste("LD1 (", round(100*(prop.lda[2]), digits=2), "%)", sep="")
		png(filename=titlelda2, units="in", width=5, height=5, res=300)						#Save to file
		print(
		ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = groups, fill=groups, shape = groups), size = 2.5) + 
			labs(x = Xlabel, y = Ylabel) +
			coord_fixed()	+
			scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette) +
			scale_shape_manual(values=myshapes) +
			theme(panel.background=element_rect(fill = "white", colour = "black"), legend.key=element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
		dev.off()
	}
	
	
	
	# Now repeat everything but with the same prior probabilities for each group
	
	
	
	onesvect = rep(1, times=numgroups)
	
	
	cat("\n", " ---------- Forced equal prior probabilities for all groups -----------", "\n", file="podarcisLDAs_out.txt", sep="", append=TRUE)
	cat("\n", grouping, "\n", "------", "\n", file="podarcisLDAs_out.txt", sep="", append=TRUE)
	cat("Levels:", groupnames, "\n", file="podarcisLDAs_out.txt", sep=" ", append=TRUE)
	phenopodLDA = lda(tgroup ~ ., data=tvar, prior=c(onesvect)/numgroups)
	
	priortxt <- paste ("Prior probabilities of class membership for each group", sep="") 
	cat("", "", priortxt, "", phenopodLDA$prior, "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)
	phenopodLDA$means		# group specific means for each covariate
	


	svdtxt <- paste ("Singular Values of Stadard Deviation for each LD")
	moretxt <-paste ("> Number of LDs = number of groups - 1 = ", numgroups, "-1 = ", numgroups-1, sep="")
	cat(moretxt, "", svdtxt, "", phenopodLDA$svd, "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE) #Singular values that gives the ratio of the between- and within-group standard deviations on the linear discriminant variables
	
	btw_gvar = phenopodLDA$svd^2/sum(phenopodLDA$svd^2)		# Use the singular values to compute the amount of the between-group variance that is explained by each linear discriminant
	
	# Check how many between group variance explains each linear discriminant
	btwgvar_txt <- ("Betwen-group variance explained by each LD")
	cat(btwgvar_txt, "", btw_gvar, "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)
	
	scaltxt <- ("Linear combination coefficients (scaling) for each LD")		# linear combination coefficients (scaling) for each linear discriminant (number of groups - 1 linear discriminants)
	scaling <- phenopodLDA$scaling
	scaling
	varrg <- (1:nvar)
	cat (scaltxt, "", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)
	for (asc in varrg)	{
		bsc=asc+numvar
		csc=bsc+numvar
		if (numgroups==2)	{ cat(varnames[asc], scaling[asc], "\n", file="podarcisLDAs_out.txt", sep="\t", append=TRUE) }
		else if (numgroups==3)	{ cat(varnames[asc], scaling[asc], scaling[bsc], "\n", file="podarcisLDAs_out.txt", sep="\t", append=TRUE) }
		else if (numgroups>3)	{ cat(varnames[asc], scaling[asc], scaling[bsc], scaling[csc], "\n", file="podarcisLDAs_out.txt", sep="\t", append=TRUE) }
		}
	cat("\n", "\n", "_____________________________________", "\n", "\n", "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)


	# LDA with Leave-one-out cross validation


	podlda = lda(tgroup ~ ., data=tvar, prior=c(onesvect)/numgroups, CV=TRUE)	#We assign the same prior probability to each group and set the CrossValidation
	head(podlda$class)
	head(podlda$posterior, numgroups-1)		#posterior probabilities for the classes
	
	cat("\n", " ---- Leave-One-Out Cross Validation analysis ---- ", "\n", file="podarcisLDAs_out.txt", sep="", append=TRUE)
	#cat(podlda$posterior, "\n\n", file="podarcisLDAs_out.txt", sep="", append=TRUE)
	posteriortxt <- ("posterior probabilities for each sample and group")
	posteriorprob <- podlda$posterior
	write.table (posteriorprob, file="podarcisLDAs_out.txt",sep="\t", append=TRUE)
	cat("\n", "\n", "__________________________", "\n", "\n", file="podarcisLDAs_out.txt", sep="\n", append=TRUE)

	

	if (numgroups==2)	{
		palette (c("purple3", "goldenrod1", "red1", "black", "green3", "orange", "magenta", "gray", "blue"))
		mypalette <- c("purple3", "goldenrod1", "red1", "black", "green3", "orange", "magenta", "gray", "blue")
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
	
	#PLOT
	if (numgroups==2)	{
		titlelda1 <- paste(grouping, "LDA_LDscores_pergroup_equalpriorP", ".png", sep="")
		png(filename=titlelda1)						#Save to file
		plot(phenopodLDA, type="both", col="purple3")
		dev.off()
	} else if (numgroups>2)	{
		lda <- phenopodLDA
		plda <- predict(object = lda, newdata = podinfo)
		dataset = data.frame(groups = podinfo[,grouping], lda = plda$x)
		prop.lda = phenopodLDA$svd^2/sum(phenopodLDA$svd^2)
		titlelda2 <- paste(grouping, "LDA", "_equalpriorP", ".png", sep="")
		Xlabel <- paste("LD1 (", round(100*(prop.lda[1]), digits=2), "%)", sep="")
		Ylabel <- paste("LD1 (", round(100*(prop.lda[2]), digits=2), "%)", sep="")
		
		str(lda)
		
		png(filename=titlelda2, units="in", width=5, height=5, res=300)						#Save to file
		print(
		ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, fill = groups, shape = groups), size = 2.5) + 
			labs(x = Xlabel, y = Ylabel) +
			coord_fixed()	+
			scale_fill_manual(values=mypalette) + scale_colour_manual(values=mypalette) +
			scale_shape_manual(values=myshapes)
			)
		dev.off()
	
	
	
	
	
	################################ COMPARE BOTH LDA and PCA

		#PCA
		
	if (numgroups==2)	{
		palette (c("purple3", "goldenrod1", "red1", "black", "green3", "orange", "magenta", "gray", "blue"))
		mypalette <- c("purple3", "goldenrod1", "red1", "black", "green3", "orange", "magenta", "gray", "blue")
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
	
		
		
		
		phenopodPCA <- prcomp(scale(varinfo))
		summary(phenopodPCA)

		pca <- phenopodPCA

		prop.pca = phenopodPCA$sdev^2/sum(phenopodPCA$sdev^2)

		prop.lda = phenopodLDA$svd^2/sum(phenopodLDA$svd^2)

		plda <- predict(object = lda,
						newdata = podinfo)

		dataset = data.frame(groups = podinfo[,grouping],
							 pca = pca$x, lda = plda$x)


		p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, fill= groups, shape = groups), size = 2.5) + 
		  labs(x = paste("LD1 (", round(100*(prop.lda[1]), digits=2), "%)", sep=""),
			   y = paste("LD1 (", round(100*(prop.lda[2]), digits=2), "%)", sep="")) +
			   scale_fill_manual(values=mypalette) + scale_shape_manual(values=myshapes)

		p2 <- ggplot(dataset) + geom_point(aes(pca.PC1, pca.PC2, fill=groups, shape = groups), size = 2.5) +
		  labs(x = paste("PC1 (", round(100*(prop.pca[1]), digits=2), "%)", sep=""),
			   y = paste("PC2 (", round(100*(prop.pca[2]), digits=2), "%)", sep="")) +
			   scale_fill_manual(values=mypalette) + scale_shape_manual(values=myshapes)

		comparetitle <- paste(grouping, "_LDAvsPCA", ".png", sep="")
		png(filename=comparetitle, units="in", width=5, height=5, res=300)						#Save to file
		grid.arrange(p1, p2)
		dev.off()
	}
}




beep(3)


########################################


