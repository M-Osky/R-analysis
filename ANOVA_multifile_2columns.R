#ANOVAs

# Do simple ANOVAs from multiple files

#Change working directorysetwd
setwd("D:/Dropbox/MOSKY/CURRO/DATA/Meassurements/Analysis/ANOVA")
getwd()			#Check

library(beepr)

inputfiles <- list.files(pattern='csv')							#take all csv files from wd
inputfiles														#just checking

#input file must be
#		#Comment or Title
#		group,measurements
#		A,220
#		B,250
#...etc

for (inputfile in inputfiles) {
	#grab data
	measurespod = read.table(inputfile, sep = ',', header =TRUE)	#read the data into a table, read headers
	measurespod														#chek
	
	variablename <- colnames (measurespod)
	
	variable <- measurespod[[2]]									#data from measurements
	variable
	group <- measurespod[[1]]										#group to analise
	group
	
	aovpod = aov(variable~group,data=measurespod)					#do the analysis of variance
	
	variablename[2]
	print(model.tables(aovpod,"means"),digits=4)					#report the means and the number of subjects/cell
	
	out <- capture.output(summary(aovpod))							#save results to file
	cat(inputfile, out, file="data_measures_anovas.txt", sep="\n", append=TRUE)
	cat(" \n--------------------", " "," ", file="data_measures_anovas.txt",sep="\n", append=TRUE)
	out = ""														#Empty
	imagename <- paste(inputfile, ".png", sep="")
	png(filename=imagename)											#save plot
	boxplot(variable~group,data=measurespod)						#graphical summary
	dev.off()
}

cat(" \n\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n\n", " ", " ", " ", "\n", file="data_measures_anovas.txt", sep="\n", append=TRUE)

beep(3)



