#ANOVAs from 1 file

#Do all the ANOVAs from all the possible combinations of group vs measurements in one input file

#SO FAR I CAN'T MAKE IT WORK... I'LL KEEP TRYING IN A NEAR FUTURE

#Change working directorysetwd
setwd("D:/Dropbox/MOSKY/CURRO/DATA/Meassurements/Analysis/ANOVA")
getwd()			#Check

library(beepr)

measurespod = read.table("allphenotlog.csv", sep = ',', header =TRUE)	#read the data into a table, read headers
measurespod														#chek

#input file must be
#		#Comment or Title
#		tag,sex,pop,group,measurement1,measurement2,measurement3,measurement4,measurement5,measurement6...measurement15
#		001,M,A,MA,04,08,15,16,23,42...n
#		002,F,B,FB,13,12,24,42,15,18...n
#...etc


variablenames <- colnames (measurespod)
variablenames

variables <- measurespod[c(5:15)]									#data from measurements
variables
group <- measurespod[[4]]										#group1 to analise
group
pop <- measurespod[[3]]										#group2 to analise
pop
sex <- measurespod[[2]]										#group3 to analise
sex


for (variable in variables) {
	#grab data

	
	aovpod = aov(variable~group,data=measurespod)					#do the analysis of variance
	
	variablename[2]
	print(model.tables(aovpod,"means"),digits=4)					#report the means and the number of subjects/cell
	
	out <- capture.output(summary(aovpod))							#save results to file
	cat(inputfile, out, file="data_mesures_anovas.txt", sep="\n", append=TRUE)
	cat(" \n--------------------", " "," ", file="data_mesures_anovas.txt",sep="\n", append=TRUE)
	out = ""														#Empty
	imagename <- paste(inputfile, ".png", sep="")
	png(filename=imagename)											#save plot
	boxplot(variable~group,data=measurespod)						#graphical summary
	dev.off()
}

cat(" \n\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n\n", " ", " ", " ", "\n", file="data_mesures_anovas.txt", sep="\n", append=TRUE)

beep(3)

