#R sucks!


# quick R script to calculate average, sicula585x39k_m05R7r6h6DP420_1p0p_0m0pSD and confidence intervals from multiple files


rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets

#Check working directory
getwd()
#setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/F1crosses/2017/mine/Parent_F1resid190402/subset_2ndround/Fem/male_var")
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/PopGenet/genpopR/sicula596x39k/sicula585x39K_noutfix/pi/NOLDHWe_max19")



#set extension or pattern to filter out files
ext="*.pi"

#set column separator
separator = "\t"

#set column in which the data is at
columndata=3

#set if file has headers
columnlabels=TRUE

#number of characters to take from the beginning of file names in order to ID the data
poplength=3



############################		OPEN INPUT DATA FILE		##############################
filelist <- list.files(path = ".", pattern = ext)		  #change extension if needed 
filelist

#check if file agrees
examplefile <- filelist[1]
measurespod = read.table(examplefile, sep = separator, header = columnlabels)	#read the data into a table, read headers
str(measurespod)				#chek the data
head(measurespod)


#run if it does


#save first part of file name and keep the rest (without extension) for the output file
numfiles = length(filelist)
myrange = c(1:numfiles)


for (filenum in myrange) {
	inputname=filelist[filenum]
	inputname
	fileID = substr(inputname, 1, poplength)
	inputname
	noext = gsub('^(.*)\\..*$', '\\1', inputname)
	longword = nchar(noext)
	outID = substr(inputname, poplength+1, longword)
	outID
	fileout = paste(outID, ".txt", sep="")



	#input file
	onefile = read.table(inputname, sep = separator, header = columnlabels)	#read the data into a table, read headers
	alldata = as.numeric(onefile[,columndata])
	
	#summary indexes
	average= mean(alldata)
	standev = sd(alldata)
	numb = length(alldata)
	CIup = average + (1.96 * (standev / (sqrt(length(alldata)))))
	CIdown = average - (1.96 * (standev / (sqrt(length(alldata)))))
	
	quantiles = quantile(alldata)
	minimum <- as.numeric(as.vector(quantiles[1]))
	firstq <- as.numeric(as.vector(quantiles[2]))
	medianq <- as.numeric(as.vector(quantiles[3]))
	thirdq <- as.numeric(as.vector(quantiles[4]))
	maximum <- as.numeric(as.vector(quantiles[5]))
	
	#save in a dataframe
	if (filenum == 1) {
		allindexes = data.frame(fileID, numb, minimum, firstq, medianq, thirdq, maximum, average, standev, CIdown, CIup)
	} else {
		allindexes[filenum,] <- data.frame(fileID, numb, minimum, firstq, medianq, thirdq, maximum, average, standev, CIdown, CIup)
	}
	
	
}

names(allindexes) <- c("ID", "n", "min", "Q1", "Median", "Q3", "Max", "Mean", "SD", "CI95down", "CI95up")
head(allindexes)
write.table (format (allindexes, digits=5),file=fileout, sep="\t", quote=FALSE, row.names=FALSE)

cat("\n\nFinished!\n\n")
















