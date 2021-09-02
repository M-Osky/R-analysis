#R sucks

#Script to calculate linear regression in multiple files, implementing now to plot the data and do quadratic regression


#Check working directory
getwd()
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/nat_pop_Iva/pheno-env/jun21_RDAnew/correldces/regress_right")




rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets

# setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements")
# getwd()

##		FILE LIST, INPUT PATH AND PATTERN
filelist <- list.files(path = ".", pattern = ".csv$")		  #Change "." for the path of a directory with your data files "folder/subdirectory/here"; You can also make pattern more specific to filter out some files: "_final_version.csv"
filelist		  #check that all the fileas are alright

#open one file to check
wrk <- read.delim(filelist[[1]], sep = ',', header =TRUE)		  #comma separated file with column headers
str(wrk)		  #check structure of the file
head(wrk)		  #check file


#Individual tags if any must be in the left most column
#Number of Columns with Stuff (to ignore) -including individual tags and factor/grouping variables, before the first numerical variable to analyze-
otherstuff=3

# "otherstuff" must be the same for all files, all the other columns must be variables to correlate
# does not matter how many variables there are in each file
# Last column must be the independent variable to correlate against (in all files)

# you want to log10 transfrom variables before regression?
logtrans = "yes"
logtrans = "no"

# you want also to try a quadratic regression?
quadratic = "yes"
quadratic = "no"

# run and check
coltags <- names(wrk)
maxcol <- ncol (wrk)								#number of columns
lastgrup = otherstuff 							#column with the last grouping variable
firstvar = lastgrup + 1										# column with the first variable
allvar = maxcol - lastgrup
varnames = coltags[c(firstvar:(maxcol-1))]
varnum = length(varnames)
indepvar = coltags[[maxcol]]

cat ("\n\nFile checked has", varnum, "variables:", varnames, "\n will be correlated against ", indepvar, "\n\n")






####################################################################################
################## RUN EVERYTHING BELLOW IF IT SEEMS ALRIGHT ########################
####################################################################################


rawtime = Sys.time()
hour<- sub(".* ", "", rawtime)
onlydate<-sub(" .*", "", rawtime)
printime <- gsub(":", "-", hour)

numfiles = length(filelist)		  # extract how many files there are
myrange = c(1:numfiles)			  #range from 1 to number of files

file1 = paste("summary_", onlydate, "_", printime, "_RegressionOut.txt", sep="")
cat ( "Regression analysis summary table\n", "File_name\ttype\tvar\tIndepVar\ttransform?\tR^2\tslope\tpval_x\tpval_x^2", sep="\n", file=file1, append=FALSE)

options(scipen = 9999)
#k=1

for (k in myrange) {
	csvfile <- filelist[[k]]		  #choose one of the file names in each loop
	# OPEN INPUT FILE
	cat("Processing ", csvfile, "\n")		  #Print on screen
	measurespod <- read.delim(filelist[[k]], sep = ',', header =TRUE)
	
	coltags <- names(measurespod)
	#coltags
	maxcol <- ncol (measurespod)								#number of columns
	lastgrup = otherstuff							#column with the last grouping variable
	firstvar = lastgrup + 1										# column with the first variable
	allvar = maxcol - lastgrup
	varnames = coltags[c(firstvar:(maxcol-1))]
	varrange = c(firstvar:(maxcol-1))
	#varrange
	varnum = length(varnames)
	indepvar = coltags[[maxcol]]
	xcorr = as.vector(measurespod[[maxcol]])
	#xcorr
	#measurespod[[maxcol]]
	file2 <- gsub(".csv", "_LinearResiduals.csv", csvfile)
	#cat ( coltags[1], varnames, indepvar, sep="\t", file=file2, append=FALSE)
	
	tableout <- data.frame(measurespod[c(1:otherstuff)])
	#tableout

	#extract the columns with the independent variable for the correlations
	#varrange
	cat("\nlinear regression", sep="\n")
	
	for (i in varrange) {
		
		#i = 14
		varname <- coltags[i]
		
		cat(varname, sep="\n")
		
		
		correlatetable <- measurespod[c(1,i,maxcol)]
		#correlatetable <- as.data.frame(Indiv=factor(measurespod[[1]]), numvar=double(measurespod[[i]]), indepvar=(measurespod[[maxcol]]))
		names(correlatetable) <- c("IDs", "NumVar", "IndepVar")
		#correlatetable
		
		#LINEAR REGRESSION
		linreg <- lm(NumVar~IndepVar, data=correlatetable)
		
		#save results
		regress <- summary(linreg)
		out <- capture.output(regress)
		r2 <- round(regress$r.squared, digits=4)
		#regress
		pval <- regress$coefficients[8]
		intercept <- round(regress$coefficients[1], digits=4)
		slope <- round(regress$coefficients[2], digits=4)
		#cat ("", paste(title1, ":  R^2 = ", r2, ";  p-value = ", pval, sep=""), "", out, "\n##################\n\n", sep="\n", file="linregr_out.txt", append=TRUE)
		cat(csvfile, "linear", varname, indepvar, "raw", r2, slope, pval, "NA", sep="\t", file=file1, append=TRUE)
		cat("", sep="\n", file=file1, append=TRUE)
		
		
		
		
		#save graphic
		
		textx=max(correlatetable$IndepVar)
		texty=min(correlatetable$NumVar)
		
		endtag = paste("_", varname, ".png", sep="")
		outname <- gsub(".csv", endtag, csvfile)
		outfileGraphic <- paste("scatter_", outname, sep="") 
		png(outfileGraphic, width = 1300, height = 1300, res=144)
		plot(correlatetable$IndepVar, correlatetable$NumVar, pch=16, xlab=indepvar, ylab=varname, )
		abline(lm(correlatetable$NumVar~correlatetable$IndepVar))
		text(textx, texty, paste("R^2 = ", r2, " (p-val=", round(pval, digits=4), ")", sep=""), pos=2, cex=0.7)
		dev.off()
		
		
		
		#residuals
		residuallist <- as.data.frame(round(linreg$residuals, digits=3))
		#residuallist
		
		#residuallist
		tableout <-as.data.frame(c(tableout, residuallist))

		#tableout
		
		rm(correlatetable)
		rm(residuallist)

	}
	
	#str(tableout)
	#tableout <-as.data.frame(c(tableout, xcorr))
	
	#finalcol = length(varrange)+1
	tablefin = tableout
	#str(finalcol)
	tablefin$indepvar <- xcorr
	#finalcol
	#tablefin <- as.data.frame(tableout[c(1:finalcol)], measurespod[[maxcol]])
	#tablefin <- data.frame()
	
	#names(finalcol) <- c("IDs", coltags[c(firstvar:(maxcol))])
	names(tablefin) <- coltags
	#names(tableout) <- c("IDs", coltags[c(firstvar:(maxcol-1))])
	#print(tablefin)
	
	#str(tablefin)
	
	write.csv(tablefin, quote=FALSE, file = file2, row.names=FALSE)
	#write.csv(finalcol, quote=FALSE, file = file2, row.names=FALSE)
	
	
	#cat(tableout, sep="\t", 
	
	rm(tableout)
	
	cat ("\nLinear regression saved for all.\n")
	
	
	
	# QUADRATIC REGRESSION
	
	if (quadratic == "yes") {
		
		tableout <- data.frame(measurespod[c(1:otherstuff)])
		cat("\nQuadratic regression", sep="\n")
		file3 <- gsub(".csv", "_QuadratResiduals.csv", csvfile)
		
		for (i in varrange) {
			
			#i = 14
			varname <- coltags[i]
			
			cat(varname, sep="\n")
			
			
			correlatetable <- measurespod[c(1,i,maxcol)]
			names(correlatetable) <- c("IDs", "NumVar", "IndepVar")
			
			#power of two:
			correlatetable$Indep2 <- correlatetable$IndepVar^2
			
			quareg <- lm(NumVar~IndepVar+Indep2, data=correlatetable)
			
			
			#save results
			regress <- summary(quareg)
			out <- capture.output(regress)
			r2 <- round(regress$r.squared, digits=4)
			pval1 <- regress$coefficients[11]
			pval2 <- regress$coefficients[12]
			intercept <- round(regress$coefficients[1], digits=4)
			#slope <- round(regress$coefficients[2], digits=4)
			cat(csvfile, "quadratic", varname, indepvar, "raw", r2, "NA", pval1, pval2, sep="\t", file=file1, append=TRUE)
			cat("", sep="\n", file=file1, append=TRUE)
			
			#residuals
			residuallist <- as.data.frame(round(quareg$residuals, digits=3))
			
			#residuallist
			tableout <-as.data.frame(c(tableout, residuallist))

			#tableout
			
			rm(correlatetable)
			rm(residuallist)

		}
		
		tablefin = tableout
		tablefin$indepvar <- xcorr
		names(tablefin) <- coltags
		write.csv(tablefin, quote=FALSE, file = file3, row.names=FALSE)
		
		cat ("\nQuadratic regression saved for all.\n")
		rm(tableout)
	
	}
	
	cat ("\nDone!\n\n")
	
}

cat ("\n\n\nAll residual calculated from raw data!\n\n\n")

#i=15
if (logtrans == "yes") {
	for (k in myrange) {
		csvfile <- filelist[[k]]		  #choose one of the file names in each loop
		# OPEN INPUT FILE
		cat("Now log 10 transforming ", csvfile, "\n")		  #Print on screen
		measurespod <- read.delim(filelist[[k]], sep = ',', header =TRUE)
		
		coltags <- names(measurespod)
		#coltags
		maxcol <- ncol (measurespod)								#number of columns
		lastgrup = otherstuff							#column with the last grouping variable
		firstvar = lastgrup + 1										# column with the first variable
		allvar = maxcol - lastgrup
		varnames = coltags[c(firstvar:(maxcol-1))]
		varrange = c(firstvar:(maxcol-1))
		#varrange
		varnum = length(varnames)
		indepvar = coltags[[maxcol]]
		xcorr = as.vector(measurespod[[maxcol]])
		#xcorr
		#measurespod[[maxcol]]
		
		#generate output file
		
		if (length(grep("raw", csvfile)) > 0) {
			newcsv <- (gsub("^(.*?)raw(.*?)$", "\\1log10\\2", csvfile))
		} else if (length(grep("RAW", csvfile)) > 0) {
			newcsv <- (gsub("^(.*?)raw(.*?)$", "\\1LOG10\\2", csvfile))
		} else {
			newcsv <- (gsub("^(.*)(.csv)$", "\\1_log10\\2", csvfile))
		}
		
		
		
		file2 <- gsub(".csv", "_LinearResiduals.csv", newcsv)
		#cat ( coltags[1], varnames, indepvar, sep="\t", file=file2, append=FALSE)
		
		tableout <- data.frame(measurespod[c(1:otherstuff)])
		#tableout

		#extract the columns with the independent variable for the correlations
		#varrange
		cat("\nlinear regression", sep="\n")
		
		for (i in varrange) {
			
			#i = 14
			varname <- paste("log10_", coltags[i], sep="")
			
			indepx = paste("log10_", coltags[maxcol], sep="")
			
			cat(varname, sep="\n")
			
			
			
			logvar = log10(measurespod[i])
			logind = log10(measurespod[maxcol])
			idtags = measurespod[1]
			
			correlatetable <- as.data.frame(c(idtags, logvar, logind))
			
			#correlatetable <- as.data.frame(Indiv=factor(measurespod[[1]]), numvar=double(measurespod[[i]]), indepvar=(measurespod[[maxcol]]))
			names(correlatetable) <- c("IDs", "NumVar", "IndepVar")
			#correlatetable
			
			#LINEAR REGRESSION
			linreg <- lm(NumVar~IndepVar, data=correlatetable)

			#save results
			regress <- summary(linreg)
			out <- capture.output(regress)
			r2 <- round(regress$r.squared, digits=4)
			#regress
			pval <- regress$coefficients[8]
			intercept <- round(regress$coefficients[1], digits=4)
			slope <- round(regress$coefficients[2], digits=4)
			#cat ("", paste(title1, ":  R^2 = ", r2, ";  p-value = ", pval, sep=""), "", out, "\n##################\n\n", sep="\n", file="linregr_out.txt", append=TRUE)
			cat(csvfile, "linear", varname, indepx, "log10", r2, slope, pval, "NA", sep="\t", file=file1, append=TRUE)
			cat("", sep="\n", file=file1, append=TRUE)
			
			
			#save graphic
			
			textx=max(correlatetable$IndepVar)
			texty=min(correlatetable$NumVar)
			
			endtag = paste("_", varname, ".png", sep="")
			outname <- gsub(".csv", endtag, newcsv)
			outfileGraphic <- paste("scatter_", outname, sep="") 
			png(outfileGraphic, width = 1300, height = 1300, res=144)
			plot(correlatetable$IndepVar, correlatetable$NumVar, pch=16, xlab=indepx, ylab=varname, )
			abline(lm(correlatetable$NumVar~correlatetable$IndepVar))
			text(textx, texty, paste("R^2 = ", r2, " (p-val=", round(pval, digits=4), ")", sep=""), pos=2, cex=0.7)
			dev.off()
			
			
			#residuals
			residuallist <- as.data.frame(round(linreg$residuals, digits=3))
			#residuallist
			
			#residuallist
			tableout <-as.data.frame(c(tableout, residuallist))

			#tableout
			
			rm(correlatetable)
			rm(residuallist)

		}
		
		#str(tableout)
		#tableout <-as.data.frame(c(tableout, xcorr))
		
		#finalcol = length(varrange)+1
		tablefin = tableout
		#str(finalcol)
		tablefin$indepvar <- xcorr
		#finalcol
		#tablefin <- as.data.frame(tableout[c(1:finalcol)], measurespod[[maxcol]])
		#tablefin <- data.frame()
		coltags
		#names(finalcol) <- c("IDs", coltags[c(firstvar:(maxcol))])
		names(tablefin) <- coltags
		#names(tableout) <- c("IDs", coltags[c(firstvar:(maxcol-1))])
		#print(tablefin)
		
		#str(tablefin)
		
		write.csv(tablefin, quote=FALSE, file = file2, row.names=FALSE)
		#write.csv(finalcol, quote=FALSE, file = file2, row.names=FALSE)
		
		
		#cat(tableout, sep="\t", 
		
		rm(tableout)
		
		cat ("\nLinear regression saved for all.\n")
		
		
		
		# QUADRATIC REGRESSION
		
		if (quadratic == "yes") {
			
			tableout <- data.frame(measurespod[c(1:otherstuff)])
			cat("\nQuadratic regression", sep="\n")
			file3 <- gsub(".csv", "_QuadratResiduals.csv", newcsv)
			
			for (i in varrange) {
				
				#i = 14
				varname <- coltags[i]
				
				varname <- paste("log10_", coltags[i], sep="")
				
				indepx = paste("log10_", coltags[maxcol], sep="")
				
				
				cat(varname, sep="\n")
				
				
				#transform and check if there are negative values
				logvar = log10(measurespod[i])
				logind = log10(measurespod[maxcol])
				idtags = measurespod[1]
				
				logind
				
				minim <- min(logind)
				
				
				#if there are negative values in independent variable make them all positive
				if(minim < -0.1) {
					toadd=0.1
					negmin=abs(minim)+toadd		   #minimum negative value, this will added to all the values, so there is no negative values
					#negmin
					logvar = logvar+negmin
					logind = logind+negmin
					
				} else if (minim < 0) {
					toadd=0.01
					negmin=abs(minim)+toadd		   #minimum negative value, this will added to all the values, so there is no negative values
					#negmin
					logvar = logvar+negmin
					logind = logind+negmin
					
				} else { negmin = 0 }
				
				
				
				correlatetable <- as.data.frame(c(idtags, logvar, logind))
				names(correlatetable) <- c("IDs", "NumVar", "IndepVar")
				
				#power of two:
				correlatetable$Indep2 <- correlatetable$IndepVar^2
				
				quareg <- lm(NumVar~IndepVar+Indep2, data=correlatetable)
				
				
				#save results
				regress <- summary(quareg)
				out <- capture.output(regress)
				r2 <- round(regress$r.squared, digits=4)
				pval1 <- regress$coefficients[11]
				pval2 <- regress$coefficients[12]
				intercept <- round(regress$coefficients[1], digits=4)
				#slope <- round(regress$coefficients[2], digits=4)
				cat(csvfile, "quadratic", varname, indepx, "log10", r2, "NA", pval1, pval2, sep="\t", file=file1, append=TRUE)
				cat("", sep="\n", file=file1, append=TRUE)
				
				
				
				#residuals
				residuallist <- as.data.frame(round(quareg$residuals, digits=3))
				
				#residuallist
				tableout <-as.data.frame(c(tableout, residuallist))

				#tableout
				
				rm(correlatetable)
				rm(residuallist)

			}
			
			tablefin = tableout
			tablefin$indepvar <- xcorr
			names(tablefin) <- coltags
			write.csv(tablefin, quote=FALSE, file = file3, row.names=FALSE)
			
			rm(tableout)
			cat ("\nQuadratic regression saved for all.\n")
		
		
		}
		
		
		
		cat ("\nAll done!\n\n\n")
		
	}
}


cat ("\n\n      FINISHED!\n\n\n")



