#R sucks!

#Pst as an aproximationn to Qst (phenotypic distance) and comparison with Fst.
#using Bayesian GLMM approach.



# You'll need to install this packages and load them in this order

#install.packages("MCMCglmm")
#library("dyplr")
library("MCMCglmm")
#matrix, coda, ape

######################################



# The input file is:
#		tag, indiv, sex, year, pop, measure1, m2, m3, m4, m5 ... m16
#		A001, 001,	M,	  A,	MA,		04,		08,15,16,  23 ... 42
#		B002, 002,	F,	  B,	FB,		02,		71,82,81,  82 ... 84
#...etc
# Columns with variables to be used as fixed or random effects, in the example there are three (sex, year, and pop, from 3 to 5)
# Columns with quantitative variables, as many as desired (measure1 to m16).
#ALL FILES ANALYSED MUST HAVE columns with no-usable data, columns with fixed/random effects, and first column with numerical response variable in the same position
#but number of numerical response variables does not need to be the same as far as they are all to the right of first column with numerical response variable parsed
# (In the example tabs and spaces are there only for the data to match the position of the header)



############		SET UP THE INFORMATION ABOUT THE INPUT DATA FILE		###############


rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets


#Change working directory if needed
getwd()			#Check
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/Pst/MCMCglmm")
getwd()

#define an ending all files must share
#filtetail="*.csv"
filtetail="*.csv"
#

#
filelist <- list.files(path = ".", pattern = filtetail)		  #pattern = ".csv" to read all the files with the same format
filelist
#check one file
filecheck = read.table(filelist[1], sep = ',', header =TRUE)	#read the data into a table, read headers
str(filecheck)															#chek
head(filecheck)															#chek
#


#Set up the position of columns with variables to be used as fixed effects (account for their effect)
# any position as far as they are at the left of quantitative variables to be used as response variables
#minimum total number 1 maximum total number 3
#if more than one parse between quotation marks, space separated, no commas.
#efixed = "3 4"
#efixed = 5
efixed = 4

#Set up the position of the column with the variable to be used as random effect (measure differences among its factors)
# any position as far as it is at the left of quantitative variables to be used as response variables
#only analysis of one random effect implemented
#erandom = 5
erandom = 3

#position of first column with a numerical response variable to analyse
#from this column till the right end all columns must be numerical response variables to analyse
firstvar=6


#MCMC
#iterations (default 13k)
iterations=15000
#take a sample every how many iterations?
samplefreq=10
#define burn-in (number of the iterations to be discarded)
discard=500
#prior c / h^2 if unknnow is recommended start with 1
ch2ratio = 1
#Fst global value
Fst=0.62718
#Fst lower 95& CI limit
upFst=0.62914




##################		RUN THIS CHUNK OF CODE FIRST TO CHECK IF ALL SEEMS ALRIGHT		###################
################	BELOW THIS DON'T CHANGE ANYTHING UNLESS YOU KNOW WHAT YOU ARE DOING!	###########

#function to count matches on a string
string.counter<-function(strings, pattern){  
  counts<-NULL
  for(i in 1:length(strings)){
    counts[i]<-length(attr(gregexpr(pattern,strings[i])[[1]], "match.length")[attr(gregexpr(pattern,strings[i])[[1]], "match.length")>0])
  }
return(counts)
}
#if more than one column prsed, split them
spacefix = string.counter(strings=efixed, pattern=" ")
if (spacefix>0) {
	colfix= as.numeric(strsplit(efixed, " ")[[1]])
} else {
	colfix= as.numeric(efixed)
}
spacernd = string.counter(strings=erandom, pattern=" ")
if (spacernd>0) {
	colrnd= as.numeric(strsplit(erandom, " ")[[1]])
} else {
	colrnd= as.numeric(erandom)
}
numfix=length(colfix)
numrnd=length(colrnd)
#count and take info from dataset
maxcol <- ncol (filecheck)	#number of columns
numvar = maxcol - (firstvar - 1)
varinfo <- filecheck[c(firstvar:maxcol)]									#take a subset including the data of the measured variables
rndinfo <- filecheck[colrnd]									#take a subset with the random effect variables to test
fixinfo <- filecheck[colfix]									#take a subset with the fixed effect variables to test
rndnames <- names (rndinfo)
fixnames <- names (fixinfo)
varnames <- names(varinfo)

#Before runing check this:
cat ("Check that this is alright:\n\t", maxcol, "columns in TOTAL
\t", numfix, " column(s) with variables analysed as fixed effects:", fixnames, "
\t", numrnd, " column(s) with variables to analyse as random effects:", rndnames,"\n and ", numvar, "numerical response variables:", varnames, "\n")
#######################################################################




				##########################
				##### RUN ANALYSIS!!! ####
				##########################


#lets create  a function to calculate the model

if (numfix > 3 || numfix == 0 || numrnd != 1 || Fst==666 || upFst==666 ) {
			cat ("\n\n\n\n\tERROR!\n\tNumber of variables parsed as \"random effects\" (", numrnd, ") or as \"fixed effects\" (", numfix, "), may be wrong.\n\tOr may be you forgot to parse a proper value for Fst?.\n\tCheck the script help: Pst_mcmc_evoleco.R --help\n\n\n\n")
			stop()
			break()
}

getmode <- function(v) {
	uniqv <- unique(v)
	uniqv[which.max(tabulate(match(v, uniqv)))]
}

myrange = c(1:length(filelist))

filenum=1		   #debug

for (filenum in myrange) {
	
	
	
	inputname=filelist[filenum]
	filename = gsub('^(.+?)_.*$', '\\1', inputname)
	
	#prepare outputs
	output_analysis = paste(filename, "_GLMMresults.txt", sep="")
	tablename = paste(filename, "_tablevals.txt", sep="")
	tableprint = data.frame(variable=character(0), fixed=character(0), random=character(0), saved_chains=integer(0), Vb_mean=double(0), Vb_mode=double(0), Vb_low=double(0), Vb_upp=double(0), Vw_mean=double(0), Vw_mode=double(0), Vw_low=double(0), Vw_upp=double(0), Pst=double(0), Fst=double(0), selec=character(0), ch2_crit=double(), stringsAsFactors = FALSE)
	
	#read
	measurespod = read.table(inputname, sep = ',', header =TRUE)	#read the data into a table, read headers
	
	cat("\n\nAnalysing ", filename, "\n\n")
	cat("#################################################################################\n ###    ", inputname, "    ### \n#################################################################################\n-f ", Fst, " -t ", upFst, " -e \"", efixed, "\" -r \"", erandom, "\" -v ", firstvar, " -i ", iterations, " -b ", discard, " -s ", samplefreq, " -x \"", filtetail, "\" -c ", ch2ratio, "\n\n", sep="", file=output_analysis)
	
	#count and take info from dataset
	maxcol <- ncol (measurespod)	#number of columns
	numvar = maxcol - (firstvar - 1)
	varrng = c(firstvar:maxcol)
	varinfo <- measurespod[varrng]									#take a subset including the data of the measured variables
	rndinfo <- measurespod[colrnd]									#take a subset with the random effect variables to test
	fixinfo <- measurespod[colfix]									#take a subset with the fixed effect variables to test
	rndnames <- names (rndinfo)
	fixnames <- names (fixinfo)
	allnames <- names(measurespod)
	
	#varrng
	
	#i=7		#debug
	
	
	skip=0
	chosen = allnames[varrng]
	cat("Analysing ", chosen, "\n")
	
	chosenplus <-paste(chosen, collapse="+")
	
	#str(measurespod)		   #debug
	
	#select only the columns to analyse
	#try adding the individual ID "id"
	#dataset <- measurespod[c(1, varrng, colfix, colrnd)]
	
	
	#I have an error with more than one response variables, it seems you need to account for residual covariation of the fixed variable (trait)
	#then got an error about " ill-conditioned G/R structure" and on internet seems like that could be fixed by scaling data
	dataset <- data.frame(measurespod[1], scale(measurespod[c(varrng)]), measurespod[c(colfix, colrnd)])
	head(dataset)
	
	keepnames <- names(dataset)
	printfix <-paste(allnames[colfix], collapse="+")
	printrnd <-paste(allnames[colrnd], collapse="+")
	
	cat("Performing MCMCglmm with", printfix, "as fixed effect(s), and", printrnd, "as random effect.\n\n", sep=" ")
	cat("\n\n-----------------------------------------------------------------------------------------------------------\n", chosenplus, "\nGLMM analysis with", allnames[colfix], "as fixed effects, and", allnames[colrnd], "as random effect.\n-----------------------------------------------------------------------------------------------------------\n\n", sep=" ", file=output_analysis, append=TRUE)
	
	
	families = rep("gaussian", times=3)
	#families = rep("gaussian", times=numvar)
	#HLgth+HWdth+HHgth+LwJaL+LwJaO+SnLgh+ILLgh+HLLgh+FLLgh+BHgth+BWdth+BiteF+BMs+SVLgh
	
	# TEST 1!!!!
	if (numfix == 1) {
	
	
		newnames <- c("id", chosen, "fixed_effect", "random_effect")
		names(dataset) <- newnames
		head(dataset)
		
		randomtest <- MCMCglmm(HLgth+HWdth+HHgth ~ fixed_effect, random = ~random_effect, family=families, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		#Error in MCMCglmm(HLgth + HWdth + HHgth ~ fixed_effect, random = ~random_effect,  :
		#R-structure miss-specified: each residual must be unique to a data point
		randomtest <- MCMCglmm(HLgth+HWdth+HHgth ~ fixed_effect, random = ~random_effect, family=families, rcov = ~us(trait):units, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		#crashes
		
		randomtest <- MCMCglmm(HLgth+HWdth+HHgth ~ fixed_effect, random = ~us(random_effect):id, family=families, rcov = ~us(trait):units, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		#Error in MCMCglmm(HLgth + HWdth + HHgth ~ fixed_effect, random = ~us(random_effect):id,  : 
		#Singular G/R structure: use proper priors
		
		randomtest <- MCMCglmm(HLgth+HWdth+HHgth ~ fixed_effect, random = ~us(trait):random_effect, family=families, rcov = ~us(trait):id, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		#MCMC iteration = 0
		#Error in MCMCglmm(HLgth + HWdth + HHgth ~ fixed_effect, random = ~us(trait):random_effect,  : 
		#  R-structure 2 is ill-conditioned: use proper priors if you haven't or rescale data if you have
		
		randomtest <- MCMCglmm(HLgth+HWdth+HHgth ~ fixed_effect, random = ~random_effect, family=families, rcov = ~us(trait):id, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		
		
		randomtest <- MCMCglmm(cbind(HLgth, HWdth) ~ fixed_effect, random = ~us(random_effect):id, family=c("gaussian", "gaussian"), rcov = ~us(fixed_effect):units, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		randomtest <- MCMCglmm(HLgth+HWdth ~ fixed_effect-1, random = ~us(random_effect):id, family=c("gaussian", "gaussian"), rcov = ~us(fixed_effect):units, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		randomtest <- MCMCglmm(HLgth+HWdth ~ fixed_effect-1, random = ~us(random_effect):id, family=c("gaussian", "gaussian"), rcov = ~us(trait):units, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		randomtest <- MCMCglmm(HLgth+HWdth ~ fixed_effect-1, random = ~us(trait):id, family=c("gaussian", "gaussian"), rcov = ~us(trait):units, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		randomtest <- MCMCglmm(HLgth+HWdth+HHgth ~ fixed_effect, random = ~random_effect, family=c("gaussian", "gaussian", "gaussian"), rcov = ~us(trait):id, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
		
	} else if (numfix == 2) {
		newnames <- c(chosen, "fixed_effect1", "fixed_effect2", "random_effect")
		names(dataset) <- c(chosen, "fixed_effect1", "fixed_effect2", "random_effect")
		randomtest <- MCMCglmm(response_var ~ fixed_effect1+fixed_effect2, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	} else if (numfix == 3) {
		newnames <- c(chosen, "fixed_effect1", "fixed_effect2", "fixed_effect3", "random_effect")
		names(dataset) <- c("response_var", "fixed_effect1", "fixed_effect2", "fixed_effect3", "random_effect")
		randomtest <- MCMCglmm(response_var ~ fixed_effect1+fixed_effect2+fixed_effect3, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	}
	
	head(dataset)
	
	
	
	
	#str(dataset)
	cat("variable names:\n", newnames, "\n", keepnames, "\n\n", sep="\t", file=output_analysis, append=TRUE)
	summary(randomtest)
	#get relevant infor from summary
	infotest=capture.output(summary(randomtest))
	realsize= as.numeric(strsplit(infotest[4], "\\s+")[[1]][5])
	randomsize=as.numeric(strsplit(infotest[11], "\\s+")[[1]][5])
	cat("First replicate -> ", sep="", file=output_analysis, append=TRUE)
	
	
	
	if (randomsize*2 < realsize && randomsize < 1000) {
		cat("\n\nEffective sample size is not high (", randomsize, ") and smaller than MCMC sample size (", realsize, ")\nFirst run will not be used for the final estimation.\n\n", sep="")
		cat("Effective sample size is not high (", randomsize, ") and smaller than MCMC sample size (", realsize, ")\nThis run will not be used for the final estimation.\n\n", sep="", file=output_analysis, append=TRUE)
		skip=1
	} else {
		cat("Effective sample size is big and not much smaller than MCMC sample size: Good levels of autocorrelation\n\n", sep="", file=output_analysis, append=TRUE)
		allvarw1 = as.vector(randomtest$VCV[,"units"])
		allvarb1 = as.vector(randomtest$VCV[,"random_effect"])
		skip=0
	}
	
	cat(infotest, sep="\n", file=output_analysis, append=TRUE)
	
	
	
	# TEST 2!!!!
	if (numfix == 1) {
		randomtest <- MCMCglmm(response_var ~ fixed_effect, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	} else if (numfix == 2) {
		randomtest <- MCMCglmm(response_var ~ fixed_effect1+fixed_effect2, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	} else if (numfix == 3) {
		randomtest <- MCMCglmm(response_var ~ fixed_effect1+fixed_effect2+fixed_effect3, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	}
	summary(randomtest)
	#get relevant infor from summary
	infotest=capture.output(summary(randomtest))
	realsize= as.numeric(strsplit(infotest[4], "\\s+")[[1]][5])
	randomsize=as.numeric(strsplit(infotest[11], "\\s+")[[1]][5])
	cat("\n\n\n-----------------------------------------------------------------------------------------------------------\nSecond replicate -> ", sep="", file=output_analysis, append=TRUE)
	
	if (randomsize*2 < realsize && randomsize < 1000) {
		cat("\n\nEffective sample size is not high (", randomsize, ") and smaller than MCMC sample size (", realsize, ")\nSecond run will not be used for the final estimation.\n\n", sep="")
		cat("Effective sample size is not high (", randomsize, ") and smaller than MCMC sample size (", realsize, ")\nThis run will not be used for the final estimation.\n\n", sep="", file=output_analysis, append=TRUE)
		allvarb = as.vector("")
		allvarb = as.vector("")
		skip=skip+2
	} else {
		cat("Effective sample size is big and not much smaller than MCMC sample size: Good levels of autocorrelation\n\n", sep="", file=output_analysis, append=TRUE)
		allvarw2 = as.vector(randomtest$VCV[,"units"])
		allvarb2 = as.vector(randomtest$VCV[,"random_effect"])
	}
	cat(infotest, sep="\n", file=output_analysis, append=TRUE)
	
	
	# TEST 3!!!!
	if (numfix == 1) {
		randomtest <- MCMCglmm(response_var ~ fixed_effect, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	} else if (numfix == 2) {
		randomtest <- MCMCglmm(response_var ~ fixed_effect1+fixed_effect2, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	} else if (numfix == 3) {
		randomtest <- MCMCglmm(response_var ~ fixed_effect1+fixed_effect2+fixed_effect3, random = ~random_effect, data = dataset, nitt=iterations, thin=samplefreq, burnin=discard, pr=TRUE, pl=TRUE, )
	}
	summary(randomtest)
	#get relevant infor from summary
	infotest=capture.output(summary(randomtest))
	realsize= as.numeric(strsplit(infotest[4], "\\s+")[[1]][5])
	randomsize=as.numeric(strsplit(infotest[11], "\\s+")[[1]][5])
	cat("\n\n\n-----------------------------------------------------------------------------------------------------------\nThird replicate ->", sep="", file=output_analysis, append=TRUE)
	
	if (randomsize*2 < realsize && randomsize < 1000) {
		cat("\n\nEffective sample size is not high (", randomsize, ") and smaller than MCMC sample size (", realsize, ")\nThird run will not be used for the final estimation.\n\n", sep="")
		cat("Effective sample size is not high (", randomsize, ") and smaller than MCMC sample size (", realsize, ")\nThis run will not be used for the final estimation.\n\n", sep="", file=output_analysis, append=TRUE)
		if (skip==3) {
			allrunb = 0
			allrunw = 0
			skip = 42
		} else if (skip == 1) {
			allrunb = data.frame(allvarb2)
			allrunw = data.frame(allvarw2)
		} else if (skip == 2) {
			allrunb = data.frame(allvarb1)
			allrunw = data.frame(allvarw1)
		} else if (skip == 0) {
			allrunb = data.frame(allvarb1, allvarb2)
			allrunw = data.frame(allvarw1, allvarw2)
		}
		
		
		
	} else {
		cat("Effective sample size is big and not much smaller than MCMC sample size: Good levels of autocorrelation\n\n", sep="", file=output_analysis, append=TRUE)
		allvarw3 = as.vector(randomtest$VCV[,"units"])
		allvarb3 = as.vector(randomtest$VCV[,"random_effect"])
		if (skip==3) {
			allrunb = data.frame(allvarb3)
			allrunw = data.frame(allvarw3)
		} else if (skip == 1) {
			allrunb = data.frame(allvarb2, allvarb3)
			allrunw = data.frame(allvarw2, allvarw3)
		} else if (skip == 2) {
			allrunb = data.frame(allvarb1, allvarb3)
			allrunw = data.frame(allvarw1, allvarw3)
		} else if (skip == 0) {
			allrunb = data.frame(allvarb1, allvarb2, allvarb3)
			allrunw = data.frame(allvarw1, allvarw2, allvarw3)
		}
	}
	
	
	cat(infotest, sep="\n", file=output_analysis, append=TRUE)
	
	
	
	
	###now check convergence
	
	#head(allrunb)
	#head(allrunw)
	runnames = names(allrunb)
	#runnames
	numchains = length(runnames)
	#head(allrunb)
	
	if (skip != 42 && numchains > 1) {
		names(allrunb) <- rep("postVarB", times=numchains)
		chainrg = c(1:numchains)
		chainrg
		n=1
		for (n in chainrg) {
			thischain = allrunb[[n]]
			chainw = allrunw[[n]]
			namechain = runnames[n]
			#namechain
			#thischain
			
			
			
			rawchain<-mcmc(thischain,thin=samplefreq)
			if(n==1) {
				#save the data
				allvarb = as.vector(thischain)
				allvarw = as.vector(chainw)
				#create a new list
				chainlist = list("cacafuti"=rawchain)
				names(chainlist) <- namechain
				#str(chainlist)
			} else { 
				#join chains inn the same vector for later
				allvarb = as.vector(c(allvarb, as.vector(thischain)))
				allvarw = as.vector(c(allvarw, as.vector(chainw)))
				#add chains to the list
				chainlist[[namechain]] = rawchain
			}
			
		}
		
		#str(chainlist)
		combined <- mcmc.list(chainlist)
		#str(combined)
		#Check convergence
		convergtest <- gelman.diag(combined)
		
		#extract values
		point <- convergtest$psrf[1]
		pointprint <- round(point, 4)
		upper <- convergtest$psrf[2]
		upperprint <- round(upper, 4)
		
		
		if(upper<1.1 && point<1.1) {
			cat("\nLog post within and between chain variances are quite similar, this is a sign of convergence!\n")
			cat("\n\n\nGelman and Rubin's convergence diagnostic:.\n\nLog post within (", pointprint, ") and between chain variances (", upperprint, ") are quite similar.\nBoth are below 1.1 and similar to 1, this is a sign of convergence!\n\n", sep="", file=output_analysis, append=TRUE)
		} else { 
			cat("\nWarning:Log post within and between chain variance are different, this is a sign of NO convergence!\n")
			cat("\n\n\nGelman and Rubin's convergence diagnostic:.\n\nWARNING: Log post within (", pointprint, ") and between chain variances (", upperprint, ") are disimilar or above/different from 1.1.\nThis is a sign of NO convergence!\n\n", sep="", file=output_analysis, append=TRUE)
		}
		
		names(allrunb) <- runnames
		#cat(allrunb,sep="\t")
		saveresults=1
		
	} else if (skip != 42 && numchains == 1) {
		
		cat("\nOnly one succesful run, no convergence analysis possible.\n")
		cat("\nOnly one succesful run, no convergence analysis possible.\n", sep="", file=output_analysis, append=TRUE)
		allvarw = as.vector(allrunw[1])
		allvarb = as.vector(allrunb[1])
		cat(allvarn)
		cat("\n\n")
		cat(allvarb)
		cat("\n\n")
		saveresults=1
		
	} else if (skip == 42) { 
		
		cat("\nThere was NO succesful runs for ", chosen , "\nTry to run it again with different parameters...\n\n\n", sep="")
		cat("\nThere was NO succesful runs for ", chosen , "\nTry to run it again with different parameters...\n\n\n", sep="", file=output_analysis, append=TRUE)
		saveresults=0
		numchains=0
	} else {
		cat("\nWell, this is akward... This message was not supposed to appear.\nThere is some problem with variable numchains (", numchains,"), and parameter skip (", skip, ").\nPlease debug me...\n\n\n", sep="")
		cat("\nWell, this is akward... This message was not supposed to appear.\nThere is some problem with variable numchains (", numchains,"), and parameter skip (", skip, ").\nPlease debug me...\n\n\n", sep="", file=output_analysis, append=TRUE)
	}
	
	
	
	
	if (saveresults==1) {
		
		#extract the mode (you will need first to join the output from all runs)
		modevarw = getmode(allvarw)
		meanw = mean(allvarw)
		modevarb = getmode(allvarb)
		meanb=mean(allvarb)
		modevarb
		modevarw
		
		
		#extract the quantiles of the distribution  (you will need first to join the output from all runs)
		sortedvalb = sort(allvarb)
		credible_int_b = quantile(sortedvalb, probs=c(0.05,0.95))
		minb = credible_int_b[[1]]
		maxb = credible_int_b[[2]]
		sortedvalw = sort(allvarw)
		credible_int_w = quantile(sortedvalw, probs=c(0.05,0.95))
		minw=credible_int_w[[1]]
		maxw=credible_int_w[[2]]
		
		#calculate Pst and c /h2 crit
		Pst = (ch2ratio*modevarb) / ((ch2ratio*modevarb)+(2*modevarw))
		Pst=round(Pst, digits=6)
		if (Pst>Fst) { select="divergent"} else if (Pst < Fst) { select="stabilizing"  } else { select="neutral" } 
		ch2crit = (2*maxb*upFst)/(minb*(1-upFst))
		
		#generate an output table
		inputline = nrow(tableprint)+1
		tableprint[inputline,1] <- chosen
		tableprint[inputline,2] <- as.character(printfix)
		tableprint[inputline,3] <- as.character(printrnd)
		tableprint[inputline,4] <- numchains
		tableprint[inputline,5] <- round(meanb, digits=6)
		tableprint[inputline,6] <- round(modevarb, digits=6)
		tableprint[inputline,7] <- round(minb, digits=6)
		tableprint[inputline,8] <- round(maxb, digits=6)
		tableprint[inputline,9] <- round(meanw, digits=6)
		tableprint[inputline,10] <- round(modevarw, digits=6)
		tableprint[inputline,11] <- round(minw, digits=6)
		tableprint[inputline,12] <- round(maxw, digits=6)
		tableprint[inputline,13] <- round(Pst, digits=6)
		tableprint[inputline,14] <- round(Fst, digits=6)
		tableprint[inputline,15] <- as.character(select)
		tableprint[inputline,16] <- round(ch2crit, digits=6)
		#tableprint[nrow(tableprint)+1,] <- addvector
		
	} else {
		inputline = nrow(tableprint)+1
		tableprint[inputline,1] <- chosen
		tableprint[inputline,2] <- as.character(printfix)
		tableprint[inputline,3] <- as.character(printrnd)
		tableprint[inputline,4] <- numchains
		tableprint[inputline,5] <- "N/A"
		tableprint[inputline,6] <- "N/A"
		tableprint[inputline,7] <- "N/A"
		tableprint[inputline,8] <- "N/A"
		tableprint[inputline,9] <- "N/A"
		tableprint[inputline,10] <- "N/A"
		tableprint[inputline,11] <- "N/A"
		tableprint[inputline,12] <- "N/A"
		tableprint[inputline,13] <- "N/A"
		tableprint[inputline,14] <- "N/A"
		tableprint[inputline,15] <- "N/A"
		tableprint[inputline,16] <- "N/A"
		
	}
	
	
	
	
	write.table(tableprint, file=tablename, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=TRUE)
	
	
}		#for each file

cat("\n\n\nDONE!\n\n\n\n")
