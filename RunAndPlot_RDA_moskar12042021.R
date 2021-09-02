# R script 
# for for running RDA analysis using on smaller datasets (i.e. tens of varibles, not thousands) vol.4
# Will run basic or partial RDA analysis, extract and plot main results, 
# run permutation test for model, axis and marginal term significance, 
# and perform forward or backward explanatory variable selection for the full model specified in the analysis.

# INPUT FILE MUST BE:
# individual data -> grouping variables -> response variables -> explanatory variables -> constraining variables (no missing values allowed):
#	ind  group1  group2   ... groupn  res1  res2   ...  resn   exp1  exp2  ...  expn  cons1 cons2 ... consn  
#	1	  	04	    M       ...   1     1.5   1.6    ...  2.3    3.2   2.7   ...  4.2   3.5   2.8   ... 1.7
# 2     07      P       ...   2     2.8   5.4    ...  4.2    5.7   5.8   ...  3.5   4.7   2.1   ... 2.9
#	3		  02	    E       ...   2     8.2   8.1    ...  8.2    0.8   1.3   ...  8.4   2.8   0.7   ... 4.3
# ...etc

#needed packages:
# install.packages("vegan")
# install.packages("RColorBrewer")
# install.packages("ggplot2")
# install.packages("formattable")
# install.packages("gridExtra")

library(vegan)
library(RColorBrewer)
library(ggplot2)
library(formattable)
library(gridExtra)








rm(list = ls())		#Remove all objects
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n") # I like clean spaces when trying the same analysis with slight changes or differents datasets


#set working directory
setwd("D:/Dropbox/MOSKY/CURRO/PODARCIS/Measurements/nat_pop_Iva/pheno-env/jun21_RDAnew/forward_select/together")




#set parameters:
input_file = "PCindivMAL_8ecolvar_allMEMforwardGeoEcoREV.csv"     #set the input file name (comma separated)


#import input file
mydata <- read.table(input_file, header = TRUE, sep=",")
str(mydata)






nametag="forwardall_geo-ecol"                       #somedescriptor of the analysis a part from filename and full/partial
stuff = 4                                #set the number of individual or non relevant columns at the left of the file
groupingv = 1                                #set the number of grouping variables to test
responsev = 5                               #set the number of response variables
predictv = 5                                #set the number of predictor (explanatory) variables

constrainv = 0                               #if you do NOT want to run partial analysis set number of variables to be used as constraints as 0
constrainv = 3                               #if you want to run partial analysis set number of variables from dataset to be used as constraints (> 0)

#set the type or types of plots to print
scalebiplot = 3                              #scale predictor and response scores for triplot
scalebiplot =c(1, 2)                              #scale predictor and response scores for triplot

selectvarforw = FALSE                             #perform forward explanatory variable selection (options: TRUE or FALSE)
selectvarback = FALSE                             #perform backward explanatory variable selection (options: TRUE or FALSE)




############## NO NEED TO SET ANYTHING UNDER THIS LINE ###############

###### extract info about variables ######
options(scipen=999)                 #supress scientific notation

id = stuff

firstgroup = id+1
lastgroup = id+groupingv
colgroup = c(firstgroup:lastgroup)
groupnames = colnames(mydata[colgroup])

firstresponse = id+groupingv+1
lastresponse = id+groupingv+responsev
responses = mydata[,firstresponse:lastresponse]
colresp = c(firstresponse:lastresponse)
respnames = colnames(mydata[colresp])

firstpredict = id+groupingv+responsev+1
lastpredict = id+groupingv+responsev+predictv
predictors = mydata[,firstpredict:lastpredict]
colpred = c(firstpredict:lastpredict)
prednames = colnames(mydata[colpred])


if (constrainv > 0) {
	
	firstconstrain = id+groupingv+responsev+predictv+1
	lastconstrain = id+groupingv+responsev+predictv+constrainv
	constraints = mydata[,firstconstrain:lastconstrain]
	colconstrain = c(firstconstrain:lastconstrain)
	constranames = colnames(mydata[colconstrain])

	explainv = mydata[,firstpredict:lastconstrain]
	cat("\n\nPartial RDA ", nametag, " will be performed in", input_file, "\nPredictor variables:", prednames, "\nResponse variables:", names(responses), "\nConstrains (its/their effect will be removed):", constranames, "\nGrouping variable(s):", groupnames, "\n\n", sep=" ")
	filename = gsub('^(.*)\\..*$', '\\1', input_file)
	outputname = paste(filename, "_partlRDA", sep="")
	typerda="Partial"

} else {
	cat("\n\nFull RDA ", nametag, " will be performed in", input_file, "\nPredictor variables:", names(predictors), "\nResponse variables:", names(responses), "\nGrouping variable(s):", groupnames, "\n\n", sep=" ")
	filename = gsub('^(.*)\\..*$', '\\1', input_file)
	#prepare outputs
	outputname = paste(filename, "_basicRDA", sep="")
	typerda="Full"
}




############## CHECK AND RUN!! ###############




####### run RDA analysis ########
if (constrainv > 0) {
  
  form <- formula(paste("responses~", paste(prednames, collapse = "+"), "+ Condition(", paste(constranames, collapse = "+"),")"))
  myrda <- rda(form, data = explainv, scale=TRUE)
  
} else if (constrainv == 0) {
  
  form <- formula(paste("responses~", paste(prednames, collapse = "+")))
  myrda <- rda(form, data=predictors, scale=TRUE)
}

summary(myrda)


####### extract results #######

#prepare the output file with the results of linear regression
filename_res = paste(outputname, "_", nametag,"_out.txt", sep="")

cat("Redundancy analysis (RDA) on:", input_file, "------------------------------------------------------------------------------------", file=filename_res, sep="\n", append=TRUE)
cat("Response variables: ", respnames, "\n", file=filename_res, sep="\n", append=TRUE)
cat("Predictor variables: ", prednames, "\n", file=filename_res, sep="\n", append=TRUE)

if (constrainv > 0) {
  
  cat("Constraining variables: ", constranames, "\n", file=filename_res, sep="\n", append=TRUE)
  cat("------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

} else if (constrainv == 0) {
  
  cat("No constraining variables specified for the analysis.", file=filename_res, sep="\n", append=TRUE)
  cat("------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
}

#extract main results
capture.output(myrda, file=filename_res, append=TRUE)
cat("------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

#str(myrda)

numelem = length(myrda$CA$eig)+length(myrda$CCA$eig)
if (numelem > 10) { numelem = 10}


#extract scores
filename_scores = paste(outputname, "_", nametag, "_scores.txt", sep="")

rdascores_sites <- scores(myrda, display="sites", choices=c(1:numelem))
rdascores_species <- scores(myrda, display="species", choices=c(1:numelem))
rdascores_bp <- scores(myrda, display="bp", choices=c(1:numelem))
options(warn=-1) 
cat(filename, " scores\n\n\nsites\n", sep="", file=filename_scores)
write.table(rdascores_sites, filename_scores, sep="\t", append=TRUE) 
cat("\n\n\nspecies\n", sep="", file=filename_scores, append=TRUE)
write.table(rdascores_species, filename_scores, sep="\t", append=TRUE) 
cat("\n\n\nbp\n", sep="", file=filename_scores, append=TRUE)
write.table(rdascores_bp, filename_scores, sep="\t", append=TRUE) 
options(warn=0) 

#extract regression coefficients
mylrcoeff <- coef(myrda)
cat("\n\n\nLinear Regression Coefficients\n", sep="", file=filename_scores, append=TRUE)
options(warn=-1) 
write.table(mylrcoeff, filename_scores, sep="\t", append=TRUE)
options(warn=0) 


#variance explained
myeigenvalues <- summary(eigenvals(myrda, model = "constrained"))

myxlab <- paste0("RDA1 (", percent(myeigenvalues[2,1]), ")")
myylab <- paste0("RDA2 (", percent(myeigenvalues[2,2]), ")")

capture.output(myeigenvalues, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

 
#extract R2 values
myR2 <- RsquareAdj(myrda)$r.squared
myR2adj <- RsquareAdj(myrda)$adj.r.squared

cat("R2:", myR2, "\n", file=filename_res, sep=" ", append=TRUE)
cat("adjusted R2:", myR2adj, "\n", file=filename_res, sep=" ", append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

myrdaR2 <- data.frame(score=c("R2", "R2adj"), value=c(myR2, myR2adj))

title_rplot <- paste0(outputname, "_", nametag, "_R2.png")
png(filename=title_rplot, units="in", width=5, height=5, res=300)
gp <- ggplot(data=myrdaR2, aes(x=score, y=value)) +
  geom_bar(stat = "identity", color="black", fill=c("forestgreen","darkorange")) +
  scale_x_discrete(limits=myrdaR2$score, labels = c("R2","R2adj")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits=c(0,1)) +
  theme_bw() +
  labs(title=expression(paste(R^2," and adjusted ",R^2," scores"))) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(size=12, face="bold"),
        plot.title = element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(gp)

dev.off()


myrda


#permutation tests for RDA results
myaov <- anova.cca(myrda, step=1000)
myaovax <- anova.cca(myrda, by='axis', step=1000)
myaovmg <- anova.cca(myrda, by='margin', step=1000)

cat("ANOVA permutation test for RDA model:", "\n", file=filename_res, sep="\n", append=TRUE)
capture.output(myaov, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
cat("ANOVA significance for each constrained axis:", "\n", file=filename_res, sep="\n", append=TRUE)
capture.output(myaovax, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
cat("ANOVA significance for marginal effect of predictors:", "\n", file=filename_res, sep="\n", append=TRUE)
capture.output(myaovmg, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)








png(paste(outputname, "_", nametag, "AOV_model.png", sep=""))
a1 <- tableGrob(myaov)
grid.arrange(a1)
dev.off()

png(paste(outputname, "_", nametag, "AOV_axis.png", sep=""))
a2 <- tableGrob(myaovax)
grid.arrange(a2)
dev.off()

png(paste(outputname, "_", nametag, "AOV_marginal_effect.png", sep=""))
a3 <- tableGrob(myaovmg)
grid.arrange(a3)
dev.off()


#goodness of fit
mygoodnessoffit <- goodness(myrda)

cat("Goodness of fit:", "\n", file=filename_res, sep="\n", append=TRUE)
capture.output(mygoodnessoffit, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)


#variance inflation factors
myvif <- vif.cca(myrda)

cat("Variance inflation factors (VIF):", "\n", file=filename_res, sep="\n", append=TRUE)
capture.output(myvif, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)



#forward explanatory variable selection
if (selectvarforw == TRUE) {
  cat("\n\nForward selection")
  if (constrainv > 0) {
  
    cat(" with constrains.\n\n")
    mod0 <- rda(responses~1, data = explainv)
  
  } else if (constrainv == 0) {
    cat(" with no constrains.\n\n")
    mod0 <- rda(responses~1, data = predictors)
  }
  
  cat("\nForward ordistep selection of explanatory variables:", "\n", file=filename_res, sep="\n", append=TRUE)
  capture.output(step.forward <- ordistep(mod0, scope=formula(myrda), direction="forward", permutations=how(nperm=1000)), file=filename_res, append=TRUE)
  myR2adj_ordi <- (RsquareAdj(step.forward))$adj.r.squared
  
  cat("\n", "\n", file=filename_res, sep=" ", append=TRUE)
  cat("Forward ordistep selection adjusted R2:", myR2adj_ordi, file=filename_res, sep=" ", append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

  cat("\nForward ordiR2step selection of explanatory variables:", "\n", file=filename_res, sep="\n", append=TRUE)
  
  capture.output(step2.forward <- ordiR2step(mod0, scope=formula(myrda), direction="forward", R2scope=TRUE, permutations=how(nperm=1000)), file=filename_res, append=TRUE)
  myR2adj_ordi <- (RsquareAdj(step2.forward))$adj.r.squared
  cat("\n", "\n", file=filename_res, sep="\n", append=TRUE)
  cat("Forward ordiR2step ANOVA of the selected model\n", file=filename_res, sep=" ", append=TRUE)
  capture.output(step2.forward$anova, file=filename_res, append=TRUE)
  cat("\n\n\nForward ordiR2step selected model adjusted R2:", myR2adj_ordi, file=filename_res, sep=" ", append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
  
}






#backward explanatory variable selection
if (selectvarback == TRUE) {
  
  myordi <- ordistep(myrda, permutations=how(nperm=1000))
  myR2adj_ordi<- (RsquareAdj(myordi))$adj.r.squared
  
  cat("Backward ordistep selection of explanatory variables:", "\n", file=filename_res, sep="\n", append=TRUE)
  capture.output(myordi, file=filename_res, append=TRUE)
  cat("\n", "\n", file=filename_res, sep=" ", append=TRUE)
  cat("Backward selection adjusted R2:", myR2adj_ordi, file=filename_res, sep=" ", append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
}



#plot
if (length(myeigenvalues)>3) {
  
  
  ####### plot eigenvalues #######
  title_eig <- paste(outputname, "_", nametag, "_RDA_eigenvalues.png", sep="")
  png(filename=title_eig, units="in", width=7, height=5, res=300)
  
  screeplot(myrda, main="RDA eigenvalues")

  dev.off()


  ####### plot triplot #######
  
  for (numscale in scalebiplot) {
  
  numscale
  
  file_scattplot <- paste0(outputname, "_", nametag, "_RDA_triplot-scale", numscale, ".png")
  png(filename=file_scattplot, units="in", width=7, height=7, res=300)
  title_plot=paste(typerda, " RDA ", nametag, " biplot (scaling = ", numscale, ")", sep="")
  plot(myrda, scaling=numscale, type="n", xlab=myxlab, ylab=myylab, font.lab=2, main=title_plot)
  points(myrda, scaling=numscale, display="sites", pch=20, cex=1, col="gray55")                 #add individuals scores
  points(myrda, scaling=numscale, display="species", pch=20, cex=1.5, col="red")                #add phenotypic measures 
  text(myrda, scaling=numscale, display="species", col="red", cex=0.8, font=2, pos=3)           #add phenotypic measures
  #arrows(0, 0, rdascores_species[,1], rdascores_species[,2], lty=1, length=0.05, col="red")       #try adding arrows for species (not working at this point)
  text(myrda, scaling=numscale, display="bp", col="black", cex=1, font=2)                       #add predictor biplot scores

  dev.off()

  }

  ####### plot grouping scores #######
  for ( i in groupnames) {
  
    #extract group levels
	printfact = as.factor(mydata[[i]])
    tags <- levels(printfact)
    numgroups <- length(tags)
   
    #set colours
    if (numgroups<=4)	{
      mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3")
      myshapes=c(22,21,23,24)
    } else if (numgroups>4 & numgroups<=11) {	
      mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3", "black", "red1", "goldenrod1", "purple3", "blue", "magenta", "cyan")
      myshapes=c(21:25,21:25,21)
    } else if (numgroups==12) {	
      mypalette <- brewer.pal(12,"Paired")
      myshapes=c(21:25,21:25,21:22)
    } else if (numgroups==14) {	
      mypalette <- colorRampPalette(brewer.pal(name="Paired", n = 12))(14)
      myshapes=c(21:25,21:25,21:24)
    } else if (numgroups==20) {	
      mypalette <- colorRampPalette(brewer.pal(name="Paired", n = 12))(20)
      myshapes=c(21:25,21:25,21:25,21:25)
    }else if (numgroups==26) {	
      mypalette <- colorRampPalette(brewer.pal(name="Paired", n = 12))(26)
      myshapes=c(21:25,21:25,21:25,21:25,21:25,21)
    }

    #plot
    title_scattplot <- paste0(outputname, "_", nametag, "_", i, "_RDA", ".png")
    png(filename=title_scattplot, units="in", width=7, height=5, res=300)
  
    plot(myrda, type = "n", xlab=myxlab, ylab=myylab, font.lab=2)
    with(mydata, points(myrda, display = "sites", col = mypalette[printfact], pch = 21, bg=mypalette[printfact]))
    with(mydata, legend("topleft", legend = tags, bty = "n", col = mypalette, pch = 21, pt.bg = mypalette, text.font=2, cex=0.55))
  
    dev.off()
  
  }
  
} else if (length(myeigenvalues)<=3) {
  
  print("Can't plot, only one RDA eigenvector obtained.")
}



cat("\n\n\nDONE!!\n\n\n")

