#Example of a quick R script that can be used for running RMark.
#It is quite specific for my data so may not be useful, but I keep it here as a reference
#POPAN CHORRILLO 11: Two plots (j and i); two files each (naz and ALL) + All entire area


getwd()
#Change it (set) to the directory including the input file
setwd("D:/Dropbox/WORKING NOW/Ecol/Datos MarkRecap/MARK/R")
#Check that is all right
getwd()


#Mark library

library(RMark)


e11naz.data=convert.inp("files/e11naz.inp", group.df= data.frame(sex=c("Female", "Male")), covariates="Span")
str(e11naz.data)

str(e11naz.data)
Span <- e11naz.data$Span
span <- factor(e11naz.data$Span, labels = c("fresh", "good", "old", "worn"))
e11naz.data$span=span
#Check again
str(e11naz.data)
e11naz.data

e11naz.pd=process.data(e11naz.data, model="POPAN", groups="sex", age.unit=0.35, begin.time=110705, time.intervals=c(3,1,1,2,3,3,3))
str(e11naz.pd)


#Now let's create the design data library to set up all the parameters of the chosen model
e11naz.dd=make.design.data(e11naz.pd)
str(e11naz.dd)



#which parameters it includes?
names(e11naz.dd)

e11naz.dd$Phi
str(e11naz.dd$Phi)
e11naz.dd$p
e11naz.dd$pent



#we need to adjust age so the maximum is 3 (if > 3); 0 if < 0

print(e11naz.dd$p)

#add the levels of "age" including the round numbers, just in case they don't exist
#levp <- levels(e11naz.dd$p$age)
#levels(e11naz.dd$p$age) <- c("0", "1", "2", "3" , levp )
#levels(e11naz.dd$p$age)

#Now round up the values
print(e11naz.dd$p)
i=substr(e11naz.dd$p$age,1,4)<=0.54
e11naz.dd$p$Age[i]=0
i=((substr(e11naz.dd$p$age,1,4)>=0.55)&(substr(e11naz.dd$p$age,1,4)<=1.71))
e11naz.dd$p$Age[i]=1
i=((substr(e11naz.dd$p$age,1,4)>=1.72)&(substr(e11naz.dd$p$age,1,4)<=2.88))
e11naz.dd$p$Age[i]=2
i=substr(e11naz.dd$p$age,1,4)>=2.89
e11naz.dd$p$Age[i]=3
print(e11naz.dd$p)



print(e11naz.dd$Phi)

i=substr(e11naz.dd$Phi$age,1,4)<=0.54
e11naz.dd$Phi$Age[i]=0
i=((substr(e11naz.dd$Phi$age,1,4)>=0.55)&(substr(e11naz.dd$Phi$age,1,4)<=1.71))
e11naz.dd$Phi$Age[i]=1
i=((substr(e11naz.dd$Phi$age,1,4)>=1.72)&(substr(e11naz.dd$Phi$age,1,4)<=2.88))
e11naz.dd$Phi$Age[i]=2
i=substr(e11naz.dd$Phi$age,1,4)>=2.89
e11naz.dd$Phi$Age[i]=3
print(e11naz.dd$Phi)



print(e11naz.dd$pent)

i=substr(e11naz.dd$pent$age,1,4)<=0.54
e11naz.dd$pent$Age[i]=0
i=((substr(e11naz.dd$pent$age,1,4)>=0.55)&(substr(e11naz.dd$pent$age,1,4)<=1.71))
e11naz.dd$pent$Age[i]=1
i=((substr(e11naz.dd$pent$age,1,4)>=1.72)&(substr(e11naz.dd$pent$age,1,4)<=2.88))
e11naz.dd$pent$Age[i]=2
i=substr(e11naz.dd$pent$age,1,4)>=2.89
e11naz.dd$pent$Age[i]=3
print(e11naz.dd$pent)



#Now let's add the days with bad weather, the not good ones as "wind" will include wind, cloudy, rain and storm, and the really bad ones (we weren't even able to sample) as "storm" will only include the days with an actual storm or really strong dangerous winds
#We'll add them to Phi (survival), p (capture) and PENT so they will affect to the next day after their occurrence, for recapture (p) the same day (or previous).

e11naz.dd$Phi$windPhi=0 #We add the new "event"
e11naz.dd$Phi$windPhi[e11naz.dd$Phi$time==110709  | e11naz.dd$Phi$time==110710]=1 #set the event to affect they day after it happened
e11naz.dd$Phi$stormPhi=0 #We add the new "event"
e11naz.dd$Phi$stormPhi[e11naz.dd$Phi$time==110710  | e11naz.dd$Phi$time==110712  | e11naz.dd$Phi$time==110718]=2 ## If the storm happened the last day it doesn't affect survival in terms of our data
e11naz.dd$Phi



#pent
e11naz.dd$pent$windb=0 #We add the new "event"
e11naz.dd$pent$windb[e11naz.dd$pent$time==110709  | e11naz.dd$pent$time==110710]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
e11naz.dd$pent$stormb=0 #We add the new "event"
e11naz.dd$pent$stormb[e11naz.dd$pent$time==110712  | e11naz.dd$pent$time==110721]=2
e11naz.dd$pent

#p
e11naz.dd$p$windp=0 #We add the new "event"
e11naz.dd$p$windp[e11naz.dd$p$time==110709  | e11naz.dd$p$time==110710]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
e11naz.dd$p$stormp=0 #We add the new "event"
e11naz.dd$p$stormp[e11naz.dd$p$time==110712]=2
e11naz.dd$p

#effort
e11naz.dd$p$effort=1 #We add the new "event"
e11naz.dd$p$effort[e11naz.dd$p$time==110709 | e11naz.dd$p$time==110710 | e11naz.dd$p$time==110715 | e11naz.dd$p$time==110718]=2
e11naz.dd$p

names(e11naz.dd$p)
names(e11naz.dd$Phi)
names(e11naz.dd$pent)




#We can set a function which will include all the factors that we specify
#and then runs all the possible combinations 

#we'll need to allow R to write a lot more lines than by deffault
options(max.print=99999999)

log <- file("loge11naz.txt", open="wt")
sink(log, type="message")

#Run models for each possible combination of the preselected relevant factors affecting the parameters
e11naz.models=function()
{
  
  p.001=list(formula=~time)
  p.002=list(formula=~time+Age)
  p.003=list(formula=~time*Age)
  p.005=list(formula=~time+windp+stormp)
  p.013=list(formula=~effort+time)
  p.014=list(formula=~effort+time+Age)
  
  Phi.009=list(formula=~windPhi)
  Phi.012=list(formula=~stormPhi)
  
  pent.007=list(formula=~Time)
  pent.008=list(formula=~Time+Age)
  pent.013=list(formula=~Age)
  pent.014=list(formula=~Age+windb)
  pent.015=list(formula=~Age+windb+stormb)
  pent.016=list(formula=~windb)

  N.001=list(formula=~sex)
  
  cml=create.model.list("POPAN")
  results=mark.wrapper(cml,data=e11naz.pd, ddl=e11naz.dd,adjust=FALSE)
  return(results)
}

#Now the function has the info, but we need to retrieve the results

e11naz.results=e11naz.models()

sink(type="message")
close(log)


#Con Dropbox y modo normal fallan modelos aleatoriamente, sin DropBox y en modo Administrador No
#OK Dropbox seems to interfere with the analysis don't run analysis inside a Dropbox folder unless it is switch off
#Models

e11naz.results


export.MARK(e11naz.pd, "e11best", e11naz.results, replace = FALSE, ind.covariates = "all")

model.average(e11naz.results)

#export model to MARK to check if the models are good (C_hat )
#release.gof(dipper.proc) #First GOF to check if they are good (non-significant)
#For some reason is not working, let's try in MARK

#export.chdata(e11naz.pd, "e11nazMARK", replace="TRUE") #create an input
#Open "New" in MARK and select the created inp file
#export.model(e11naz.results) #export model (or list of models)

#Click on Browse once, then go to Output>Append and select the last temp file

#Go to Test> Program RELEASE GOF
#check GOF (summary test 2 and 3)

#Test>BOOTSTRAP GOF -> RUN
#Simulations>View simulations; we can click on the calculator to get the mean
#Observed C-hat, method 1: Observed deviance of our model / mean deviance of simulated (bootstrap)
#Observed C-hat, method 2: (model deviance / deviance degrees of freedom) / mean simulated (bootstrap) c-hat ## Observed Deviance DF = Select model, click button "view output", look for "DEVIANCE Degrees of Freedom" in the document
#From both methods, choose the correction for the observed c-hat more different from 1

#models.adjusted=adjust.chat(1.346, e11naz.results)
#models.adjusted


#We sould be ready now to look at the parametres estimated from the best models 


##########################################################
##########################################################

#e11all


e11all.data=convert.inp("files/e11all.inp", group.df= data.frame(sex=c("Female", "Male")), covariates="Span")
str(e11all.data)

str(e11all.data)
Span <- e11all.data$Span
span <- factor(e11all.data$Span, labels = c("fresh", "good", "old", "worn"))
e11all.data$span=span
#Check again
str(e11all.data)
e11all.data

e11all.pd=process.data(e11all.data, model="POPAN", groups="sex", age.unit=0.35, begin.time=110705, time.intervals=c(3,1,1,2,3,3,3))
str(e11all.pd)


#Now let's create the design data library to set up all the parameters of the chosen model
e11all.dd=make.design.data(e11all.pd)
str(e11all.dd)



#which parameters it includes?
names(e11all.dd)

e11all.dd$Phi
str(e11all.dd$Phi)
e11all.dd$p
e11all.dd$pent



#we need to adjust age so the maximum is 3 (if > 3); 0 if < 0

print(e11all.dd$p)

#add the levels of "age" including the round numbers, just in case they don't exist
#levp <- levels(e11all.dd$p$age)
#levels(e11all.dd$p$age) <- c("0", "1", "2", "3" , levp )
#levels(e11all.dd$p$age)

#Now round up the values
print(e11all.dd$p)
i=substr(e11all.dd$p$age,1,4)<=0.54
e11all.dd$p$Age[i]=0
i=((substr(e11all.dd$p$age,1,4)>=0.55)&(substr(e11all.dd$p$age,1,4)<=1.71))
e11all.dd$p$Age[i]=1
i=((substr(e11all.dd$p$age,1,4)>=1.72)&(substr(e11all.dd$p$age,1,4)<=2.88))
e11all.dd$p$Age[i]=2
i=substr(e11all.dd$p$age,1,4)>=2.89
e11all.dd$p$Age[i]=3
print(e11all.dd$p)



print(e11all.dd$Phi)

i=substr(e11all.dd$Phi$age,1,4)<=0.54
e11all.dd$Phi$Age[i]=0
i=((substr(e11all.dd$Phi$age,1,4)>=0.55)&(substr(e11all.dd$Phi$age,1,4)<=1.71))
e11all.dd$Phi$Age[i]=1
i=((substr(e11all.dd$Phi$age,1,4)>=1.72)&(substr(e11all.dd$Phi$age,1,4)<=2.88))
e11all.dd$Phi$Age[i]=2
i=substr(e11all.dd$Phi$age,1,4)>=2.89
e11all.dd$Phi$Age[i]=3
print(e11all.dd$Phi)



print(e11all.dd$pent)

i=substr(e11all.dd$pent$age,1,4)<=0.54
e11all.dd$pent$Age[i]=0
i=((substr(e11all.dd$pent$age,1,4)>=0.55)&(substr(e11all.dd$pent$age,1,4)<=1.71))
e11all.dd$pent$Age[i]=1
i=((substr(e11all.dd$pent$age,1,4)>=1.72)&(substr(e11all.dd$pent$age,1,4)<=2.88))
e11all.dd$pent$Age[i]=2
i=substr(e11all.dd$pent$age,1,4)>=2.89
e11all.dd$pent$Age[i]=3
print(e11all.dd$pent)



#Now let's add the days with bad weather, the not good ones as "wind" will include wind, cloudy, rain and storm, and the really bad ones (we weren't even able to sample) as "storm" will only include the days with an actual storm or really strong dangerous winds
#We'll add them to Phi (survival), p (capture) and PENT so they will affect to the next day after their occurrence, for recapture (p) the same day (or previous).

e11all.dd$Phi$windPhi=0 #We add the new "event"
e11all.dd$Phi$windPhi[e11all.dd$Phi$time==110709  | e11all.dd$Phi$time==110710]=1 #set the event to affect they day after it happened
e11all.dd$Phi$stormPhi=0 #We add the new "event"
e11all.dd$Phi$stormPhi[e11all.dd$Phi$time==110710  | e11all.dd$Phi$time==110712  | e11all.dd$Phi$time==110718]=2 ## If the storm happened the last day it doesn't affect survival in terms of our data
e11all.dd$Phi



#pent
e11all.dd$pent$windb=0 #We add the new "event"
e11all.dd$pent$windb[e11all.dd$pent$time==110709  | e11all.dd$pent$time==110710]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
e11all.dd$pent$stormb=0 #We add the new "event"
e11all.dd$pent$stormb[e11all.dd$pent$time==110712  | e11all.dd$pent$time==110721]=2
e11all.dd$pent

#p
e11all.dd$p$windp=0 #We add the new "event"
e11all.dd$p$windp[e11all.dd$p$time==110709  | e11all.dd$p$time==110710]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
e11all.dd$p$stormp=0 #We add the new "event"
e11all.dd$p$stormp[e11all.dd$p$time==110712]=2
e11all.dd$p

#effort
e11all.dd$p$effort=1 #We add the new "event"
e11all.dd$p$effort[e11all.dd$p$time==110709 | e11all.dd$p$time==110710 | e11all.dd$p$time==110715 | e11all.dd$p$time==110718]=2
e11all.dd$p

names(e11all.dd$p)
names(e11all.dd$Phi)
names(e11all.dd$pent)




#We can set a function which will include all the factors that we specify
#and then runs all the possible combinations 

#we'll need to allow R to write a lot more lines than by deffault
options(max.print=99999999)

log <- file("loge11all.txt", open="wt")
sink(log, type="message")

#Run models for each possible combination of the preselected relevant factors affecting the parameters
e11all.models=function()
{
  
  p.002=list(formula=~sex+time)
  p.003=list(formula=~sex+time+Span)
  p.007=list(formula=~sex+Age+time)
  
  
  Phi.001=list(formula=~time)
  Phi.002=list(formula=~time+windPhi)
  Phi.003=list(formula=~time+windPhi+stormPhi)
  Phi.006=list(formula=~time*Span)
  
  pent.003=list(formula=~time*Age)
  pent.010=list(formula=~Time+stormb+windb)
  pent.014=list(formula=~windb)
  
  N.001=list(formula=~sex)
  
  
  cml=create.model.list("POPAN")
  results=mark.wrapper(cml,data=e11all.pd, dd=e11all.dd,adjust=FALSE)
  return(results)
}

#Now the function has the info, but we need to retrieve the results

e11all.results=e11all.models()

sink(type="message")
close(log)


#Con Dropbox y modo normal fallan modelos aleatoriamente, sin DropBox y en modo Administrador No
#OK Dropbox seems to interfere with the analysis don't run analysis inside a Dropbox folder unless it is switch off
#Models

e11all.results

export.MARK(e11all.pd, "e11allbest", e11all.results, replace = FALSE, ind.covariates = "all")

#export model to MARK to check if the models are good (C_hat )
#release.gof(dipper.proc) #First GOF to check if they are good (non-significant)
#For some reason is not working, let's try in MARK

export.chdata(e11all.pd, "e11allMARK", replace="TRUE") #create an input
#Open "New" in MARK and select the created inp file
export.model(e11all.results) #export model (or list of models)

#Click on Browse once, then go to Output>Append and select the last temp file

#Go to Test> Program RELEASE GOF
#check GOF (summary test 2 and 3)

#Test>BOOTSTRAP GOF -> RUN
#Simulations>View simulations; we can click on the calculator to get the mean
#Observed C-hat, method 1: Observed deviance of our model / mean deviance of simulated (bootstrap)
#Observed C-hat, method 2: (model deviance / deviance degrees of freedom) / mean simulated (bootstrap) c-hat ## Observed Deviance DF = Select model, click button "view output", look for "DEVIANCE Degrees of Freedom" in the document
#From both methods, choose the correction for the observed c-hat more different from 1

#models.adjusted=adjust.chat(1.346, e11all.results)
#models.adjusted


#We sould be ready now to look at the parametres estimated from the best models 



####################################################
####################################################

#f11all

f11all.data=convert.inp("files/f11all.inp", group.df= data.frame(sex=c("Female", "Male")), covariates="Span")
str(f11all.data)

str(f11all.data)
Span <- f11all.data$Span
span <- factor(f11all.data$Span, labels = c("fresh", "good", "old", "worn"))
f11all.data$span=span
#Check again
str(f11all.data)
f11all.data

f11all.pd=process.data(f11all.data, model="POPAN", groups="sex", age.unit=0.35, begin.time=110709, time.intervals=c(3,3,3,3))
str(f11all.pd)


#Now let's create the design data library to set up all the parameters of the chosen model
f11all.dd=make.design.data(f11all.pd)
str(f11all.dd)



#which parameters it includes?
names(f11all.dd)

f11all.dd$Phi
str(f11all.dd$Phi)
f11all.dd$p
f11all.dd$pent



#we need to adjust age so the maximum is 3 (if > 3); 0 if < 0

print(f11all.dd$p)

#add the levels of "age" including the round numbers, just in case they don't exist
#levp <- levels(f11all.dd$p$age)
#levels(f11all.dd$p$age) <- c("0", "1", "2", "3" , levp )
#levels(f11all.dd$p$age)

#Now round up the values
print(f11all.dd$p)
i=substr(f11all.dd$p$age,1,4)<=0.54
f11all.dd$p$Age[i]=0
i=((substr(f11all.dd$p$age,1,4)>=0.55)&(substr(f11all.dd$p$age,1,4)<=1.71))
f11all.dd$p$Age[i]=1
i=((substr(f11all.dd$p$age,1,4)>=1.72)&(substr(f11all.dd$p$age,1,4)<=2.88))
f11all.dd$p$Age[i]=2
i=substr(f11all.dd$p$age,1,4)>=2.89
f11all.dd$p$Age[i]=3
print(f11all.dd$p)



print(f11all.dd$Phi)

i=substr(f11all.dd$Phi$age,1,4)<=0.54
f11all.dd$Phi$Age[i]=0
i=((substr(f11all.dd$Phi$age,1,4)>=0.55)&(substr(f11all.dd$Phi$age,1,4)<=1.71))
f11all.dd$Phi$Age[i]=1
i=((substr(f11all.dd$Phi$age,1,4)>=1.72)&(substr(f11all.dd$Phi$age,1,4)<=2.88))
f11all.dd$Phi$Age[i]=2
i=substr(f11all.dd$Phi$age,1,4)>=2.89
f11all.dd$Phi$Age[i]=3
print(f11all.dd$Phi)



print(f11all.dd$pent)

i=substr(f11all.dd$pent$age,1,4)<=0.54
f11all.dd$pent$Age[i]=0
i=((substr(f11all.dd$pent$age,1,4)>=0.55)&(substr(f11all.dd$pent$age,1,4)<=1.71))
f11all.dd$pent$Age[i]=1
i=((substr(f11all.dd$pent$age,1,4)>=1.72)&(substr(f11all.dd$pent$age,1,4)<=2.88))
f11all.dd$pent$Age[i]=2
i=substr(f11all.dd$pent$age,1,4)>=2.89
f11all.dd$pent$Age[i]=3
print(f11all.dd$pent)



#Now let's add the days with bad weather, the not good ones as "wind" will include wind, cloudy, rain and storm, and the really bad ones (we weren't even able to sample) as "storm" will only include the days with an actual storm or really strong dangerous winds
#We'll add them to Phi (survival), p (capture) and PENT so they will affect to the next day after their occurrence, for recapture (p) the same day (or previous).

f11all.dd$Phi$windPhi=0 #We add the new "event"
f11all.dd$Phi$windPhi[f11all.dd$Phi$time==110709]=1 #set the event to affect they day after it happened
f11all.dd$Phi$stormPhi=0 #We add the new "event"
f11all.dd$Phi$stormPhi[f11all.dd$Phi$time==110709 | f11all.dd$Phi$time==110712 | f11all.dd$Phi$time==110718]=2 ## If the storm happened the last day it doesn't affect survival in terms of our data
f11all.dd$Phi



#pent
f11all.dd$pent$windb=0 #We add the new "event"
f11all.dd$pent$windb[f11all.dd$pent$time==110712]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
f11all.dd$pent$stormb=0 #We add the new "event"
f11all.dd$pent$stormb[f11all.dd$pent$time==110712 | f11all.dd$pent$time==110721]=2
f11all.dd$pent

#p
f11all.dd$p$windp=0 #We add the new "event"
f11all.dd$p$windp[f11all.dd$p$time==110709]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
f11all.dd$p$stormp=0 #We add the new "event"
f11all.dd$p$stormp[f11all.dd$p$time==110712]=2
f11all.dd$p

#effort
f11all.dd$p$effort=1 #We add the new "event"
f11all.dd$p$effort[f11all.dd$p$time==110709 | f11all.dd$p$time==110715 | f11all.dd$p$time==110718]=2
f11all.dd$p

names(f11all.dd$p)
names(f11all.dd$Phi)
names(f11all.dd$pent)



#We can set a function which will include all the factors that we specify
#and then runs all the possible combinations 

#we'll need to allow R to write a lot more lines than by deffault
options(max.print=99999999)

log <- file("logf11all.txt", open="wt")
sink(log, type="message")

#Run models for each possible combination of the preselected relevant factors affecting the parameters
f11all.models=function()
{
  
  p.001=list(formula=~time)
  p.002=list(formula=~time+Age)
  p.003=list(formula=~time*Age)
  p.004=list(formula=~time+effort)
  p.007=list(formula=~Age)
  p.008=list(formula=~Age+effort)
  p.011=list(formula=~Age+stormp)
  
  Phi.005=list(formula=~time+Span)
  Phi.006=list(formula=~time*Span)
  Phi.007=list(formula=~time+Age+Span)
  
  pent.001=list(formula=~1)
  pent.002=list(formula=~sex)
  pent.003=list(formula=~sex+windb)
  pent.004=list(formula=~sex+windb+stormb)
  pent.005=list(formula=~sex+Span)
  pent.006=list(formula=~sex+Span+windb)
  pent.007=list(formula=~Span)
  pent.009=list(formula=~Span+Time+sex)
  
  N.001=list(formula=~sex)
  
  
  cml=create.model.list("POPAN")
  results=mark.wrapper(cml,data=f11all.pd, ddl=f11all.dd,adjust=FALSE)
  return(results)
}

#Now the function has the info, but we need to retrieve the results

f11all.results=f11all.models()

sink(type="message")
close(log)


#Con Dropbox y modo normal fallan modelos aleatoriamente, sin DropBox y en modo Administrador No
#OK Dropbox seems to interfere with the analysis don't run analysis inside a Dropbox folder unless it is switch off
#Models

f11all.results

export.MARK(f11all.pd, "f11allbest", f11all.results, replace = FALSE, ind.covariates = "all")

#export model to MARK to check if the models are good (C_hat )
#release.gof(dipper.proc) #First GOF to check if they are good (non-significant)
#For some reason is not working, let's try in MARK

export.chdata(f11all.pd, "f11allMARK", replace="TRUE") #create an input
#Open "New" in MARK and select the created inp file
export.model(f11all.results) #export model (or list of models)

#Click on Browse once, then go to Output>Append and select the last temp file

#Go to Test> Program RELEASE GOF
#check GOF (summary test 2 and 3)

#Test>BOOTSTRAP GOF -> RUN
#Simulations>View simulations; we can click on the calculator to get the mean
#Observed C-hat, method 1: Observed deviance of our model / mean deviance of simulated (bootstrap)
#Observed C-hat, method 2: (model deviance / deviance degrees of freedom) / mean simulated (bootstrap) c-hat ## Observed Deviance DF = Select model, click button "view output", look for "DEVIANCE Degrees of Freedom" in the document
#From both methods, choose the correction for the observed c-hat more different from 1

#models.adjusted=adjust.chat(1.346, f11all.results)
#models.adjusted


#We sould be ready now to look at the parametres estimated from the best models 




#####################################################
#####################################################

#f11naz

f11naz.data=convert.inp("files/f11naz.inp", group.df= data.frame(sex=c("Female", "Male")), covariates="Span")
str(f11naz.data)

str(f11naz.data)
Span <- f11naz.data$Span
span <- factor(f11naz.data$Span, labels = c("fresh", "good", "old", "worn"))
f11naz.data$span=span
#Check again
str(f11naz.data)
f11naz.data

f11naz.pd=process.data(f11naz.data, model="POPAN", groups="sex", age.unit=0.35, begin.time=110709, time.intervals=c(3,3,3,3))
str(f11naz.pd)


#Now let's create the design data library to set up all the parameters of the chosen model
f11naz.dd=make.design.data(f11naz.pd)
str(f11naz.dd)



#which parameters it includes?
names(f11naz.dd)

f11naz.dd$Phi
str(f11naz.dd$Phi)
f11naz.dd$p
f11naz.dd$pent



#we need to adjust age so the maximum is 3 (if > 3); 0 if < 0

print(f11naz.dd$p)

#add the levels of "age" including the round numbers, just in case they don't exist
#levp <- levels(f11naz.dd$p$age)
#levels(f11naz.dd$p$age) <- c("0", "1", "2", "3" , levp )
#levels(f11naz.dd$p$age)

#Now round up the values
print(f11naz.dd$p)
i=substr(f11naz.dd$p$age,1,4)<=0.54
f11naz.dd$p$Age[i]=0
i=((substr(f11naz.dd$p$age,1,4)>=0.55)&(substr(f11naz.dd$p$age,1,4)<=1.71))
f11naz.dd$p$Age[i]=1
i=((substr(f11naz.dd$p$age,1,4)>=1.72)&(substr(f11naz.dd$p$age,1,4)<=2.88))
f11naz.dd$p$Age[i]=2
i=substr(f11naz.dd$p$age,1,4)>=2.89
f11naz.dd$p$Age[i]=3
print(f11naz.dd$p)



print(f11naz.dd$Phi)

i=substr(f11naz.dd$Phi$age,1,4)<=0.54
f11naz.dd$Phi$Age[i]=0
i=((substr(f11naz.dd$Phi$age,1,4)>=0.55)&(substr(f11naz.dd$Phi$age,1,4)<=1.71))
f11naz.dd$Phi$Age[i]=1
i=((substr(f11naz.dd$Phi$age,1,4)>=1.72)&(substr(f11naz.dd$Phi$age,1,4)<=2.88))
f11naz.dd$Phi$Age[i]=2
i=substr(f11naz.dd$Phi$age,1,4)>=2.89
f11naz.dd$Phi$Age[i]=3
print(f11naz.dd$Phi)



print(f11naz.dd$pent)

i=substr(f11naz.dd$pent$age,1,4)<=0.54
f11naz.dd$pent$Age[i]=0
i=((substr(f11naz.dd$pent$age,1,4)>=0.55)&(substr(f11naz.dd$pent$age,1,4)<=1.71))
f11naz.dd$pent$Age[i]=1
i=((substr(f11naz.dd$pent$age,1,4)>=1.72)&(substr(f11naz.dd$pent$age,1,4)<=2.88))
f11naz.dd$pent$Age[i]=2
i=substr(f11naz.dd$pent$age,1,4)>=2.89
f11naz.dd$pent$Age[i]=3
print(f11naz.dd$pent)



#Now let's add the days with bad weather, the not good ones as "wind" will include wind, cloudy, rain and storm, and the really bad ones (we weren't even able to sample) as "storm" will only include the days with an actual storm or really strong dangerous winds
#We'll add them to Phi (survival), p (capture) and PENT so they will affect to the next day after their occurrence, for recapture (p) the same day (or previous).


f11naz.dd$Phi$windPhi=0 #We add the new "event"
f11naz.dd$Phi$windPhi[f11naz.dd$Phi$time==110709]=1 #set the event to affect they day after it happened
f11naz.dd$Phi$stormPhi=0 #We add the new "event"
f11naz.dd$Phi$stormPhi[f11naz.dd$Phi$time==110709 | f11naz.dd$Phi$time==110712 | f11naz.dd$Phi$time==110718]=2 ## If the storm happened the last day it doesn't affect survival in terms of our data
f11naz.dd$Phi



#pent
f11naz.dd$pent$windb=0 #We add the new "event"
f11naz.dd$pent$windb[f11naz.dd$pent$time==110712]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
f11naz.dd$pent$stormb=0 #We add the new "event"
f11naz.dd$pent$stormb[f11naz.dd$pent$time==110712 | f11naz.dd$pent$time==110721]=2
f11naz.dd$pent

#p
f11naz.dd$p$windp=0 #We add the new "event"
f11naz.dd$p$windp[f11naz.dd$p$time==110709]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
f11naz.dd$p$stormp=0 #We add the new "event"
f11naz.dd$p$stormp[f11naz.dd$p$time==110712]=2
f11naz.dd$p

#effort
f11naz.dd$p$effort=1 #We add the new "event"
f11naz.dd$p$effort[f11naz.dd$p$time==110709 | f11naz.dd$p$time==110715 | f11naz.dd$p$time==110718]=2
f11naz.dd$p

names(f11naz.dd$p)
names(f11naz.dd$Phi)
names(f11naz.dd$pent)



#We can set a function which will include all the factors that we specify
#and then runs all the possible combinations 

#we'll need to allow R to write a lot more lines than by deffault
options(max.print=99999999)



log <- file("logf11naz.txt", open="wt")
sink(log, type="message")

#Run models for each possible combination of the preselected relevant factors affecting the parameters
f11naz.models=function()
{
  
  p.001=list(formula=~time+effort)
  p.002=list(formula=~time+effort+stormp)
  p.003=list(formula=~time+stormp)
  p.004=list(formula=~time+stormp+windp)
  p.005=list(formula=~time)
  p.007=list(formula=~effort+stormp)
  p.008=list(formula=~effort+stormp+windp)
  
  Phi.007=list(formula=~Span)
  Phi.008=list(formula=~Span+windPhi)
  
  pent.001=list(formula=~1)
  pent.010=list(formula=~Span)
  
  N.001=list(formula=~sex)
  
  
  cml=create.model.list("POPAN")
  results=mark.wrapper(cml,data=f11naz.pd, ddl=f11naz.dd,adjust=FALSE)
  return(results)
}

#Now the function has the info, but we need to retrieve the results

f11naz.results=f11naz.models()

sink(type="message")
close(log)


#Con Dropbox y modo normal fallan modelos aleatoriamente, sin DropBox y en modo Administrador No
#OK Dropbox seems to interfere with the analysis don't run analysis inside a Dropbox folder unless it is switch off
#Models

f11naz.results

export.MARK(f11naz.pd, "f11nazbest", f11naz.results, replace = FALSE, ind.covariates = "all")

#export model to MARK to check if the models are good (C_hat )
#release.gof(dipper.proc) #First GOF to check if they are good (non-significant)
#For some reason is not working, let's try in MARK

export.chdata(f11naz.pd, "f11nazMARK", replace="TRUE") #create an input
#Open "New" in MARK and select the created inp file
export.model(f11naz.results) #export model (or list of models)

#Click on Browse once, then go to Output>Append and select the last temp file

#Go to Test> Program RELEASE GOF
#check GOF (summary test 2 and 3)

#Test>BOOTSTRAP GOF -> RUN
#Simulations>View simulations; we can click on the calculator to get the mean
#Observed C-hat, method 1: Observed deviance of our model / mean deviance of simulated (bootstrap)
#Observed C-hat, method 2: (model deviance / deviance degrees of freedom) / mean simulated (bootstrap) c-hat ## Observed Deviance DF = Select model, click button "view output", look for "DEVIANCE Degrees of Freedom" in the document
#From both methods, choose the correction for the observed c-hat more different from 1

#models.adjusted=adjust.chat(1.346, f11naz.results)
#models.adjusted


#We sould be ready now to look at the parametres estimated from the best models 




#####################################################
#####################################################

#cho11all

cho11all.data=convert.inp("files/cho11all.inp", group.df= data.frame(sex=c("Female", "Male")), covariates="Span")
str(cho11all.data)

str(cho11all.data)
Span <- cho11all.data$Span
span <- factor(cho11all.data$Span, labels = c("fresh", "good", "old", "worn"))
cho11all.data$span=span
#Check again
str(cho11all.data)
cho11all.pd=process.data(cho11all.data, model="POPAN", groups="sex", age.unit=0.35, begin.time=110705, time.intervals=c(3,1,1,2,3,3,3))
str(cho11all.pd)


#Now let's create the design data library to set up all the parameters of the chosen model
cho11all.dd=make.design.data(cho11all.pd)
str(cho11all.dd)



#which parameters it includes?
names(cho11all.dd)

cho11all.dd$Phi
str(cho11all.dd$Phi)
cho11all.dd$p
cho11all.dd$pent



#we need to adjust age so the maximum is 3 (if > 3); 0 if < 0

print(cho11all.dd$p)

#add the levels of "age" including the round numbers, just in case they don't exist
#levp <- levels(cho11all.dd$p$age)
#levels(cho11all.dd$p$age) <- c("0", "1", "2", "3" , levp )
#levels(cho11all.dd$p$age)

#Now round up the values
print(cho11all.dd$p)
i=substr(cho11all.dd$p$age,1,4)<=0.54
cho11all.dd$p$Age[i]=0
i=((substr(cho11all.dd$p$age,1,4)>=0.55)&(substr(cho11all.dd$p$age,1,4)<=1.71))
cho11all.dd$p$Age[i]=1
i=((substr(cho11all.dd$p$age,1,4)>=1.72)&(substr(cho11all.dd$p$age,1,4)<=2.88))
cho11all.dd$p$Age[i]=2
i=substr(cho11all.dd$p$age,1,4)>=2.89
cho11all.dd$p$Age[i]=3
print(cho11all.dd$p)



print(cho11all.dd$Phi)

i=substr(cho11all.dd$Phi$age,1,4)<=0.54
cho11all.dd$Phi$Age[i]=0
i=((substr(cho11all.dd$Phi$age,1,4)>=0.55)&(substr(cho11all.dd$Phi$age,1,4)<=1.71))
cho11all.dd$Phi$Age[i]=1
i=((substr(cho11all.dd$Phi$age,1,4)>=1.72)&(substr(cho11all.dd$Phi$age,1,4)<=2.88))
cho11all.dd$Phi$Age[i]=2
i=substr(cho11all.dd$Phi$age,1,4)>=2.89
cho11all.dd$Phi$Age[i]=3
print(cho11all.dd$Phi)



print(cho11all.dd$pent)

i=substr(cho11all.dd$pent$age,1,4)<=0.54
cho11all.dd$pent$Age[i]=0
i=((substr(cho11all.dd$pent$age,1,4)>=0.55)&(substr(cho11all.dd$pent$age,1,4)<=1.71))
cho11all.dd$pent$Age[i]=1
i=((substr(cho11all.dd$pent$age,1,4)>=1.72)&(substr(cho11all.dd$pent$age,1,4)<=2.88))
cho11all.dd$pent$Age[i]=2
i=substr(cho11all.dd$pent$age,1,4)>=2.89
cho11all.dd$pent$Age[i]=3
print(cho11all.dd$pent)



#Now let's add the days with bad weather, the not good ones as "wind" will include wind, cloudy, rain and storm, and the really bad ones (we weren't even able to sample) as "storm" will only include the days with an actual storm or really strong dangerous winds
#We'll add them to Phi (survival), p (capture) and PENT so they will affect to the next day after their occurrence, for recapture (p) the same day (or previous).

cho11all.dd$Phi$windPhi=0 #We add the new "event"
cho11all.dd$Phi$windPhi[cho11all.dd$Phi$time==110709  | cho11all.dd$Phi$time==110710]=1 #set the event to affect they day after it happened
cho11all.dd$Phi$stormPhi=0 #We add the new "event"
cho11all.dd$Phi$stormPhi[cho11all.dd$Phi$time==110710  | cho11all.dd$Phi$time==110712  | cho11all.dd$Phi$time==110718]=2 ## If the storm happened the last day it doesn't affect survival in terms of our data
cho11all.dd$Phi



#pent
cho11all.dd$pent$windb=0 #We add the new "event"
cho11all.dd$pent$windb[cho11all.dd$pent$time==110709  | cho11all.dd$pent$time==110710]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
cho11all.dd$pent$stormb=0 #We add the new "event"
cho11all.dd$pent$stormb[cho11all.dd$pent$time==110712  | cho11all.dd$pent$time==110721]=2
cho11all.dd$pent

#p
cho11all.dd$p$windp=0 #We add the new "event"
cho11all.dd$p$windp[cho11all.dd$p$time==110709  | cho11all.dd$p$time==110710]=1 #set the event to affect the day it happened, the first day we don't recapture, so doesn't make sense
cho11all.dd$p$stormp=0 #We add the new "event"
cho11all.dd$p$stormp[cho11all.dd$p$time==110712]=2
cho11all.dd$p

#effort
cho11all.dd$p$effort=1 #We add the new "event"
cho11all.dd$p$effort[cho11all.dd$p$time==110709 | cho11all.dd$p$time==110710 | cho11all.dd$p$time==110715 | cho11all.dd$p$time==110718]=2
cho11all.dd$p

names(cho11all.dd$p)
names(cho11all.dd$Phi)
names(cho11all.dd$pent)



#We can set a function which will include all the factors that we specify
#and then runs all the possible combinations 

#we'll need to allow R to write a lot more lines than by deffault
options(max.print=99999999)

log <- file("logcho11all.txt", open="wt")
sink(log, type="message")

#Run models for each possible combination of the preselected relevant factors affecting the parameters
cho11all.models=function()
{
  
  p.004=list(formula=~sex+Span+time)
  
  Phi.001=list(formula=~time)
  Phi.002=list(formula=~time+Span)
  Phi.004=list(formula=~time+Age)
  Phi.005=list(formula=~time+Age+Span)
  Phi.006=list(formula=~time*Age)
  Phi.012=list(formula=~Time*Age)
  Phi.014=list(formula=~Span+Age)
  Phi.016=list(formula=~Age)
  
  pent.001=list(formula=~time)
  pent.003=list(formula=~time*Age)
  pent.005=list(formula=~time+windb)
  pent.006=list(formula=~time+windb+stormb)
  pent.008=list(formula=~Time+Age)
  pent.009=list(formula=~Time*Age)
  pent.010=list(formula=~Time+Age+windb)
  pent.011=list(formula=~Time+windb)
  pent.012=list(formula=~Time+windb+stormb)
  pent.013=list(formula=~Age)
  pent.014=list(formula=~Age+windb)
  pent.016=list(formula=~windb+stormb)
  
  N.001=list(formula=~sex)
  
  
  cml=create.model.list("POPAN")
  results=mark.wrapper(cml,data=cho11all.pd, ddl=cho11all.dd,adjust=FALSE)
  return(results)
}

#Now the function has the info, but we need to retrieve the results

cho11all.results=cho11all.models()

sink(type="message")
close(log)


#Con Dropbox y modo normal fallan modelos aleatoriamente, sin DropBox y en modo Administrador No
#OK Dropbox seems to interfere with the analysis don't run analysis inside a Dropbox folder unless it is switch off
#Models

cho11all.results

export.MARK(cho11all.pd, "cho11allbest", cho11all.results, replace = FALSE, ind.covariates = "all")

#export model to MARK to check if the models are good (C_hat )
#release.gof(dipper.proc) #First GOF to check if they are good (non-significant)
#For some reason is not working, let's try in MARK

export.chdata(cho11all.pd, "cho11allMARK", replace="TRUE") #create an input
#Open "New" in MARK and select the created inp file
export.model(cho11all.results) #export model (or list of models)

#Click on Browse once, then go to Output>Append and select the last temp file

#Go to Test> Program RELEASE GOF
#check GOF (summary test 2 and 3)

#Test>BOOTSTRAP GOF -> RUN
#Simulations>View simulations; we can click on the calculator to get the mean
#Observed C-hat, method 1: Observed deviance of our model / mean deviance of simulated (bootstrap)
#Observed C-hat, method 2: (model deviance / deviance degrees of freedom) / mean simulated (bootstrap) c-hat ## Observed Deviance DF = Select model, click button "view output", look for "DEVIANCE Degrees of Freedom" in the document
#From both methods, choose the correction for the observed c-hat more different from 1

#models.adjusted=adjust.chat(1.346, cho11all.results)
#models.adjusted


#We sould be ready now to look at the parametres estimated from the best models 
                         
