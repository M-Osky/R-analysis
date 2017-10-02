#Partial Least Squares Regression First trial

getwd()
#Change it (set) to the directory including the input file
setwd("D:/Dropbox/WORKING NOW/Ecol/REGRESSION PLS/R")
#Check that is all right
getwd()

library ("pls")

#to not print more decimal possitions than needed
options(digits=4)
#we'll need to allow R to write a lot more lines than by deffault
options(max.print=99999999)

#input file

allys12 <- read.table("input/allys12.txt", header = TRUE)
summary(allys12)

#split in diffeent matrix with groups of variables

Y <- as.matrix(allys12[,2])
X <- as.matrix(allys12[,3:4])
C <- as.matrix(allys12[,5:8]) 
A <- as.matrix(allys12[,9:28]) 
B <- as.matrix(allys12[,29:35]) 

#View(Y)
#View(X)
#View(C)
#View(A)
X


#mydata <-data.frame(y=Y, x1=a, x2=b, x3=c)

#pls.options
#pls.options(plsralg = "oscorespls")

#Fit a PLSR model 10 components, 
#Include leave-one-out (LOO) cross-validatio test

allys12cv.all <- plsr(Y ~ X+A+B+C, ncomp = 4, method = "oscorespls", scale=FALSE, validation="LOO")
allys12cv.X <- plsr(Y ~ X, ncomp = 2, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.A <- plsr(Y ~ A, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.B <- plsr(Y ~ B, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.C <- plsr(Y ~ C, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.AX <- plsr(Y ~ A+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.AB <- plsr(Y ~ A+B, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.AC <- plsr(Y ~ A+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.BX <- plsr(Y ~ B+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.BC <- plsr(Y ~ B+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.CX <- plsr(Y ~ C+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.ABC <- plsr(Y ~ A+B+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.ABX <- plsr(Y ~ A+B+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.ACX <- plsr(Y ~ A+X+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allys12cv.BCX <- plsr(Y ~ X+B+C, ncomp = 2, method = "oscorespls", scale=TRUE, validation="LOO")

summary(allys12cv.all)

#N from different sex summed
sumy12 <- read.table("input/sumy12.txt", header = TRUE)
summary(sumy12)

sY <- as.matrix(sumy12[,2])
sX <- as.matrix(sumy12[,3])
sC <- as.matrix(sumy12[,4:7]) 
sA <- as.matrix(sumy12[,8:27]) 
sB <- as.matrix(sumy12[,28:34]) 

sA

sumy12cv.all <- plsr(sY ~ sX+sA+sB+sC, ncomp =5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.X <- plsr(sY ~ sX, ncomp = 1, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.A <- plsr(sY ~ sA, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.B <- plsr(sY ~ sB, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.C <- plsr(sY ~ sC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.AX <- plsr(sY ~ sA+sX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.AB <- plsr(sY ~ sA+sB, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.AC <- plsr(sY ~ sA+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.BX <- plsr(sY ~ sB+sX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.BC <- plsr(sY ~ sB+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.CX <- plsr(sY ~ sC+sX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.ABC <- plsr(sY ~ sA+sB+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.ABX <- plsr(sY ~ sA+sB+sX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.ACX <- plsr(sY ~ sA+sX+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumy12cv.BCX <- plsr(sY ~ sX+sB+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")

#Only males
maley12 <- read.table("input/maley12.txt", header = TRUE)
summary(maley12)

mY <- as.matrix(maley12[,2])
mX <- as.matrix(maley12[,3])
mC <- as.matrix(maley12[,4:7]) 
mA <- as.matrix(maley12[,8:27]) 
mB <- as.matrix(maley12[,28:34]) 

maley12cv.all <- plsr(mY ~ mX+mA+mB+mC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.X <- plsr(mY ~ mX, ncomp = 1, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.A <- plsr(mY ~ mA, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.B <- plsr(mY ~ mB, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.C <- plsr(mY ~ mC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.AX <- plsr(mY ~ mA+mX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.AB <- plsr(mY ~ mA+mB, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.AC <- plsr(mY ~ mA+mC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.BX <- plsr(mY ~ mB+mX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.BC <- plsr(mY ~ mB+mC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.CX <- plsr(mY ~ mC+mX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.ABC <- plsr(mY ~ mA+mB+mC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.ABX <- plsr(mY ~ mA+mB+mX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.ACX <- plsr(mY ~ mA+mX+mC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maley12cv.BCX <- plsr(mY ~ mX+mB+mC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")


#Only females
femy12 <- read.table("input/femy12.txt", header = TRUE)
summary(femy12)

fY <- as.matrix(femy12[,2])
fX <- as.matrix(femy12[,3])
fC <- as.matrix(femy12[,4:7]) 
fA <- as.matrix(femy12[,8:27]) 
fB <- as.matrix(femy12[,28:34]) 



femy12cv.all <- plsr(fY ~ fX+fA+fB+fC, ncomp = 2, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.X <- plsr(fY ~ fX, ncomp = 1, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.A <- plsr(fY ~ fA, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.B <- plsr(fY ~ fB, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.C <- plsr(fY ~ fC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.AX <- plsr(fY ~ fA+fX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.AB <- plsr(fY ~ fA+fB, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.AC <- plsr(fY ~ fA+fC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.BX <- plsr(fY ~ fB+fX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.BC <- plsr(fY ~ fB+fC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.CX <- plsr(fY ~ fC+fX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.ABC <- plsr(fY ~ fA+fB+fC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.ABX <- plsr(fY ~ fA+fB+fX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.ACX <- plsr(fY ~ fA+fX+fC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femy12cv.BCX <- plsr(fY ~ fX+fB+fC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")

#Eliminados componentes que no aportan varianza (explicada)


#log <- file("all12cvOUT.txt", open="wt")
#sink(log, type="output")


summary(allys12cv.all)
summary(allys12cv.X)
summary(allys12cv.A) 
summary(allys12cv.B) 
summary(allys12cv.C) 
summary(allys12cv.AX)
summary(allys12cv.AB)
summary(allys12cv.AC)
summary(allys12cv.BX)
summary(allys12cv.BC)
summary(allys12cv.CX)
summary(allys12cv.ABX)
summary(allys12cv.ABC)
summary(allys12cv.ACX)
summary(allys12cv.BCX)

summary(sumy12cv.all)
summary(sumy12cv.X)
summary(sumy12cv.A) 
summary(sumy12cv.B) 
summary(sumy12cv.C) 
summary(sumy12cv.AX)
summary(sumy12cv.AB)
summary(sumy12cv.AC)
summary(sumy12cv.BX)
summary(sumy12cv.BC)
summary(sumy12cv.CX)
summary(sumy12cv.ABC)
summary(sumy12cv.ABX)
summary(sumy12cv.ACX)
summary(sumy12cv.BCX)

summary(maley12cv.all)
summary(maley12cv.X)
summary(maley12cv.A) 
summary(maley12cv.B) 
summary(maley12cv.C) 
summary(maley12cv.AX)
summary(maley12cv.AB)
summary(maley12cv.AC)
summary(maley12cv.BX)
summary(maley12cv.BC)
summary(maley12cv.CX)
summary(maley12cv.ABC)
summary(maley12cv.ABX)
summary(maley12cv.ACX)
summary(maley12cv.BCX)

summary(femy12cv.all)
summary(femy12cv.X)
summary(femy12cv.A) 
summary(femy12cv.B) 
summary(femy12cv.C) 
summary(femy12cv.AX)
summary(femy12cv.AB)
summary(femy12cv.AC)
summary(femy12cv.BX)
summary(femy12cv.BC)
summary(femy12cv.CX)
summary(femy12cv.ABC)
summary(femy12cv.ABX)
summary(femy12cv.ACX)
summary(femy12cv.BCX)

#sink(type="message")
#close(log)







plot(RMSEP(allys12cv.all), legendpos = "topright")
plot(RMSEP(allys12cv.X), legendpos = "topright")
plot(RMSEP(allys12cv.A), legendpos = "topright") 
plot(RMSEP(allys12cv.B), legendpos = "topright") 
plot(RMSEP(allys12cv.C), legendpos = "topright") 
plot(RMSEP(allys12cv.AX), legendpos = "topright")
plot(RMSEP(allys12cv.AB), legendpos = "topright")
plot(RMSEP(allys12cv.AC), legendpos = "topright")
plot(RMSEP(allys12cv.BX), legendpos = "topright")
plot(RMSEP(allys12cv.BC), legendpos = "topright")
plot(RMSEP(allys12cv.CX), legendpos = "topright")
plot(RMSEP(allys12cv.ABC), legendpos = "topright")
plot(RMSEP(allys12cv.ABX), legendpos = "topright")
plot(RMSEP(allys12cv.ACX), legendpos = "topright")
plot(RMSEP(allys12cv.BCX), legendpos = "topright")





plot(RMSEP(sumy12cv.all), legendpos = "topright")
plot(RMSEP(sumy12cv.X), legendpos = "topright")
plot(RMSEP(sumy12cv.A), legendpos = "topright") 
plot(RMSEP(sumy12cv.B), legendpos = "topright") 
plot(RMSEP(sumy12cv.C), legendpos = "topright") 
plot(RMSEP(sumy12cv.AX), legendpos = "topright")
plot(RMSEP(sumy12cv.AB), legendpos = "topright")
plot(RMSEP(sumy12cv.AC), legendpos = "topright")
plot(RMSEP(sumy12cv.BX), legendpos = "topright")
plot(RMSEP(sumy12cv.BC), legendpos = "topright")
plot(RMSEP(sumy12cv.CX), legendpos = "topright")
plot(RMSEP(sumy12cv.ABC), legendpos = "topright")
plot(RMSEP(sumy12cv.ABX), legendpos = "topright")
plot(RMSEP(sumy12cv.ACX), legendpos = "topright")
plot(RMSEP(sumy12cv.BCX), legendpos = "topright")





plot(RMSEP(maley12cv.all), legendpos = "topright")
plot(RMSEP(maley12cv.X), legendpos = "topright")
plot(RMSEP(maley12cv.A), legendpos = "topright") 
plot(RMSEP(maley12cv.B), legendpos = "topright") 
plot(RMSEP(maley12cv.C), legendpos = "topright") 
plot(RMSEP(maley12cv.AX), legendpos = "topright")
plot(RMSEP(maley12cv.AB), legendpos = "topright")
plot(RMSEP(maley12cv.AC), legendpos = "topright")
plot(RMSEP(maley12cv.BX), legendpos = "topright")
plot(RMSEP(maley12cv.BC), legendpos = "topright")
plot(RMSEP(maley12cv.CX), legendpos = "topright")
plot(RMSEP(maley12cv.ABC), legendpos = "topright")
plot(RMSEP(maley12cv.ABX), legendpos = "topright")
plot(RMSEP(maley12cv.ACX), legendpos = "topright")
plot(RMSEP(maley12cv.BCX), legendpos = "topright")


plot(RMSEP(femy12cv.all), legendpos = "topright")
plot(RMSEP(femy12cv.X), legendpos = "topright")
plot(RMSEP(femy12cv.A), legendpos = "topright") 
plot(RMSEP(femy12cv.B), legendpos = "topright") 
plot(RMSEP(femy12cv.C), legendpos = "topright") 
plot(RMSEP(femy12cv.AX), legendpos = "topright")
plot(RMSEP(femy12cv.AB), legendpos = "topright")
plot(RMSEP(femy12cv.AC), legendpos = "topright")
plot(RMSEP(femy12cv.BX), legendpos = "topright")
plot(RMSEP(femy12cv.BC), legendpos = "topright")
plot(RMSEP(femy12cv.CX), legendpos = "topright")
plot(RMSEP(femy12cv.ABC), legendpos = "topright")
plot(RMSEP(femy12cv.ABX), legendpos = "topright")
plot(RMSEP(femy12cv.ACX), legendpos = "topright")
plot(RMSEP(femy12cv.BCX), legendpos = "topright")




#######################################################



#Sin datos de 2012

allys <- read.table("input/allys.txt", header = TRUE)
summary(allys)

#split in diffeent matrix with groups of variables

Y <- as.matrix(allys[,2])
X <- as.matrix(allys[,3:4])
C <- as.matrix(allys[,5:8]) 
A <- as.matrix(allys[,9:28]) 
B <- as.matrix(allys[,29:35]) 

#View(Y)
#View(X)
#View(C)
#View(A)
B


#mydata <-data.frame(y=Y, x1=a, x2=b, x3=c)

#pls.options
#pls.options(plsralg = "oscorespls")

#Fit a PLSR model 10 components, 
#Include leave-one-out (LOO) cross-validatio test

allyscv.all <- plsr(Y ~ X+A+B+C, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.X <- plsr(Y ~ X, ncomp = 1, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.A <- plsr(Y ~ A, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.B <- plsr(Y ~ B, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.C <- plsr(Y ~ C, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.AX <- plsr(Y ~ A+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.AB <- plsr(Y ~ A+B, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.AC <- plsr(Y ~ A+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.BX <- plsr(Y ~ B+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.BC <- plsr(Y ~ B+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.CX <- plsr(Y ~ C+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.ABC <- plsr(Y ~ A+B+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.ABX <- plsr(Y ~ A+B+X, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.ACX <- plsr(Y ~ A+X+C, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
allyscv.BCX <- plsr(Y ~ X+B+C, ncomp = 2, method = "oscorespls", scale=TRUE, validation="LOO")


#N from different sex summed
sumy <- read.table("input/sumy.txt", header = TRUE)
summary(sumy)

sY <- as.matrix(sumy[,2])
sX <- as.matrix(sumy[,3])
sC <- as.matrix(sumy[,4:7]) 
sA <- as.matrix(sumy[,8:27]) 
sB <- as.matrix(sumy[,28:34]) 

sB

sumycv.all <- plsr(sY ~ sX+sA+sB+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.X <- plsr(sY ~ sX, ncomp = 1, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.A <- plsr(sY ~ sA, ncomp = 8, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.B <- plsr(sY ~ sB, ncomp = 6, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.C <- plsr(sY ~ sC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.AX <- plsr(sY ~ sA+sX, ncomp = 9, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.AB <- plsr(sY ~ sA+sB, ncomp = 13, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.AC <- plsr(sY ~ sA+sC, ncomp = 10, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.BX <- plsr(sY ~ sB+sX, ncomp = 7, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.BC <- plsr(sY ~ sB+sC, ncomp = 10, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.CX <- plsr(sY ~ sC+sX, ncomp = 4, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.ABC <- plsr(sY ~ sA+sB+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.ABX <- plsr(sY ~ sA+sB+sX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.ACX <- plsr(sY ~ sA+sX+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
sumycv.BCX <- plsr(sY ~ sX+sB+sC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")

#Only males
maley <- read.table("input/maley.txt", header = TRUE)
summary(maley)

mY <- as.matrix(maley[,2])
mX <- as.matrix(maley[,3])
mC <- as.matrix(maley[,4:7]) 
mA <- as.matrix(maley[,8:27]) 
mB <- as.matrix(maley[,28:34]) 

mB

maleycv.all <- plsr(mY ~ mX+mA+mB+mC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.X <- plsr(mY ~ mX, ncomp = 1, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.A <- plsr(mY ~ mA, ncomp = 8, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.B <- plsr(mY ~ mB, ncomp = 6, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.C <- plsr(mY ~ mC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.AX <- plsr(mY ~ mA+mX, ncomp = 9, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.AB <- plsr(mY ~ mA+mB, ncomp = 14, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.AC <- plsr(mY ~ mA+mC, ncomp = 10, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.BX <- plsr(mY ~ mB+mX, ncomp = 7, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.BC <- plsr(mY ~ mB+mC, ncomp = 10, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.CX <- plsr(mY ~ mC+mX, ncomp = 4, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.ABC <- plsr(mY ~ mA+mB+mC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.ABX <- plsr(mY ~ mA+mB+mX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.ACX <- plsr(mY ~ mA+mX+mC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
maleycv.BCX <- plsr(mY ~ mX+mB+mC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")



#Only fems
femy <- read.table("input/femy.txt", header = TRUE)
summary(femy)

fY <- as.matrix(femy[,2])
fX <- as.matrix(femy[,3])
fC <- as.matrix(femy[,4:7]) 
fA <- as.matrix(femy[,8:27]) 
fB <- as.matrix(femy[,28:34]) 

fB

femycv.all <- plsr(fY ~ fX+fA+fB+fC, ncomp = 2, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.X <- plsr(fY ~ fX, ncomp = 1, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.A <- plsr(fY ~ fA, ncomp = 8, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.B <- plsr(fY ~ fB, ncomp = 6, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.C <- plsr(fY ~ fC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.AX <- plsr(fY ~ fA+fX, ncomp = 9, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.AB <- plsr(fY ~ fA+fB, ncomp = 14, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.AC <- plsr(fY ~ fA+fC, ncomp = 10, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.BX <- plsr(fY ~ fB+fX, ncomp = 7, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.BC <- plsr(fY ~ fB+fC, ncomp = 10, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.CX <- plsr(fY ~ fC+fX, ncomp = 4, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.ABC <- plsr(fY ~ fA+fB+fC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.ABX <- plsr(fY ~ fA+fB+fX, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.ACX <- plsr(fY ~ fA+fX+fC, ncomp = 5, method = "oscorespls", scale=TRUE, validation="LOO")
femycv.BCX <- plsr(fY ~ fX+fB+fC, ncomp = 3, method = "oscorespls", scale=TRUE, validation="LOO")



#log <- file("allcvOUT.txt", open="wt")
#sink(log, type="output")



summary(allyscv.all)
summary(allyscv.X)
summary(allyscv.A) 
summary(allyscv.B) 
summary(allyscv.C) 
summary(allyscv.AX)
summary(allyscv.AB)
summary(allyscv.AC)
summary(allyscv.BX)
summary(allyscv.BC)
summary(allyscv.CX)
summary(allyscv.ABC)
summary(allyscv.ABX)
summary(allyscv.ACX)
summary(allyscv.BCX)

summary(sumycv.all)
summary(sumycv.X)
summary(sumycv.A) 
summary(sumycv.B) 
summary(sumycv.C) 
summary(sumycv.AX)
summary(sumycv.AB)
summary(sumycv.AC)
summary(sumycv.BX)
summary(sumycv.BC)
summary(sumycv.CX)
summary(sumycv.ABC)
summary(sumycv.ABX)
summary(sumycv.ACX)
summary(sumycv.BCX)

summary(maleycv.all)
summary(maleycv.X)
summary(maleycv.A) 
summary(maleycv.B) 
summary(maleycv.C) 
summary(maleycv.AX)
summary(maleycv.AB)
summary(maleycv.AC)
summary(maleycv.BX)
summary(maleycv.BC)
summary(maleycv.CX)
summary(maleycv.ABC)
summary(maleycv.ABX)
summary(maleycv.ACX)
summary(maleycv.BCX)


#sink(type="message")
#close(log)


log <- file("R2all.txt", open="wt")
sink(log, type="output")

R2(allys12cv.CX)
R2(allys12cv.C)
R2(allys12cv.X)
R2(allys12cv.A)
R2(allys12cv.B)
R2(allys12cv.BC)
R2(allys12cv.ABC)
R2(allys12cv.BCX)
R2(allys12cv.all)

R2(sumy12cv.CX)
R2(sumy12cv.C)
R2(sumy12cv.X)
R2(sumy12cv.A)
R2(sumy12cv.B)
R2(sumy12cv.BC)
R2(sumy12cv.ABC)
R2(sumy12cv.BCX)
R2(sumy12cv.all)

R2(maley12cv.CX)
R2(maley12cv.C)
R2(maley12cv.X)
R2(maley12cv.A)
R2(maley12cv.B)
R2(maley12cv.BC)
R2(maley12cv.ABC)
R2(maley12cv.BCX)
R2(maley12cv.all)


R2(allyscv.ABC)
R2(allyscv.BCX)
R2(allyscv.all)


R2(sumycv.ABC)
R2(sumycv.BCX)
R2(sumycv.all)


R2(maleycv.BCX)
R2(maleycv.ABC)
R2(maleycv.all)


R2(femycv.ABC)
R2(femycv.BCX)
R2(femycv.all)

R2(femy12cv.CX)
R2(femy12cv.C)
R2(femy12cv.X)
R2(femy12cv.A)
R2(femy12cv.B)
R2(femy12cv.BC)
R2(femy12cv.ABC)
R2(femy12cv.BCX)
R2(femy12cv.all)

sink(type="output")
close(log)



plot(RMSEP(allyscv.all), legendpos = "topright")
plot(RMSEP(allyscv.X), legendpos = "topright")
plot(RMSEP(allyscv.A), legendpos = "topright") 
plot(RMSEP(allyscv.B), legendpos = "topright") 
plot(RMSEP(allyscv.C), legendpos = "topright") 
plot(RMSEP(allyscv.AX), legendpos = "topright")
plot(RMSEP(allyscv.AB), legendpos = "topright")
plot(RMSEP(allyscv.AC), legendpos = "topright")
plot(RMSEP(allyscv.BX), legendpos = "topright")
plot(RMSEP(allyscv.BC), legendpos = "topright")
plot(RMSEP(allyscv.CX), legendpos = "topright")
plot(RMSEP(allyscv.ABC), legendpos = "topright")
plot(RMSEP(allyscv.ABX), legendpos = "topright")
plot(RMSEP(allyscv.ACX), legendpos = "topright")
plot(RMSEP(allyscv.BCX), legendpos = "topright")





plot(RMSEP(sumycv.all), legendpos = "topright")
plot(RMSEP(sumycv.X), legendpos = "topright")
plot(RMSEP(sumycv.A), legendpos = "topright") 
plot(RMSEP(sumycv.B), legendpos = "topright") 
plot(RMSEP(sumycv.C), legendpos = "topright") 
plot(RMSEP(sumycv.AX), legendpos = "topright")
plot(RMSEP(sumycv.AB), legendpos = "topright")
plot(RMSEP(sumycv.AC), legendpos = "topright")
plot(RMSEP(sumycv.BX), legendpos = "topright")
plot(RMSEP(sumycv.BC), legendpos = "topright")
plot(RMSEP(sumycv.CX), legendpos = "topright")
plot(RMSEP(sumycv.ABC), legendpos = "topright")
plot(RMSEP(sumycv.ABX), legendpos = "topright")
plot(RMSEP(sumycv.ACX), legendpos = "topright")
plot(RMSEP(sumycv.BCX), legendpos = "topright")




plot(RMSEP(maleycv.all), legendpos = "topright")
plot(RMSEP(maleycv.X), legendpos = "topright")
plot(RMSEP(maleycv.A), legendpos = "topright") 
plot(RMSEP(maleycv.B), legendpos = "topright") 
plot(RMSEP(maleycv.C), legendpos = "topright") 
plot(RMSEP(maleycv.AX), legendpos = "topright")
plot(RMSEP(maleycv.AB), legendpos = "topright")
plot(RMSEP(maleycv.AC), legendpos = "topright")
plot(RMSEP(maleycv.BX), legendpos = "topright")
plot(RMSEP(maleycv.BC), legendpos = "topright")
plot(RMSEP(maleycv.CX), legendpos = "topright")
plot(RMSEP(maleycv.ABC), legendpos = "topright")
plot(RMSEP(maleycv.ABX), legendpos = "topright")
plot(RMSEP(maleycv.ACX), legendpos = "topright")
plot(RMSEP(maleycv.BCX), legendpos = "topright")

#Refited models allys12 X; C; CX; malesy12 C; allys X; C and CX with the ncomp
# (Modified above in the original lines)

#To see the power of each model we need to calculate Q^2 = cross-validated R^2




#VIP
### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
### `pls' package.

## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.

VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}

#Tying to figure it out


vipallys12.cx  <- VIP(allys12cv.CX)
vipallys12.c   <- VIP(allys12cv.C)
vipallys12.x   <- VIP(allys12cv.X)
vipallys12.a   <- VIP(allys12cv.A)
vipallys12.b   <- VIP(allys12cv.B)
vipallys12.bc  <- VIP(allys12cv.BC)
vipallys12.abc  <- VIP(allys12cv.ABC)
vipallys12.bcx <- VIP(allys12cv.BCX)
vipallys12.all <- VIP(allys12cv.all)

vipsumy12.cx  <- VIP(sumy12cv.CX)
vipsumy12.c   <- VIP(sumy12cv.C)
vipsumy12.x   <- VIP(sumy12cv.X)
vipsumy12.a   <- VIP(sumy12cv.A)
vipsumy12.b   <- VIP(sumy12cv.B)
vipsumy12.bc  <-  VIP(sumy12cv.BC)
vipsumy12.abc  <-  VIP(sumy12cv.ABC)
vipsumy12.bcx <-  VIP(sumy12cv.BCX)
vipsumy12.all <- VIP(sumy12cv.all)

vipmaley12.cx  <- VIP(maley12cv.CX)
vipmaley12.c   <- VIP(maley12cv.C)
vipmaley12.x   <- VIP(maley12cv.X)
vipmaley12.a   <- VIP(maley12cv.A)
vipmaley12.b   <- VIP(maley12cv.B)
vipmaley12.bc  <- VIP(maley12cv.BC)
vipmaley12.abc  <- VIP(maley12cv.ABC)
vipmaley12.bcx <- VIP(maley12cv.BCX)
vipmaley12.all <- VIP(maley12cv.all)


vipallys.cx  <- VIP(allyscv.CX)
vipallys.c   <- VIP(allyscv.C)
vipallys.x   <- VIP(allyscv.X)
vipallys.a   <- VIP(allyscv.A)
vipallys.b   <- VIP(allyscv.B)
vipallys.bc  <- VIP(allyscv.BC)
vipallys.abc  <- VIP(allyscv.ABC)
vipallys.bcx <- VIP(allyscv.BCX)
vipallys.all <- VIP(allyscv.all)

vipsumy.cx  <- VIP(sumycv.CX)
vipsumy.c   <- VIP(sumycv.C)
vipsumy.x   <- VIP(sumycv.X)
vipsumy.a   <- VIP(sumycv.A)
vipsumy.b   <- VIP(sumycv.B)
vipsumy.bc  <- VIP(sumycv.BC)
vipsumy.abc <- VIP(sumycv.ABC)
vipsumy.bcx <- VIP(sumycv.BCX)
vipsumy.all <- VIP(sumycv.all)

vipmaley.cx  <- VIP(maleycv.CX)
vipmaley.c   <- VIP(maleycv.C)
vipmaley.x   <- VIP(maleycv.X)
vipmaley.a   <- VIP(maleycv.A)
vipmaley.b   <- VIP(maleycv.B)
vipmaley.bc  <- VIP(maleycv.BC)
vipmaley.bcx  <- VIP(maleycv.BCX)
vipmaley.abc <- VIP(maleycv.ABC)
vipmaley.all <- VIP(maleycv.all)

vipfemy.cx  <- VIP(femycv.CX)
vipfemy.c   <- VIP(femycv.C)
vipfemy.x   <- VIP(femycv.X)
vipfemy.a   <- VIP(femycv.A)
vipfemy.b   <- VIP(femycv.B)
vipfemy.bc  <- VIP(femycv.BC)
vipfemy.abc  <- VIP(femycv.ABC)
vipfemy.bcx <- VIP(femycv.BCX)
vipfemy.all <- VIP(femycv.all)

vipfemy12.cx  <- VIP(femy12cv.CX)
vipfemy12.c   <- VIP(femy12cv.C)
vipfemy12.x   <- VIP(femy12cv.X)
vipfemy12.a   <- VIP(femy12cv.A)
vipfemy12.b   <- VIP(femy12cv.B)
vipfemy12.bc  <- VIP(femy12cv.BC)
vipfemy12.abc  <- VIP(femy12cv.ABC)
vipfemy12.bcx <- VIP(femy12cv.BCX)
vipfemy12.all <- VIP(femy12cv.all)


#log <- file("VIPs.txt", open="wt")
#sink(log, type="output")
#vipallys.cx



write.table (vipallys12.cx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipallys12.c, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipallys12.x, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipallys12.a, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipallys12.b, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipallys12.bc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipallys12.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipallys12.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipallys12.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)

write.table (vipsumy12.cx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipsumy12.c, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy12.x, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy12.a, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy12.b, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy12.bc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipsumy12.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipsumy12.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipsumy12.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 

write.table (vipmaley12.cx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipmaley12.c, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley12.x, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley12.a, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley12.b, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley12.bc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipmaley12.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipmaley12.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipmaley12.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)


write.table (vipallys.cx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipallys.c, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipallys.x, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipallys.a, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipallys.b, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipallys.bc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipallys.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipallys.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipallys.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)

write.table (vipsumy.cx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipsumy.c, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy.x, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy.a, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy.b, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)   
write.table (vipsumy.bc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipsumy.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipsumy.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipsumy.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)

write.table (vipmaley.cx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipmaley.c, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley.x, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley.a, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley.b, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipmaley.bc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipmaley.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipmaley.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipmaley.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)

write.table (vipfemy.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipfemy.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipfemy.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)

write.table (vipfemy12.cx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipfemy12.c, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipfemy12.x, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipfemy12.a, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipfemy12.b, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)  
write.table (vipfemy12.bc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE) 
write.table (vipfemy12.bcx, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipfemy12.abc, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)
write.table (vipfemy12.all, sep= "; ", eol="\r\n", dec=".", na=" ", row.names=TRUE, col.names=NA, file="VIPs.txt", append=TRUE, quote=FALSE)



#sink(type="message")
#close(log)

#vipallys12.cx
#vipallys12.c
#summary(allys12cv.all)
