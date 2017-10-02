#Partial Least Squares Regression TUTORIAL

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

#example data sets, 60 sample with two response variables: Octane number and NIR spectra
  #and 401 reflectance measurements at differents nm 
data("gasoline")
data("yarn")



# First divide the dataset in two, train and test data sets

gasTrain <- gasoline[1:50, ]
gasTest <- gasoline[51:60, ]

gasTest

#typical way of fitting a PLSR model 10 components, 
  #and including leave-one-out (LOO) cross-validated predictions
gas1 <- plsr(octane ~ NIR, ncomp = 10, data = gasTrain, validation = "LOO")

#overview of the fit and validation results
summary(gas1)

#We'll use The root mean squared error of prediction (RMSEP) to validate the data
#The values of CV will give us a result for the LOO. adjCV which is the same but unbiased
#If validation = "LOO", the values of adjCV ad CV should be the same

#But usually is easier to see by plottig it. legendpos argument adds a legend at the indicated position
plot(RMSEP(gas1), legendpos = "topright")
#In this case two COMPONENTS seem to be more tha enough


#Once the number of components has been chosen
#Let's inspect different aspects of the fit by plotting predictions, scores, loadings, etc.

#Prediction Plot: cross-validated predictions versus the measured values with 2 components, 
  #aspect ratio of 1 and trend line
plot(gas1, ncomp = 2, asp = 1, line = TRUE)
#The points follow the line quite nicely 

#Now pairwise plot of the score values for the three first components
  #to look for patterns, groups or outliers in the data
plot(gas1, plottype = "scores", comps = 1:3)
#no clear indication of grouping or outliers. 

#The numbers in parentheses after the component labels are the relative amount of
  #X variance explained by each component.
#explained variances can be extracted explicitly
explvar(gas1)


# Now a loading plot for interpretation purposes, to look for known spectral peaks or profiles:
plot(gas1, "loadings", comps = 1:3, legendpos = "topleft", labels = "numbers", xlab = "nm")
abline(h = 0)
#The >labels = "numbers"< argument makes the plot function try to interpret
#the variable "names" as numbers, and use them as x axis labels.


#A fitted model can be used to predict the response values of new observations.
#The following predicts responses for the ten observations in gasTest, 
#using two components from gasTrain (gas1)
predict(gas1, ncomp = 2, newdata = gasTest)

#In this case we actualy know the true response values, so we can calculate it for gasTest
#RMSEP
RMSEP(gas1, newdata = gasTest)
#For two components with this new samples, we get 0.244, which is quite close to the cross-validated
  #estimate (0.297) <--CV values of the RMSEP obtained with the 50 obs of (gasTrain)
#0.2966 was the CV value for two components.


                      ##MANUAL##



# FORMULAS IN R

#left-hand side (response) + ~ + right-hand side (regressor/s)
  #function:> bodymass ~ age, season, food

#It is also possible to specify transformations of the variables. For instance, log(y) ~ msc(Z)
  #regression of the logarithm of "y" onto "Z" after "Z" has been transformed by multiplicative scatter (or signal) correction (MSC)
#If the transformations contain symbols to be interpreted in the formula +, * or ^
  #should be protected with the I() function
  #this:> y ~ x1 + I(x2 + x3). This specifies two regressors: x1, and the sum of x2 and x3.

#To collect the variables in a dataframe with "data.frame"
  #if v1, v2 and v3 are factors or numeric vectors
  #>mydata <- data.frame(y = v1, a = v2, b = v3)
#will result in a data frame with variables named y, a and b.

#If Z is a matrix, it has to be protected by the 'protect function' I() when using data.frame
  #>mydata <- data.frame(..., Z = I(Z)).
  #Otherwise, it will be split into separate variables for each column and there will be no variable called Z in the data frame.

#We can always add data 'a posteriori' in the data frame
  #>mydata <- data.frame(...)
  #>mydata$Z <- Z
  #this will also respect Z as a matrix

#We can use "cbind" to combine vectors and matrices into matrices on the fly in the formula
  #>cbind(y1, y2) ~ X



# FITTINNG MODELS

#The main function to fit models is plsr, is a wrapper for functio mvr using the right algorithm

#>plsr(formula, ncomp, data)

#The argument "formula" is a formula as described above
#argument "ncomp" is the number of components one wishes to fit (if blank it uses the MAX)
#"data" is the data frame containing the variables to use in the model (if variables in formula aren't already in the environment)
#The function returns a fitted model (an object of class "mvr")
#If the response term of the formula is a matrix, a multi-response model is fit

#To try different alternatives of the model, we can use the function "update"
  #>dens3 <- update(dens1, ncomp = 10)
#updating with new conditions or data

#The predictor variables are always centered, as a part of the fit algorithm.
#Scaling can be requested with the "scale" argument.
#If scale is TRUE, each variable is standardised by dividing it by its standard deviation
#if scale is a numeric vector, each variable is divided by the corresponding number



# CHOOSING THE NUMBER OF COMPONENTS (CROSS-VALIDATION)

#Cross-validation is used to determine the optimal number of components to take into account
#In the plsr function the default value for "validation" (the argument for cross-validation) is "NONE"
  #Validation = "CV" will divide the data into segments, by default 10 segments randomly selected, 
    #but also segments of consecutive objects or interleaved segments are possible through the use of the argument
    #segment.type One can also specify the segments explicitly with the argument segments. See ?mvrCv for details.
  #validation= "LOO" will provide leaveone-out cross-validation.


#the model then will then contain an element comprising the extra information
  #and extra predicted values, as well as MSEP and R2 values. 
  #and the MSEP error using no components at all
#This validation results can be visualised using the plottype = "validation" argument of the standard plotting function.

#the crossvalidation procedure is not optimal. Cross-validation errors are biased downward.
#As long as the only purpose is to select the optimal number of components, this bias may not be very important,
  #but it is not too difficult to avoid it.
  #argument scale that can be used for auto-scaling per segment.
# more elaborate methods such as MSC need explicit handling per segment.
  #For this, function "crossval" is available. It takes an mvr object and performs the cross-validation as it
    #should be done: applying the pre-treatment for each segment

    #> gas2.cv <- crossval(gas2, segments = 10)
    #> plot(MSEP(gas2.cv), legendpos = "topright")
    #> summary(gas2.cv, what = "validation")
    
    #When the scaling does not depend on the division of the data into segments (e.g., log-scaling),
    #functions crossval and mvrCv give the same results; however, crossval is much slower.



# INSPECTING THE FITTED MODELS

#The most visual mode is by plotting them
#wen using the function >plot the default plot type (argument "plottype") is the prediction plot ("predplot")
    #showing predicted versus measured values
  #We can choose other values for plottype in ">plot"

#The number of rows and columns are chosen automatically, can be changed with arguments nRows and nCols

#If we provided a "test set" with the "newdata" argument by deffault it will use the test set.
#If the model was built using cross-validation, the cross-validated predictions are used
#If neither, it will use the predictions for the training set.
#This deffaults can be overriden with the "which" argument

#If we want to know how many components are optimal we should call a validation plot ("validationplot")
#of prediction performance (RMSEP, MSEP, or R2) against the number of components

#regression coefficients can be visualised using plottype = "coef" in >plot or directly through function >coefplot
#either way with this we can plot simultaneously the regression vectors for several different numbers of components at once
  #in despite of nccom=2; ncom= 1:2 (or whatever interval of numbers)

#Scores and loadings can be plotted using functions >scoreplot and >loadingplot
  #One can indicate the number of components with the "comps" argument

#With the function >corrplot; or in >plot, the argument plottype = "correlation"
  #we can see correlations between each variable and the selected components
  #it'll show a scatter plot of two sets of scores with concentric circles, each point = 1 variable
  #The squared distance between the point and the center equals the fraction of the variance 
    #of the variable explained by the ploted components
  #The argument "radii" by deffault shows a circle corresponding with 50% of explained variance and another of 75%


#Check and Summarize

# >print shows: type of regression, form of validation, and the function called
#>summary also shows: the amount of variance explained by the model for training and validation
  #With the argument "what" we ccan choose only either training phase or validation phase


# SELECTING FIT ALGORITHM

#There is three algorithms: 
  #the deffault Kernel algorith (recommended for samples > variables)
  #the classic Orthogonal scores algorithm (NIPALS)
  #the SIMPLS algorith (better for multi-response models)
#We can change the algorithm through the "method" argument method="oscorespls" but we need to do it EVERY TIME

#We can change it for the whole session (or sessions if we change the global environment) with the pls.options
#>pls.option will show the options for PLSR, PCR and MVR that we can change 
#>pls.options(plsralg = "oscorespls") will set the orthogonal scores algorithm for plsr analysis

#We can also call each function separately with kernelpls.fit, oscorespls.fit, and simpls.fit
  #They all take arguments X, Y,ncomp, and stripped. Arguments X, Y, (as matrices) and ncomp (number of components)
  #also if argument stripped=TRUE the calculations are stripped down to the bare minimum required
    #This only returns the X means, Y means, and the regression coefficients. To speed up cross-validation procedures



























