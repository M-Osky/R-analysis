# R analysis
Some short R scripts I wrote to analyse my data.
This is just a sample of some of the work I have done with R.

Basically you have:
one-way ANOVA and t-test for Normal data
Kruskall-Wallis and Wilcoxon test for non-Normal
And multivariant analysis: PCA and LDA
They include a few lines about the input file needed and so, so should be quite straigh forward to use them.
Usually columns should be sorted as this: descriptors (not analised), grouping variables, numerical variables.

  ind.tag   extrainfo   pop   sex   var1    var2    var3    var4    var5    var6
  Lem042    LemuriaNA   LM    fem   0004    0008    0015    0016    0023    0042  
  Atl016    Atlantida   AT    fem   0003    0014    0015    0092    0062    0036
  Atl028    Atlantida   AT    mal   0013    0021    0034    0055    0089    0144

I'm still working (when I have time) with the RMark scripts, recently found out I could improve some of the models (they weren't 100% correct), as soon as I update them, they will be here. Feedback will be strongly appreciated.
Cheers!
