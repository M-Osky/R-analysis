# R ANALYSIS
Some short R scripts I wrote to analyse my data.
This is just a sample of some of the bioestatistics work I do with R either for our own research or for the students.

- The scripts done to be run in desktop R have some lines on the top describing the analysis, options and input file
- The scripts done to be run from command line (unix) have their own help information with the flag --help

# Short description

- HWtest - Analysis per population of loci out of Hardy-Weinberg equilibrium, filter out according to a p-value and percentage.
- LDA - Incomplete Linear Discriminant Analysis script, it works but needs to be adjusted and annotated properly for each case, still could be useful for some people starting
- Normality - test normality in multiple files and variables, plots distributions and tests different data transformation
- parametric tests (ANOVA) - Useful script when doing exploratory analysis of multiple files, we use it to check if the same variables behave similarly in different datasets
- Plot Fst Matrix - Plot per population pair-wise Fst hemimatrices (tested with Stampp package outputs)
- Plot Manhattan Fst - Quick script to plot per-locus Fst pair-wise values from populations (Stacks) output
- Pst_mcmc - Pst as an aproximationn to Qst (phenotypic distance) and comparison with Fst using Bayesian GLMM approach. [command line version]
- Pst_mcmc - Pst as an aproximationn to Qst (phenotypic distance) and comparison with Fst using Bayesian GLMM approach. [desktop version]
- Run and Plot RDA - not a ready to use script, needs better anotation and to be adjusted for each input file and option.
- average SD - Quick script that will output a table of summary statistics (quartiles, average, SD, etc.) from multiple files 
- basicgen Fst - Calculates pairwise Fst values with confidence intervals and p-values [command line version]
- basicgen diverindex - Calculates various allelic diversity indexes using hierfsts, diveRsity and dartR, outputs a table per population [command line version] 
- pairLDboxplot - categorize plink LD r-squared table in distance (kb) groups, useful to run in server in order to save time [command line version]
- pairLDdist - get summary tables from the scaffold sizes, distance, number of loci and LD from Plink LD analysis output, useful to run it in the server to save time [command line version]
- pairLDdist_windows - if your computer is fast enough skip the command line version of the script and directly calculate and plot from desktop R
- regression ressiduals - calculate linear and quadratic regression from multiple files, save plots, residuals, and summary table

# I'll be adding more complex scripts here once they are ready to be used by anyone (properly annotated and more broad options)


