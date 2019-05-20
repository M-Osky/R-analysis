# R analysis

Those are old R scripts I use for various purpouses.
I strongly recommend to use the ones in the parent directory, but it may be worth to check this two if you need simpler things.

Basically you have:
one-way ANOVA and t-test for Normal data
Kruskall-Wallis and Wilcoxon test for non-Normal
And multivariant analysis: PCA and LDA
They include a few lines about the input file needed and so, so should be quite straigh forward to use them.
Usually columns should be sorted as this: tags/descriptors (not analised), grouping variables, numerical variables.

		id	loccode		pop	sex	var1	var2	var3	var4	var5	var6
		Lem42	LemuriaNA	LM	fem	004	008	015	016	023	042
		Atl16	AtlantisSE	AT	fem	003	014	015	092	062	036
		Atl28	Atlantida	AT	mal	013	021	034	055	089	144

I'm still working (when I have time) with the RMark scripts, recently found out I could improve some of the models (they weren't 100% correct), as soon as I update them, they will be here. Feedback will be strongly appreciated.
Cheers!
