#Quick script to plot per-locus Fst pair-wise values from populations (Stacks) output

#populations Fst import
setwd("")
# perl one liner to extract only relevant columns:
# perl -lne '/^(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t.*?$/; print "$1\t$2\t$3\t$8\t$9"' populations.fst_PK-PM.tsv > fstlist
#if some error, try changing EOL to windows
# if some error, try changing "\t" to ","
# you may need to rename column names
dataset = read.table("fstlist", sep = ',', header =TRUE)

nrow(dataset)

#change names
str(dataset)
class(dataset)
FSTs <- as.data.frame(dataset)   			#transform to dataframe if needed
names(FSTs) <- c("locus", "pop1", "pop2", "Fst", "pvalue")

#Manhattan plot!
ggplot(FSTs, aes(x=FSTs$locus, y=FSTs$Fst))+
geom_point(alpha=0.8, size=1.3)+
theme(axis.text.x=element_blank()) +
labs(title = "Fst values per loci calculated by populations (Stacks) from PM and PK", x = "Loci", y = "Fst") +
theme(plot.title = element_text(hjust=0.5, size = 18))


alpha = 0.001
FSTs$significant <- "black"
FSTs$significant[FSTs$pvalue < alpha] <- "red3"


#Manhattan plot with color!

image=
ggplot(FSTs, aes(x=FSTs$locus, y=FSTs$Fst))+
geom_point(alpha=0.8, size=1.3, color=FSTs$significant)+
theme(axis.text.x=element_blank()) +

theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +

labs(title = "Fst values per locus calculated by populations (Stacks) from PM and PK", x = "Loci", y = "Fst") +
theme(plot.title = element_text(hjust=0.5, size = 18))
loadplotfile="PKPM_Fsts.png"
ggsave(file=loadplotfile, plot=image, width=15, height=5, dpi=350)





groups <- factor(PKPM[[2]])
groups
class(groups)
pops <- as.data.frame(groups)
length(groups)
odd_indexes<-seq(1,67,2)
poptags <- pops[odd_indexes,1]
tags <- as.factor(poptags)
tags
plot(PCApkpm$x[,1:2], col = tags, pch = 19)





# COOL
titleplot <- "PCA of PK and PM allelic frequencies"
titleplot
titlepca <- paste("PKPM_poda", "_PCA", "_autoplot", ".svg", sep="")
prop.pca = PCApkpm$sdev^2/sum(PCApkpm$sdev^2)
prop.pca
Xlabel <- paste("PC1 (", round(100*(prop.pca[1]), digits=2), "%)", sep="")
Ylabel <- paste("PC2 (", round(100*(prop.pca[2]), digits=2), "%)", sep="")

image=
autoplot(PCApkpm, data=PCAmatrix4) +
labs(x = Xlabel, y = Ylabel) + coord_fixed(ratio = 1) +
#theme(axis.line = element_line(colour = "gray70"), panel.background = element_blank(), legend.key=element_blank(), legend.background = element_rect(linetype="solid", colour ="gray70")) +
theme_bw() +
geom_point(aes(shape = tags, color=tags, fill = tags), size = 2.5) +
scale_fill_manual(values=mypalette)+
scale_colour_manual(values=mypalette) +
scale_shape_manual(values=myshapes)
