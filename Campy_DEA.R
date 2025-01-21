# RNA-Seq DE analysis of Campylobacter jejuni mapped to reference genome
# data from Manca Volk <Manca.Volk@bf.uni-lj.si>
# running @ local machine
# 21 counts files in total

# Load and dcombine all counts into an expression matrix
setwd("./pISA-Projects/_p_ExtDataAnalysis/_I_CampyBF/_S_CampyTS/_A_03_limma-R/scripts")
files <- list.files(path="../../_A_02_STAR-RNASeq/output", pattern=".ReadsPerGene.out.tab",  full.names= TRUE)

geneNames <- read.table(files[1], header = FALSE, sep = "\t", skip= 4)[,1]

counts <- NULL
colNs <- NULL
counti <- NULL

for(i in 1:length(files)) {
  counti <- read.table(files[i], header = FALSE, sep = "\t", skip= 4)[,4] #  row 4 are reverse stranded reads
  colNi <- unlist(strsplit(unlist(strsplit(files[i], ".", fixed = TRUE))[5], "/", fixed = TRUE))[4]
  
  counts <- cbind(counts, counti)
  colNs <- cbind(colNs , colNi) 
}

colNs

colnames(counts) <- as.vector(colNs)
rownames(counts) <- stringr::str_remove(geneNames, "gene:")

head(counts) # check if it makes sense!

#export table with all raw counts
write.table(counts, file="../output/Campy2022_STAR_counts.txt", sep="\t", row.names=TRUE)

# Load packages

library ("limma")
library("edgeR")
library("stringr")


colnames(counts)

# Read phenodata (sample metadata used for RNA-Seq)
phenodata <- read.table("../../../phenodata_20221005.txt",  row.names=1, header = TRUE, sep = "\t")
dim(phenodata)
head(phenodata)

#check if samples are in same order; if TRUE it's OK :)
all(rownames(phenodata)==colnames(counts))

x <- counts

# assign groups to samples
group <- as.factor(phenodata$SampleGroup)

## Create a DGEList object for limma statistical analysis and specify library size i.e. number of reads sequenced per sample
# default library sizes (sum of mapped reads)!!
y <- DGEList(counts=x, group=group)

#######################################################################################################################################
#unfiltered TMM-normalized results output
y_GSEA <- calcNormFactors(y)
y_GSEA <- cpm(y_GSEA, log=TRUE, prior.count=0.5)
y_GSEA <- cbind(rownames(y_GSEA), rep("NA", length(rownames(y_GSEA))), y_GSEA)
colnames(y_GSEA)[1:2] <- c("NAME","DESCRIPTION")

head(y_GSEA)  
write.table(y_GSEA, file="../output/RNAseq_TMM_GSEA_noDESC.txt", sep="\t", quote=FALSE, row.names=FALSE)
#add DESCRIPTIONS in Excel if needed

#########
#remove CP_Pf_Mx from object y as it will not be in contrasts
y <- y[, -21]

# assign groups to samples, also here remove CP_Pf_Mx sample
group <- as.factor(phenodata$SampleGroup[-21])

# Create design matrix
design <- model.matrix(~0+y$samples$group)
colnames(design) <- levels(y$samples$group)
rownames(design) <- rownames(y$samples)

# Filter low expressed genes
keep <- filterByExpr(y, design)

# Normalization of dataset for different library sizes
y1 <- y[keep, ]
y1 <- calcNormFactors(y1)

#######################################################################################################################################
# QC plots

#define color palette for the samples in graphs
col = as.numeric(y$samples$group) +1

#check density plot
pdf("../other/density_plot_before_lowexpr_filter.pdf", onefile=TRUE, family= "Helvetica")
nsamples <- ncol(x)

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.6), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
  
}
legend("topright", colnames(lcpm), text.col=col, cex=0.8, bty="n")
dev.off()
# comment: looks OK

## Plot QC plots using different functions e.g.:

pdf("../other/log10rawcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
opar <- par()
par(cex.lab=1, cex.axis=1)
boxplot(log(y$counts+1,10), las=2, ylab="log10(counts)", col=col, cex=0.5)
par(opar)
dev.off()

pdf("../other/log10filteredcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
opar <- par()
par(cex.lab=1, cex.axis=1)
boxplot(log(y1$counts+1,10), las=2, ylab="log10(counts)", col=col, cex=0.5)
par(opar)
dev.off()

pdf("../other/TMMnormLOGcpm_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(cpm(y1, log=TRUE, prior.count=0.5), las=2, ylab="TMM normalized log_cpm values with prior.counts 0.5", col=col)
dev.off()

#####

#density plots before and after removing low expressed genes
pdf("../other/norm_counts_raw&filtered_densityplots.pdf", onefile=TRUE, family= "Helvetica")
opar <- par()
par(mfrow=c(1,2), cex = 0.6)
nsamples <- ncol(x)

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.6), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", colnames(lcpm), text.col=col, bty="n")

lcpm <- log(as.matrix(y1),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.6), las=2, main="", xlab="")
title(main="B. AFTER REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(lcpm), text.col=col, bty="n", cex=0.6)
par(opar)
dev.off()

#####

#density plots before and after removing low expressed genes
pdf("../other/TMMnorm_counts_raw&filtered_densityplots.pdf", onefile=TRUE, family= "Helvetica")
opar <- par()
par(mfrow=c(1,2), cex = 0.6)
nsamples <- ncol(x)

lcpm <- cpm(x, log=TRUE, prior.count=0.5)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.2), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", colnames(lcpm), text.col=col, bty="n")

lcpm <- cpm(y1, log=TRUE, prior.count=0.5)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.2), las=2, main="", xlab="")
title(main="B. AFTER REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(lcpm), text.col=col, bty="n", cex=0.6)
par(opar)
dev.off()

#MDS (PCA-like) graph
pdf("../other/norm_counts_filtered_MDS.pdf", onefile=TRUE, family= "Helvetica")
plotMDS(y1, labels=colnames(y1), col = col, cex = 0.5)
legend("topright", levels(y1$samples$group), text.col=c(1:7), bty="o", cex=0.8)
dev.off()

#TMM-normalized log cpm
pdf("../other/norm_lcpm_MDS.pdf", onefile=TRUE, family= "Helvetica")
plotMDS(cpm(y1, log=TRUE, prior.count=0.5), labels=colnames(y1), col = col, cex = 0.5)
legend("topright", levels(y1$samples$group), text.col=c(1:7), bty="o", cex=0.8)
dev.off()

#######################################################################################################################################

# Define contrasts i.e. comparisons between groups:

#create design matrix
colnames(design) <- c("CJ16H", "CJ16P", "CJ24H", "CJ48H", "CJ72H")
design

contrastMatrix = makeContrasts("CJ16H-CJ16P",
                               "CJ24H-CJ16P",
                               "CJ48H-CJ16P",
                               "CJ72H-CJ16P",
                               levels=design)

###########
## limma-trend function

logCPM <- cpm(y1, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit2 <- contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit2, trend=TRUE)

## make results table
results1 <- topTable(fit2, coef=1, number=1000000, sort.by="none")
results2 <- topTable(fit2, coef=2, number=1000000, sort.by="none")
results3 <- topTable(fit2, coef=3, number=1000000, sort.by="none")
results4 <- topTable(fit2, coef=4, number=1000000, sort.by="none")

#make one expression matrix with logFCs and adj.P.vals
results <- cbind(results1[,1], results1[,5], results2[,1], results2[,5], results3[,1], results3[,5], results4[,1], results4[,5])

colnames(results) <- c(paste(colnames(contrastMatrix)[1], " ", colnames(results1[1])),
                       paste(colnames(contrastMatrix)[1], " ", colnames(results1[5])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[1])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[5])),
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[1])),
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[5])),
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[1])),
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[5]))
)

rownames(results) <- rownames(results1)
# add raw expression data
length(rownames(y1))==length(results[,1])

# have to do merge, not cbind!
results.raw <- merge(results, y1$counts, by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.raw)
write.table(results.raw, file="../output/RNAseq_logFC_padj_rawReads_TREND.txt", sep="\t", quote=TRUE, row.names=FALSE)

#instead of raw counts here I export TMM normalized cpm values with prior.counts 0.5 
results.tmm <- merge(results, cpm(y1, log=TRUE, prior.count=0.5), by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.tmm)
write.table(results.tmm, file="../output/RNAseq_logFC_padj_TMMcpm_TREND.txt", sep="\t", quote=TRUE, row.names=FALSE)




###########
## limma-voom function (alternative)

# limma voom fit for filtered RNA-seq dataset (y1)
pdf("../other/voom_mean-variance_trend.pdf", onefile=TRUE, family= "Helvetica")
v <- voom(y1,design,plot=TRUE)
dev.off()

fit <- lmFit(v, design)
fit2 = contrasts.fit(fit, contrastMatrix)

## eBayes statistics calculation
fit2 <- eBayes(fit2)
pdf("../other/SIGMA_vs_A_plot.pdf", onefile=TRUE, family= "Helvetica")
plotSA(fit2)
dev.off()

## make results table
results1 <- topTable(fit2, coef=1, number=1000000, sort.by="none")
results2 <- topTable(fit2, coef=2, number=1000000, sort.by="none")
results3 <- topTable(fit2, coef=3, number=1000000, sort.by="none")
results4 <- topTable(fit2, coef=4, number=1000000, sort.by="none")

#make one expression matrix with logFCs and adj.P.vals
results <- cbind(results1[,1], results1[,5], results2[,1], results2[,5], results3[,1], results3[,5], results4[,1], results4[,5])

colnames(results) <- c(paste(colnames(contrastMatrix)[1], " ", colnames(results1[1])),
                       paste(colnames(contrastMatrix)[1], " ", colnames(results1[5])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[1])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[5])),
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[1])),
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[5])),
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[1])),
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[5]))
)

rownames(results) <- rownames(results1)
# add raw expression data
length(rownames(y1))==length(results[,1])
# have to do merge, not cbind!

results.raw <- merge(results, y1$counts, by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.raw)
write.table(results.raw, file="../output/RNAseq_logFC_padj_rawReads_VOOM.txt", sep="\t", quote=TRUE, row.names=FALSE)

#instead of raw counts here I export TMM normalized cpm values with prior.counts 0.5 
results.tmm <- merge(results, cpm(y1, log=TRUE, prior.count=0.5), by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.tmm)
write.table(results.tmm, file="../output/RNAseq_logFC_padj_TMMcpm_VOOM.txt", sep="\t", quote=TRUE, row.names=FALSE)

sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Slovenian_Slovenia.1250  LC_CTYPE=Slovenian_Slovenia.1250    LC_MONETARY=Slovenian_Slovenia.1250 LC_NUMERIC=C                       
# [5] LC_TIME=Slovenian_Slovenia.1250    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] stringr_1.4.0 edgeR_3.36.0  limma_3.50.3 
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9        magrittr_2.0.3    splines_4.1.1     usethis_2.1.5     devtools_2.4.3    pkgload_1.2.4     lattice_0.20-44   R6_2.5.1         
# [9] rlang_1.0.2       fastmap_1.1.0     tools_4.1.1       grid_4.1.1        pkgbuild_1.3.1    sessioninfo_1.2.2 cli_3.3.0         withr_2.5.0      
# [17] ellipsis_0.3.2    remotes_2.4.2     rprojroot_2.0.3   lifecycle_1.0.1   crayon_1.5.1      brio_1.1.3        processx_3.5.3    purrr_0.3.4      
# [25] callr_3.7.0       fs_1.5.2          ps_1.6.0          testthat_3.1.3    memoise_2.0.1     glue_1.6.2        cachem_1.0.6      stringi_1.7.6    
# [33] compiler_4.1.1    desc_1.4.1        prettyunits_1.1.1 locfit_1.5-9.5   
