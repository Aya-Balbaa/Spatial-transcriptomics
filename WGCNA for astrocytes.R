library(WGCNA)
library(flashClust)
library(DESeq2)
g=cbind(gf,tgf)
g2=round(g)
batchg=factor(c(1,1,1,2,3,3, 1,1,1,2,3,3,1,2,2,3,1,2,2,3))
condition=factor(c(rep("A", 6),rep("P", 6), rep("AT", 4), rep("PT", 4)))

sampleTable <- data.frame(condition =condition, batchg=batchg)
dds2<-DESeqDataSetFromMatrix(countData=g2, 
                             colData= sampleTable, 
                             design= ~batchg+condition, tidy = FALSE)
dds=DESeq(dds2)
design4= model.matrix(~0+condition)
rlv=vst(dds, blind=FALSE)
rl=assay(rl)
rb <- limma::removeBatchEffect(rl, batchg, design=design4)

##detecting outliers samples
trb=t(rb)
sampleTree2 = hclust(dist(trb), method = "average")
##sample 6 is an outlier

rb6=rb[,-6]

madg=apply(rb6, 1, mad)
mad1=as.matrix(madg)
varb2=mad1[order(mad1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
write.csv(varb2, "mad-gfap.csv")
varb3=varb2[1:5000,]
varb3=as.matrix(varb3)
g3=c(row.names(varb3))
rbg=rb6[g3,]
trbg=t(rbg)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(trbg6, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power graph
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net = blockwiseModules(trbg, power = 7,
                       TOMType = "signed", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "gfapMouseTOM",
                       verbose = 3)
ModuleColors2 = labels2colors(net$colors)
nGenes = ncol(trbg);
nSamples = nrow(trbg);
MEs0 = moduleEigengenes(trbg6, ModuleColors2)$eigengenes


NLF= factor(c(rep("-Away",5 ), rep("-Peri-Plaque", 6), rep("Trem2-Away", 4), rep("Trem2-Peri-plaque", 4)))
design3=model.matrix(~0+NLF)

moduleTraitCor = cor(MEs, design3, use = "p");
write.csv(moduleTraitCor, "GFAP-batch-A6-5000-modultrait.csv")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.csv(moduleTraitPvalue, "GFAP-modules-pvalue.csv")
# to display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design3),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))