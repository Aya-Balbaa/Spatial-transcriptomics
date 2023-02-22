library(limma)
library(edgeR)

##microglua in NLF
AP=read.csv("AP-IBA.csv")
row.names(AP)=AP$rownames
AP=AP[,-1]
##design
group = factor(rep(c(rep("a", 6),rep("pp", 6), rep("p", 6))))
pairinfo = factor(rep(1:6,3))
design=model.matrix(~0+pairinfo+group)


d <- DGEList(AP)
d2=calcNormFactors(d)
y <- voom(d2, design, plot = T)
fit <- lmFit(y, design)
tmp <- eBayes(fit)

##to compare away vs plaque 
allDiff_pair <- topTable(
  tmp,adjust = 'BH',
  coef = 7,
  number = Inf)
write.csv(allDiff_pair,"plaque vs away microglia in NLF-.csv")

##to compare away vs peri-plaque 
allDiff_pair <- topTable(
  tmp,adjust = 'BH',
  coef = 8,
  number = Inf)
write.csv(allDiff_pair,"peri-plaque vs away microglia in NLF-.csv")

##Astrocytes
gf=read.csv("PA-GFAP-NLF.csv")
row.names(gf)=gf$rownames
gf=gf[,-1]
d <- DGEList(gf)
d2=calcNormFactors(d)
pairinfo = factor(rep(1:6,2))
group = factor(rep(c(rep("A",6 ), rep("p", 6))))
design=model.matrix(~ 0+pairinfo+group)
y <- voom(d2, design, plot = T)
fit <- lmFit(y, design)
tmp <- eBayes(fit)
allDiff_pair <- topTable(
  tmp,adjust = 'BH',
  coef = 7,
  number = Inf)
write.csv(allDiff_pair,"plaque vs away in astrocytes NLF.csv" )
