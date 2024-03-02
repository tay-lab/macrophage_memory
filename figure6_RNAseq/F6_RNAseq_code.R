#setwd()
library(limma)
library(edgeR)
library(dplyr)
library(gplots)
library(org.Mm.eg.db)
library(ggplot2)
library("RColorBrewer")
library(EGSEA)

#Preprocessing from count matrix
DGEobj = read.csv("merged_FC.csv", header= TRUE)
group=as.factor( c('CF', 'CF', 'CF', 'CL', 'CL', 'CL', 'CT', 'CT', 'CT', 
                         'FC', 'FC', 'FC', 'FF', 'FF', 'FF', 'FL', 'FL', 'FL',
                         'FP', 'FP', 'FP', 'FT', 'FT', 'FT', 'PF', 'PF', 'PF', 
                         'PL', 'PL', 'PL', 'PT', 'PT', 'PT') )

DGEobj = tibble::column_to_rownames(DGEobj[-1], var="gene")
DGEobj = DGEList(counts = DGEobj, group=group)
keep <- filterByExpr(DGEobj, min.count = 25, group=group)
DGEobj <- DGEobj[keep,,keep.lib.sizes=FALSE]
DGEobj <- calcNormFactors(DGEobj)

logCPM <- cpm(DGEobj, log=TRUE, prior.count = 3)
geneid <- rownames(logCPM)
genes <- AnnotationDbi::select(org.Mm.eg.db, keys = geneid, column = c('SYMBOL'), keytype = "ENTREZID")


rownames(logCPM) <- genes$SYMBOL
colnames(logCPM) <-  c('CF', 'CF', 'CF', 'CL', 'CL', 'CL', 'CT', 'CT', 'CT', 
                       'FC', 'FC', 'FC', 'FF', 'FF', 'FF', 'FL', 'FL', 'FL',
                       'FP', 'FP', 'FP', 'FT', 'FT', 'FT', 'PF', 'PF', 'PF', 
                       'PL', 'PL', 'PL', 'PT', 'PT', 'PT') 
write.csv(logCPM, "logCPM.csv")

design <- model.matrix(~0+group)
contr.matrix <- makeContrasts(CF = groupCF - groupFF, 
                              CL = groupCL - groupFF, 
                              CT = groupCT - groupFF, 
                              FC = groupFC - groupFF, 
                              FL = groupFL - groupFF, 
                              FP = groupFP - groupFF, 
                              FT = groupFT - groupFF, 
                              PF = groupPF - groupFF, 
                              PL = groupPL - groupFF,
                              PT = groupPT - groupFF,
                              FPFC = groupFP - groupFC,
                              levels = colnames(design))

fit <- lmFit(logCPM, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit)
dt <- decideTests(efit, adjust.method = "BH", p.value = 0.01,lfc = 1)

#IDing synergy and antagonism
deltaHistoryUP = data.frame(deltaPICTNF=double(), deltaPICLPS=double(), 
                              deltaCPGTNF=double(), deltaCPGLPS=double())

deltaHistoryDOWN = data.frame(deltaPICTNF=double(), deltaPICLPS=double(), 
                              deltaCPGTNF=double(), deltaCPGLPS=double())

for (aa in 1:length(foo2)) 
{
  temp = logCPM[foo2[aa], ]
  if ( any( dt[foo2[aa], c(7:8, 10) ] > 0 ) ) {
    deltaPICTNF = ( mean(2^temp[31:33])  - ( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) )+ mean(2^temp[13:15])  ), 1 )#mean(2^temp[13:15])
    deltaHistoryUP[aa, 1] = deltaPICTNF
    
    deltaHistoryDOWN[aa, 1] = 0

  } else if (any( dt[foo2[aa], c(7:8, 10) ] < 0)) {
    deltaPICTNF = ( mean(2^temp[31:33])  - ( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( ( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ), 1))#mean(2^temp[13:15])
    deltaHistoryDOWN[aa, 1] = deltaPICTNF
    deltaHistoryUP[aa, 1] = 0
  } else {
    deltaPICTNF = 0
    deltaHistoryUP[aa, 1] = deltaPICTNF
    deltaHistoryDOWN[aa, 1] = deltaPICTNF
  }
  
  
  if ( any( dt[foo2[aa], c(5, 8, 9) ] > 0 ) ) {
    deltaPICLPS = ( mean(2^temp[28:30])  - ( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ), 1)#mean(2^temp[13:15])
    deltaHistoryUP[aa, 2] = deltaPICLPS
    deltaHistoryDOWN[aa, 2] = 0
  } else if (any( dt[foo2[aa], c(5, 8, 9) ] < 0)) {
    deltaPICLPS = ( mean(2^temp[28:30])  - ( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( ( (mean(2^temp[25:27]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ), 1) )#mean(2^temp[13:15])
    deltaHistoryDOWN[aa, 2] = deltaPICLPS
    deltaHistoryUP[aa, 2] = 0
  } else {
    deltaPICLPS = 0
    deltaHistoryUP[aa, 2] = deltaPICLPS
    deltaHistoryDOWN[aa, 2] = deltaPICLPS
  }

  if ( any( dt[foo2[aa], c(1, 3, 7)  ] > 0 ) ) {
    deltaCPGTNF = ( mean(2^temp[7:9])  - ( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) )+ mean(2^temp[13:15]) ), 1 )#mean(2^temp[13:15])
    deltaHistoryUP[aa, 3] = deltaCPGTNF
    deltaHistoryDOWN[aa, 3] = 0
  } else if (any( dt[foo2[aa], c(1, 3, 7) ] < 0)) {
    deltaCPGTNF = ( mean(2^temp[7:9])  - ( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( ( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[22:24]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ), 1) )#mean(2^temp[13:15])
    deltaHistoryDOWN[aa, 3] = deltaCPGTNF
    deltaHistoryUP[aa, 3] = 0
  } else {
    deltaCPGTNF = 0
    deltaHistoryUP[aa, 3] = deltaCPGTNF
    deltaHistoryDOWN[aa, 3] = deltaCPGTNF
  }

  if ( any( dt[foo2[aa], c(1, 2, 5)  ] > 0 ) ) {
    deltaCPGLPS = ( mean(2^temp[4:6])  - ( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ), 1)#mean(2^temp[13:15])
    deltaHistoryUP[aa, 4] = deltaCPGLPS
    deltaHistoryDOWN[aa, 4] = 0
  } else if (any( dt[foo2[aa], c(1, 2, 5) ] < 0)) {
    deltaCPGLPS = ( mean(2^temp[4:6])  - ( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ) )/max( c( ( (mean(2^temp[1:3]) - mean(2^temp[13:15])) + (mean(2^temp[16:18]) - mean(2^temp[13:15]) ) + mean(2^temp[13:15]) ), 1 ) )#mean(2^temp[13:15])
    deltaHistoryDOWN[aa, 4] = deltaCPGLPS
    deltaHistoryUP[aa, 4] = 0
  } else {
    deltaCPGLPS = 0
    deltaHistoryUP[aa, 4] = deltaCPGLPS
    deltaHistoryDOWN[aa, 4] = deltaCPGLPS
  }
  
}


###PT hmaps
deltaPICTNF = deltaHistoryUP[deltaHistoryUP$deltaPICTNF > .25, ]
deltaPICTNF = rownames( deltaPICTNF[ order(deltaPICTNF$deltaPICTNF, decreasing = TRUE), ] )
write.csv(deltaPICTNF, file="PT_1.csv")

deltaPICTNF = deltaHistoryUP[deltaHistoryUP$deltaPICTNF < -.25, ]
deltaPICTNF = rownames( deltaPICTNF[ order(deltaPICTNF$deltaPICTNF, decreasing = FALSE), ] )
write.csv(deltaPICTNF, file="PT_2.csv")


p <- coolmap(logCPM[deltaPICTNF, c( 31:33, 25:27, 22:24, 13:15) ], linkage.row="ward.D", linkage.col="none", 
             show.dendrogram="row", xlab=NULL, ylab=NULL, breaks=seq(-2.5, 2.5, length.out=257))

###PL hmaps###
deltaPICLPS = deltaHistoryUP[deltaHistoryUP$deltaPICLPS > .25, ]
deltaPICLPS = rownames( deltaPICLPS[ order(deltaPICLPS$deltaPICLPS, decreasing = TRUE), ] )
write.csv(deltaPICLPS, file="PL_1.csv")

deltaPICLPS = deltaHistoryUP[deltaHistoryUP$deltaPICLPS < -0.25, ]
deltaPICLPS = rownames( deltaPICLPS[ order(deltaPICLPS$deltaPICLPS, decreasing = FALSE), ] )
write.csv(deltaPICLPS, file="PL_2.csv")


p <-coolmap(logCPM[deltaPICLPS, c( 28:30, 16:18, 25:27, 13:15)], linkage.row="ward.D", linkage.col="none", show.dendrogram="row",
            xlab=NULL, ylab=NULL, breaks=seq(-2.5, 2.5, length.out=257))


###CT hmaps###
deltaCPGTNF = deltaHistoryUP[deltaHistoryUP$deltaCPGTNF > 0.25, ]
deltaCPGTNF = rownames( deltaCPGTNF[ order(deltaCPGTNF$deltaCPGTNF, decreasing = TRUE), ] )
write.csv(deltaCPGTNF, file="CT_1.csv")

deltaCPGTNF = deltaHistoryUP[deltaHistoryUP$deltaCPGTNF < -0.25, ]
deltaCPGTNF = rownames( deltaCPGTNF[ order(deltaCPGTNF$deltaCPGTNF, decreasing = FALSE), ] )
write.csv(deltaCPGTNF, file="CT_2.csv")


p = coolmap(logCPM[deltaCPGTNF, c( 7:9, 1:3, 22:24, 13:15)], linkage.row="ward.D", linkage.col="none", 
            show.dendrogram="row",breaks=seq(-2.5, 2.5, length.out=257) )


###CL hmaps###
deltaCPGLPS = deltaHistoryUP[deltaHistoryUP$deltaCPGLPS > 0.25, ]
deltaCPGLPS = rownames( deltaCPGLPS[ order(deltaCPGLPS$deltaCPGLPS, decreasing = TRUE), ] )
write.csv(deltaCPGLPS, file="CL_1.csv")

deltaCPGLPS = deltaHistoryUP[deltaHistoryUP$deltaCPGLPS < -0.25, ]
deltaCPGLPS = rownames( deltaCPGLPS[ order(deltaCPGLPS$deltaCPGLPS, decreasing = FALSE), ] )
write.csv(deltaCPGLPS, file="CL_2.csv")


p = coolmap(logCPM[deltaCPGLPS, c( 4:6, 16:18,  1:3, 13:15)], linkage.row="ward.D", linkage.col="none", show.dendrogram="row",
            xlab=NULL, ylab=NULL, breaks=seq(-2.5, 2.5, length.out=257))

## Supp figure 4E
foo3 = rownames( dt[ dt[,"FPFC"] != 0, ] )
foo3 = foo3[!is.na(foo3)]
write.csv(foo3, file="FPvsFC.csv")

gs.annots = buildIdx(entrezIDs = rownames(DGEobj), species = "mouse", msigdb.gsets = "none",
                     kegg.exclude = c("Metabolism"))


test = egsea.cnt(counts = DGEobj, group = group, design = design,
            contrasts = contr.matrix, gs.annots = gs.annots, symbolsMap = genes,
            baseGSEAs = egsea.base()[-c(2, 12)], sort.by = "avg.rank",
            num.threads = 4, report = FALSE)

topSets(test, contrast = 11, gs.label = "kegg", number = 15)
plotPathway(test, contrast = 11, "NF-kappa B signaling pathway", gs.label = "kegg", file.name = "NF-kappa B signaling pathway")
plotPathway(test, contrast = 11, "TNF signaling pathway", gs.label = "kegg", file.name = "TNF signaling pathway")
coolmap(logCPM[foo3, c( 19:21, 10:12, 13:15)], linkage.row="ward", linkage.col="none", show.dendrogram="row",
        xlab=NULL, ylab=NULL, breaks=seq(-2.5, 2.5, length.out=257))
