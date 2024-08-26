#setwd("")
library(limma)
library(edgeR)
library(dplyr)
library(gplots)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- Mmusculus
library(monaLisa)
library(JASPAR2022)
library(TFBSTools)
pwms <- getMatrixSet(JASPAR2022,
                     opts = list(matrixtype = "PWM",
                                 tax_group = "vertebrates"))
library(SummarizedExperiment)
library(stringr)
library(stringi)
library(pheatmap)
set.seed(1)

DGEobj = read.csv("ATAC_FC.csv", header= TRUE)
group=as.factor( c('CL', 'CL', 'CT','CT', 'MC', 'MC','ML', 'ML','MM', 'MM', 
                   'MP', 'MP', 'MT','MT', 'PL', 'PL', 'PT', 'PT') )
                   
DGEobj = tibble::column_to_rownames(DGEobj, var="X")
rownames = rownames(DGEobj)
DGEobj = DGEList(counts = DGEobj, group=group)
keep <- filterByExpr(DGEobj, min.count = 10, group=group)
DGEobj <- DGEobj[keep,,keep.lib.sizes=FALSE]
DGEobj <- calcNormFactors(DGEobj)
logCPM <- cpm(DGEobj, log=TRUE, prior.count = 3)
logCPMRNA <- read.csv("../figure6_RNAseq/logCPM.csv")

rm(DGEobj)
colnames(logCPM) <- c('CL1', 'CL2', 'CT1','CT2', 'MC1', 'MC2','ML1', 'ML2','MM1', 'MM2', 
                      'MP1', 'MP2', 'MT1','MT2', 'PL1', 'PL2', 'PT1', 'PT2') 

write.csv(logCPM, "supp_table_3_ATAClogCPM.csv")

bedfile = import("../08mergedpeaks/All_Samples.fwp.filter.non_overlapping.bed", format = "BED")

design <- model.matrix(~0+group)
contr.matrix <- makeContrasts(MC = groupMC - groupMM, 
                              MP = groupMP - groupMM, 
                              PT = groupPT - groupMM, 
                              CT = groupCT - groupMM, 
                              PL = groupPL - groupMM, 
                              CL = groupCL - groupMM,  
                              MT = groupMT - groupMM, 
                              ML = groupML - groupMM, 
                              levels = colnames(design))
fit <- lmFit(logCPM, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit)
dt <- decideTests(efit, adjust.method = "BH", p.value = 0.01,lfc = 0.585)
summary(dt)

#MC or MP diff
MCMPdiff <- rownames( dt[dt[, 1] != 0 | dt[, 2] != 0, ]  )

p = coolmap(logCPM[MCMPdiff, c(9, 10, 5, 6, 11, 12)], 
            linkage.row="ward.D", linkage.col="none", show.dendrogram="row",
            xlab=NULL, ylab=NULL, breaks=seq(-1.5, 1.5, length.out=257) )
q = as.hclust(p$rowDendrogram)
geneclusters <- cutree(q, k=5)

coolmap(logCPM[names(geneclusters[ geneclusters==2]) , c(9, 10, 5, 6, 11, 12)],
        linkage.row="ward", show.dendrogram="row" )

grouppeaks <- unlist( unlist(strsplit(names(geneclusters[ geneclusters==4]), "[.]") ) )
grouppeaks <- as.data.frame(matrix(grouppeaks,ncol =3,byrow = T))
colnames(grouppeaks) <- c("seqnames", "start", "end")
grouppeaks$start = as.numeric(grouppeaks$start)
grouppeaks$end = as.numeric(grouppeaks$end)
peakseqs <- GRanges(seqnames = grouppeaks$seqnames, ranges = IRanges(start = grouppeaks$start, end = grouppeaks$end))

seq3 <- getSeq(genome, peakseqs3)
names(seq3) <- c(names(geneclusters[ geneclusters==3]))

macsTFs = import("Oth.Bld.20.AllAg.Macrophages.bed", format = "BED")
testTF = str_detect(macsTFs$name, "Name=Irf8%")
testTF2 = str_detect(macsTFs$name, "Name=Irf9%")
testTF3 = str_detect(macsTFs$name, "Name=Irf3%")
testTF = testTF|testTF2|testTF3
subsetTFs = macsTFs[testTF]
testTFchip <- findOverlaps(peakseqs, subsetTFs, ignore.strand = TRUE, select = "first")
peakseqsnew <- peakseqs[!is.na(testTFchip)]
length(peakseqsnew)
length(peakseqsnew)/length(peakseqs)


#General motif analysis code
peakAnno <- annotatePeak(peakseqs, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(peakAnno)
foo <- as.data.frame(peakAnno)
peakids <- GRanges(seqnames = foo$seqnames, ranges = IRanges(start = foo$start, end = foo$end),
                   strand = foo$strand)
peakseqs <- as.GRanges(peakAnno)
peakseqs <- intersect(peakseqs, peakids, ignore.strand = TRUE)

seq <- getSeq(genome, peakseqs)
names(seq) <- c(1:length(seq))
pwms <- getMatrixSet(JASPAR2022,
                    opts = list(matrixtype = "PWM",
                                tax_group = "vertebrates"))
se <- calcBinnedMotifEnrR(seqs = seq,
                          pwmL = pwms,
                          background = "genome",
                          genome = BSgenome.Mmusculus.UCSC.mm10,
                          genome.regions = NULL, # sample from full genome
                          genome.oversample = 2,
                          BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
                          verbose = TRUE)
se_expr <- se[which( tolower( se@elementMetadata@listData$motif.name) %in% tolower( logCPMRNA$X ) ), ]
sel <- assay(se_expr, "log2enr")[, 1] > 1 #& assay(se, "negLog10Padj")[, 1] > 4.0
sel[is.na(sel)] = FALSE
plotMotifHeatmaps(x = se_expr[sel,][c( 13, 12, 2, 1, 15, 26, 25, 27, 33)], which.plots = c("log2enr", "negLog10Padj"), 
                  width = 0.5, width.seqlogo = 1, maxEnr = 4, maxSig = 25, 
                  show_seqlogo = TRUE)
write.csv(cbind(assays(se_expr[sel,])$log2enr, assays(se_expr[sel,])$negLog10Padj ), "PLCL_group5_motiftable.csv")
#For PTCT group 1, FC > 1
#[c( 2, 8, 9, 19, 10, 3, 22, 36)]
#For PTCT group 2, FC > 1
#[c( 13, 14, 3, 8, 9, 17, 36, 4)]
#For PTCT group 3, FC > 1
#[c( 8, 7, 1, 10, 23, 14, 25, 27)]

#For PLCL group 1, FC > 1
#[c( 11, 12, 2, 8, 7, 14, 32, 35, 28)]
#For PLCL group 3, FC > 1
#[c( 2, 8, 9, 18, 12, 15, 21, 28, 3)]
#For PLCL group 5, FC > 1
#[c( 13, 12, 2, 1, 15, 26, 25, 27, 32)]

#for group 3, FC > 1
#[c( 14, 15, 3, 2, 10, 11, 35, 29, 28)]
#for group 2, FC > 1
#[c( 14, 15, 3, 2, 9, 10, 35, 29, 4)]
#for group 4, FC > 1
#[c( 9, 10, 1, 13, 24, 21, 17, 7, 35)]
#for group 1, FC > 1.5
#all
#for group 5, FC > 1.5
#all

#for CTvsPT group1, FC > 1.25
#[c( 3, 2, 4, 15, 12, 13, 14)]
#for CTvsPT group2, FC > 1.25
#[c( 26, 27, 4, 5, 43, 2, 32)]

#for FC 1kb, p > 5
#[c( 6, 2, 10, 30, 3, 24, 9, 21, 36, 1)]
#for FP 1kb, p > 5
#[c( 26, 11, 10, 7, 31, 3, 2, 12, 18, 4, 37, 33, 30)]
#for PTup 1kb, p > 4
#[c( 2, 7, 19, 5, 8, 20, 23, 25)]
rm(sel, se)

#pfms = getMatrixByID(JASPAR2022, c("MA0107.1", "MA0491.2", "MA0489.2", 
#                                   "MA1508.1", "MA0496.2", "MA1125.1"))


pfms = getMatrixByID(JASPAR2022, c("MA0107.1")) #MA0517.1 - STAT1::STAT2; MA0107.1 - RELA
pwms = toPWM(pfms)
res <- findMotifHits(query = pwms, subject = seq2, min.score = 6.0, 
                     method = "matchPWM")
m <- table(factor(seqnames(res), levels = names(seq2)),
           factor(res$pwmname, levels = name(pwms)))
seq <- seq2[rowSums(m) > 0, ]


pfms = getMatrixByID(JASPAR2022, c("MA0488.1")) 
#c("MA0652.1", "MA0653.1", "MA1418.1", "MA0476.1") OR c("MA0652.1", "MA0653.1", "MA1418.1", "MA0772.1")
#c("MA0838.1", "MA0102.3", "MA0466.1")
#c("MA0476.1", "MA0488.1")

pwms = toPWM(pfms)
res <- findMotifHits(query = pwms, subject = seq, min.score = 6.0, 
                     method = "matchPWM")
m <- table(factor(seqnames(res), levels = names(seq)),
           factor(res$pwmname, levels = name(pwms)))
seq <- seq[rowSums(m) > 0, ]


#comparing with transcriptomics
CT1 <- read.csv("../figure6_RNAseq/20220809_fullset/CT_1_20220809.csv") 
CT2 <- read.csv("../figure6_RNAseq/20220809_fullset/CT_2_20220809.csv")
CL1 <- read.csv("../figure6_RNAseq/20220809_fullset/CL_1_20220809.csv") 
CL2 <- read.csv("../figure6_RNAseq/20220809_fullset/CL_2_20220809.csv")
PT1 <- read.csv("../figure6_RNAseq/20220809_fullset/PT_1_20220809.csv") 
PT2 <- read.csv("../figure6_RNAseq/20220809_fullset/PT_2_20220809.csv")
PL1 <- read.csv("../figure6_RNAseq/20220809_fullset/PL_1_20220809.csv") 
PL2 <- read.csv("../figure6_RNAseq/20220809_fullset/PL_2_20220809.csv")

foo <- as.data.frame(peakAnno)
foo <- foo[grepl("Promoter", foo$annotation), ]
motifgenes <- unique( foo$SYMBOL )
foo1 <- intersect(CT1$x, motifgenes ) 
length( intersect(CT1$x, motifgenes ) )/length( CT1$x)

#PTsyn = 0.444 PTantag = 0.142, PLsyn = 0.384, PLantag = 14.6, CTsyn = 0.364, CTantag = 0.24, CLsyn = 0.265, CLantag = 0.217

#a small plot of syn vs antag
overlapdata = data.frame(
  treatment = c("polyI:C -> TNF", "polyI:C -> TNF", "CpG -> TNF", "CpG -> TNF", 
                "polyI:C -> LPS", "polyI:C -> LPS", "CpG -> LPS", "CpG -> LPS"),
  group = c("Syn", "Antag", "Syn", "Antag", "Syn", "Antag", "Syn", "Antag"),
  fraction = c(.44, .142, .364, .24, .384, .146, .265, .217)
)

ggplot(overlapdata, aes(x = reorder( group, -fraction) , y=fraction, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge2(reverse = TRUE)) +
  scale_fill_manual(values = c("brown2", "dodgerblue2", "firebrick4", "steelblue4"), 
                    labels = c("CpG -> LPS", "polyI:C -> LPS", "CpG -> TNF", "polyI:C -> TNF")) +
  theme_classic() +
  theme(legend.position="right",
        legend.key.size = unit(0.4, 'cm'), axis.text=element_text(size=11), 
        axis.title=element_text(size=14), legend.text = element_text(size=12), 
        legend.title = element_text(size=12)) +
  labs(x = "",
       y = "Fraction",
       fill = "Treatment") 


peakAnno <- annotatePeak(intersect(MCupbedfile, MPupbedfile), tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
foo <- as.data.frame(peakAnno)
foo <- foo[grepl("Promoter", foo$annotation), ]
peakids <- GRanges(seqnames = foo$seqnames, ranges = IRanges(start = foo$start, end = foo$end),
                   strand = foo$strand)
seq <- getSeq(genome, peakids)
names(seq) <- c(1:length(seq))

#Finding synergistic genes with opened peaks
res <- findMotifHits(query = pwms, subject = seq, min.score = 6.0, 
                     method = "matchPWM")
m <- table(factor(seqnames(res), levels = names(seq)),
           factor(res$pwmname, levels = name(pwms)))
seq <- seq[m[,1] > 0 & rowSums(m[, 2:3 > 0]), ]
diffopen <- unique( foo[names(seq), ]$SYMBOL )
diffopen <- intersect(CT1$x, diffopen ) 


#Correlation between conditions
CPMavgtemp = cpm(DGEobj, log = FALSE, prior.count = 3)
colnames(CPMavgtemp) <- c('CL1', 'CL2', 'CT1','CT2', 'MC1', 'MC2','ML1', 'ML2','MM1', 'MM2', 
                      'MP1', 'MP2', 'MT1','MT2', 'PL1', 'PL2', 'PT1', 'PT2') 
CPMavg = matrix( ,nrow = length(CPMavgtemp[,1]), ncol = 9)
CPMavg[, 1] = rowMeans(CPMavgtemp[, c("CL1", "CL2")])
CPMavg[, 2] = rowMeans(CPMavgtemp[, c("CT1", "CT2")])
CPMavg[, 3] = rowMeans(CPMavgtemp[, c("MC1", "MC2")])
CPMavg[, 4] = rowMeans(CPMavgtemp[, c("ML1", "ML2")])
CPMavg[, 5] = rowMeans(CPMavgtemp[, c("MM1", "MM2")])
CPMavg[, 6] = rowMeans(CPMavgtemp[, c("MP1", "MP2")])
CPMavg[, 7] = rowMeans(CPMavgtemp[, c("MT1", "MT2")])
CPMavg[, 8] = rowMeans(CPMavgtemp[, c("PL1", "PL2")])
CPMavg[, 9] = rowMeans(CPMavgtemp[, c("PT1", "PT2")])
colnames(CPMavg) = c("CL", "CT", "MC", "ML", "MM", "MP", "MT", "PL", "PT")
CPMavg = as.data.frame(CPMavg)

CPMavg_log = log2(CPMavg)
write.csv(CPMavg_log, "logCPMavg.csv")


CPMfcavg = CPMavg/CPMavg$MM
CPMfcavg = log2(CPMfcavg)
CPMfcavg = as.data.frame(CPMfcavg)

ggplot(CPMfcavg[dt[,3] > 0, ], aes(x = PT, y = MT)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) PT",
       y = "log2(FC) MT") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,3] > 0, ]$PT, CPMfcavg[dt[,3] > 0, ]$MT, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,3] > 0, ], aes(x = PT, y = MP)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) PT",
       y = "log2(FC) MP") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,3] > 0, ]$PT, CPMfcavg[dt[,3] > 0, ]$MP, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,4] > 0, ], aes(x = CT, y = MT)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) CT",
       y = "log2(FC) MT") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,4] > 0, ]$CT, CPMfcavg[dt[,4] > 0, ]$MT, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,4] > 0, ], aes(x = CT, y = MC)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) CT",
       y = "log2(FC) MC") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,4] > 0, ]$CT, CPMfcavg[dt[,4] > 0, ]$MC, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,5] > 0, ], aes(x = PL, y = ML)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) PL",
       y = "log2(FC) ML") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,5] > 0, ]$PL, CPMfcavg[dt[,5] > 0, ]$ML, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,5] > 0, ], aes(x = PL, y = MP)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) PL",
       y = "log2(FC) MP") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,5] > 0, ]$PL, CPMfcavg[dt[,5] > 0, ]$MP, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,6] > 0, ], aes(x = CL, y = ML)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) CL",
       y = "log2(FC) ML") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,6] > 0, ]$CL, CPMfcavg[dt[,6] > 0, ]$ML, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,6] > 0, ], aes(x = CL, y = MC)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) CL",
       y = "log2(FC) MC") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,6] > 0, ]$CL, CPMfcavg[dt[,6] > 0, ]$MC, 
         method = "spearman", exact=FALSE )

ggplot(CPMfcavg[dt[,6] > 0, ], aes(x = ML, y = MC)) +
  geom_point(alpha = 0.2, stroke=NA) +
  theme_classic() +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=14)) +
  labs(x = "log2(FC) CL",
       y = "log2(FC) ML") + 
  xlim(-1, 8) +
  ylim(-1, 8)
cor.test(CPMfcavg[dt[,6] > 0, ]$CL, CPMfcavg[dt[,6] > 0, ]$ML, 
         method = "spearman", exact=FALSE )
