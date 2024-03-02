#setwd("./ChEA3s")
set.seed(3)


library(ggplot2)
library(tidyverse)

logCPM = read.csv("../logCPM.csv")
logCPM = logCPM[!is.na(logCPM$X), ]
rownames(logCPM) = logCPM$X
logCPM = logCPM[, 2:34]
write.csv(logCPM, "Macrophagememory_RNAlogCPM.csv")

### PLOT FOR PT vs CT GROUP 1 GSEA ### 

PT_1 <- read.table("PT_1_CHEA3_20220809.tsv", sep = '\t', header = TRUE)
CT_1 <- read.table("CT_1_CHEA3_20220809.tsv", sep = '\t', header = TRUE)

PT_1 <- PT_1[ ,c(5, 8)]
PT_1 <- PT_1[which( tolower( PT_1$TF) %in% tolower( rownames(logCPM) ) ), ]
PT_1tops <- PT_1[c( 1:10), ]

CT_1 <- CT_1[ ,c(5, 8)]
CT_1 <- CT_1[which( tolower( CT_1$TF) %in% tolower( rownames(logCPM) ) ), ]
CT_1tops <- CT_1[c( 1:10), ]
CT_1 <-  CT_1[ match(PT_1[,1],CT_1[,1]) ,]

datapoints <- cbind( PT_1, CT_1$FET.p.value)


both = intersect(PT_1tops[,1], CT_1tops[,1])

datapoints$color <- NA
rownames(datapoints) <- datapoints[,1]
colnames(datapoints) <- c("TF", "PT_pval", "CT_pval", "color")
for (aa in rownames( datapoints) )
  
{ 
  if (aa %in% both) { datapoints[aa, 4] <- 0}
  else if (aa %in% CT_1tops[,1]) { datapoints[aa, 4] <- 2}
  else if (aa %in% PT_1tops[,1]) { datapoints[aa, 4] <- 1}
  else { datapoints[aa, 4] <- 3}
  
}


ggplot() + 
  geom_point( aes(x=-log10( datapoints$CT_pval), y=-log10( datapoints$PT_pval), color = as.factor( datapoints$color ) ), size = 3,  ) +
  scale_y_continuous(limits = c(0, 25), oob = scales::squish) +
  scale_x_continuous(limits = c(0, 25), oob = scales::squish) + 
  geom_text(data=subset(datapoints, PT_pval < 4e-07 | CT_pval < 0.008), aes(x=-log10( CT_pval)+1.5, y=-log10( PT_pval)+1, label =  TF ), size = 5) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.border = element_blank(),
        text = element_text(size = 20), 
        legend.background = element_blank(), 
        legend.key = element_rect(fill = NA, color = NA)) +
  xlab("-log(p-Value) CT") + 
  ylab("-log(p-Value) PT") +
  labs(color = "Top enriched TFs") +
  scale_color_manual(labels = c("Both", "PT only", "CT only", "Neither"), values=c("#cc00ff","#66ccff", "#ff0000", "#000000"))


### PLOT FOR PT vs CT GROUP 2 GSEA ###

PT_2 <- read.table("PT_2_CHEA3_20220809.tsv", sep = '\t', header = TRUE)
CT_2 <- read.table("CT_2_CHEA3_20220809.tsv", sep = '\t', header = TRUE)

PT_2 <- PT_2[ ,c(5, 8)]
PT_2 <- PT_2[which( tolower( PT_2$TF) %in% tolower( rownames(logCPM) ) ), ]
PT_2tops <- PT_2[c( 1:10), ]

CT_2 <- CT_2[  ,c(5, 8)]
CT_2 <- CT_2[which( tolower( CT_2$TF) %in% tolower( rownames(logCPM) ) ), ]
CT_2tops <- CT_2[c( 1:10), ]
CT_2 <-  CT_2[ match(PT_2[,1],CT_2[,1]) ,]

both = intersect(PT_2tops[,1], CT_2tops[,1])

datapoints2 <- cbind( PT_2, CT_2$FET.p.value)
datapoints2$color <- NA
rownames(datapoints2) <- datapoints2[,1]
colnames(datapoints2) <- c("TF", "PT_pval", "CT_pval", "color")
for (aa in rownames( datapoints2) )
  
{ 
  if (aa %in% both) { datapoints2[aa, 4] <- 0}
  else if (aa %in% CT_2tops[,1]) { datapoints2[aa, 4] <- 2}
  else if (aa %in% PT_2tops[,1]) { datapoints2[aa, 4] <- 1}
  else { datapoints2[aa, 4] <- 3}
  
}


ggplot() + 
  geom_point( aes(x=-log10( datapoints2$CT_pval), y=-log10( datapoints2$PT_pval), color = as.factor( datapoints2$color ) ), size = 3,  ) +
  scale_y_continuous(limits = c(0, 10), oob = scales::squish) +
  scale_x_continuous(limits = c(0, 10), oob = scales::squish) + 
  geom_text(data=subset(datapoints2, PT_pval < .00135 | CT_pval < 8.762e-4), aes(x=-log10( CT_pval)+.5, y=-log10( PT_pval)+.5, label =  TF ), size = 5) + 
  theme(panel.background = element_rect(fill = "white", size = 4, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.border = element_blank(),
        text = element_text(size = 20), 
        legend.background = element_blank(), 
        legend.key = element_rect(fill = NA, color = NA)) +
  xlab("log(p-Value) CT") + 
  ylab("log(p-Value) PT") +
  labs(color = "Top enriched TFs") +
  scale_color_manual(labels = c("Both", "PT only", "CT only", "Neither"), values=c("#cc00ff","#66ccff", "#ff0000", "#000000"))

### PLOT FOR PL vs CL GROUP 1 GSEA ### 

PL_1 <- read.table("PL_1_CHEA3_20220809.tsv", sep = '\t', header = TRUE)
CL_1 <- read.table("CL_1_CHEA3_20220809.tsv", sep = '\t', header = TRUE)

PL_1 <- PL_1[ ,c(5, 8)]
PL_1 <- PL_1[which( tolower( PL_1$TF) %in% tolower( rownames(logCPM) ) ), ]
PL_1tops <- PL_1[c( 1:10), ]

CL_1 <- CL_1[ ,c(5, 8)]
CL_1 <- CL_1[which( tolower( CL_1$TF) %in% tolower( rownames(logCPM) ) ), ]
CL_1tops <- CL_1[c( 1:10), ]
CL_1 <-  CL_1[ match(PL_1[,1],CL_1[,1]) ,]

datapoints <- cbind( PL_1, CL_1$FET.p.value)


both = intersect(PL_1tops[,1], CL_1tops[,1])

datapoints$color <- NA
rownames(datapoints) <- datapoints[,1]
colnames(datapoints) <- c("TF", "PL_pval", "CL_pval", "color")
for (aa in rownames( datapoints) )
  
{ 
  if (aa %in% both) { datapoints[aa, 4] <- 0}
  else if (aa %in% CL_1tops[,1]) { datapoints[aa, 4] <- 2}
  else if (aa %in% PL_1tops[,1]) { datapoints[aa, 4] <- 1}
  else { datapoints[aa, 4] <- 3}
  
}


ggplot() + 
  geom_point( aes(x=-log10( datapoints$CL_pval), y=-log10( datapoints$PL_pval), color = as.factor( datapoints$color ) ), size = 3,  ) +
  scale_y_continuous(limits = c(0, 45), oob = scales::squish) +
  scale_x_continuous(limits = c(0, 45), oob = scales::squish) + 
  geom_text(data=subset(datapoints, PL_pval < 1e-08 | CL_pval < 5e-08), aes(x=-log10( CL_pval)+1.5, y=-log10( PL_pval)+1, label =  TF ), size = 5) + 
  theme(panel.background = element_rect(fill = "white", size = 4, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.border = element_blank(),
        text = element_text(size = 20), 
        legend.background = element_blank(), 
        legend.key = element_rect(fill = NA, color = NA)) +
  xlab("-log(p-Value) CL") + 
  ylab("-log(p-Value) PL") +
  labs(color = "Top enriched TFs") +
  scale_color_manual(labels = c("Both", "PL only", "CL only", "Neither"), values=c("#cc00ff","#66ccff", "#ff0000", "#000000"))


### PLOT FOR PL vs CL GROUP 2 GSEA ### 

PL_2 <- read.table("PL_2_CHEA3_20220809.tsv", sep = '\t', header = TRUE)
CL_2 <- read.table("CL_2_CHEA3_20220809.tsv", sep = '\t', header = TRUE)

PL_2 <- PL_2[ ,c(5, 8)]
PL_2 <- PL_2[which( tolower( PL_2$TF) %in% tolower( rownames(logCPM) ) ), ]
PL_2tops <- PL_2[c( 1:10), ]

CL_2 <- CL_2[ ,c(5, 8)]
CL_2 <- CL_2[which( tolower( CL_2$TF) %in% tolower( rownames(logCPM) ) ), ]
CL_2tops <- CL_2[c( 1:10), ]
CL_2 <-  CL_2[ match(PL_2[,1],CL_2[,1]) ,]

datapoints <- cbind( PL_2, CL_2$FET.p.value)


both = intersect(PL_2tops[,1], CL_2tops[,1])

datapoints$color <- NA
rownames(datapoints) <- datapoints[,1]
colnames(datapoints) <- c("TF", "PL_pval", "CL_pval", "color")
for (aa in rownames( datapoints) )
  
{ 
  if (aa %in% both) { datapoints[aa, 4] <- 0}
  else if (aa %in% CL_2tops[,1]) { datapoints[aa, 4] <- 2}
  else if (aa %in% PL_2tops[,1]) { datapoints[aa, 4] <- 1}
  else { datapoints[aa, 4] <- 3}
  
}


ggplot() + 
  geom_point( aes(x=-log10( datapoints$CL_pval), y=-log10( datapoints$PL_pval), color = as.factor( datapoints$color ) ), size = 3,  ) +
  scale_y_continuous(limits = c(0, 15), oob = scales::squish) +
  scale_x_continuous(limits = c(0, 15), oob = scales::squish) + 
  geom_text(data=subset(datapoints, PL_pval < 2e-04 | CL_pval < 1e-06), aes(x=-log10( CL_pval)+1.5, y=-log10( PL_pval)+1, label =  TF ), size = 5) + 
  theme(panel.background = element_rect(fill = "white", size = 4, colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.border = element_blank(),
        text = element_text(size = 20), 
        legend.background = element_blank(), 
        legend.key = element_rect(fill = NA, color = NA)) +
  xlab("-log(p-Value) CL") + 
  ylab("-log(p-Value) PL") +
  labs(color = "Top enriched TFs") +
  scale_color_manual(labels = c("Both", "PL only", "CL only", "Neither"), values=c("#cc00ff","#66ccff", "#ff0000", "#000000"))
