library(miRNATCGA)
library(RTCGAToolbox)
library(cgdsr)
library(Cairo)
library(ssGSEA4TCGA)
library(ConsensusClusterPlus)
library(NMF)


###https://www.biostars.org/p/215175/

manifest= "gdc_manifest.2017-03-29T01-21-07.652710.tsv" #Manifest name 
x=read.table(manifest,header = T)

id= toString(sprintf('"%s"', x$id))

Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '


Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
"size":"1500"} '

Sentence= paste(Part1,id,Part2, collapse=" ") #This creates the search sentence for the command line



write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)

###### Go to Linux server and follow the instruction (https://www.biostars.org/p/215175/)

file_metadata <- read.delim('File_metadata.txt')

file_names <- list.dirs(path = "./miRNA_file/", full.names = FALSE, recursive = FALSE)

data_miRNA <- c()

for (i in file_names){
  path <-  paste('./miRNA_file/', i, '/mirnas.quantification.txt', sep = '')
  a <- read.delim(path)[,3]
  data_miRNA <- rbind(data_miRNA, a) 
  
}

miRNA_ID <- as.character(read.delim(paste('./miRNA_file/', file_names[1], '/mirnas.quantification.txt', sep = ''))[,1])
colnames(data_miRNA) <- miRNA_ID
data_miRNA <-data.frame(file_id = file_names, data_miRNA)



data <- merge(file_metadata, data_miRNA, by = 'file_id')


colnames(data)[1:15]
data_pt <- subset(data, cases_0_samples_0_sample_type == 'Primary Tumor')
data_m <- subset(data, cases_0_samples_0_sample_type == 'Metastatic')
data_n <- subset(data, cases_0_samples_0_sample_type == 'Solid Tissue Normal')
colnames(data_pt)[1:15]
data_pt <- data_pt[,c(2, 14:1894)]
colnames(data_pt)[1] <- 'ID'

#http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000088
rdate <- getFirehoseRunningDates(last = NULL)
dset <- 'BRCA'

readData <- getFirehoseData (dataset=dset, runDate=rdate[1], RNAseq2_Gene_Norm=TRUE)
mRNA <- getData(readData, 'RNASeq2GeneNorm')

readData = getFirehoseData (dataset=dset, runDate=rdate[1], Clinic = TRUE)
clin = getData(readData, "Clinical")

mRNA <- t(mRNA)

ID <- gsub('-...-...-....-..', '', rownames(mRNA))
index_normal <- grepl('TCGA-..-....-1..-...-....-..', rownames(mRNA))
index_tumor <- grepl('TCGA-..-....-0..-...-....-..', rownames(mRNA))
mRNA_tumor <- mRNA [index_tumor, ] 
mRNA_normal <- mRNA [index_normal,]

rownames (mRNA_tumor) <- ID[index_tumor]
rownames(mRNA_tumor) 

rownames (mRNA_normal) <- ID[index_normal]
rownames(mRNA_normal) 

sum(duplicated(rownames(mRNA_tumor)))

mRNA_tumor <- mRNA_tumor[duplicated(rownames(mRNA_tumor)) == FALSE,]

sum(duplicated(rownames(mRNA_normal)))

mRNA_normal <- mRNA_normal[duplicated(rownames(mRNA_normal)) == FALSE,]

##############################################################

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
studies <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
cs <- getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[22,1]

caselist <- getCaseLists(mycgds,mycancerstudy)

mycaselist = getCaseLists(mycgds,mycancerstudy)[21,1]

ID_luminalA <- gsub('-01', '', unlist(strsplit(caselist[6,5], ' ')))
ID_basal <- gsub('-01', '', unlist(strsplit(caselist[7,5], ' ')))
ID_her2 <- gsub('-01', '', unlist(strsplit(caselist[8,5], ' ')))
ID_luminalB <- gsub('-01', '', unlist(strsplit(caselist[9,5], ' ')))
ID_lobular <- gsub('-01', '', unlist(strsplit(caselist[10,5], ' ')))
ID_mixed <- gsub('-01', '', unlist(strsplit(caselist[12,5], ' ')))
ID_others <- gsub('-01', '', unlist(strsplit(caselist[13,5], ' ')))

subtype <- rep(NA, dim(data_pt)[1])
subtype [data_pt$ID %in% ID_luminalA] <- 'LuminalA'
subtype [data_pt$ID %in% ID_basal] <- 'Basal-like'
subtype [data_pt$ID %in% ID_her2] <- 'HER2'
subtype [data_pt$ID %in% ID_luminalB] <- 'LuminalB'
subtype [data_pt$ID %in% ID_lobular] <- 'Lobular'
subtype [data_pt$ID %in% ID_others] <- 'Others'

data_pt$subtype <- factor(subtype)

mRNA <- data.frame(ID = rownames(mRNA_tumor), mRNA_tumor)

colnames(data_pt)[grep('\\.30a', colnames(data_pt))]
grep('\\.30a', colnames(data_pt))
colnames(data_pt)[grep('\\.214', colnames(data_pt))]
grep('\\.214', colnames(data_pt))
df_mir30a_214 <- data_pt[,c(1,297,354,1883)]

grep('\\.34a', colnames(data_pt))
df_mir34a <- data_pt[,c(1, 500)]

data <- merge(df_mir30a_214, mRNA, by = 'ID')
summary(data$subtype)

data_la <- subset(data, subtype == 'LuminalA')
data_lb <- subset(data, subtype == 'LuminalB')
data_h <- subset(data, subtype == 'HER2')
data_b <- subset(data, subtype == 'Basal-like')
############## Basal-like ###############################


p_value <- c()

for (i in 5:dim(data_b)[2]){
  a <- cor.test(data_b[,i], data_b$hsa.mir.214, method = 'spearman')
  b <- a$p.value
  p_value <- rbind(p_value, b)
}


#tapply(data[,3], data$group_2, quantile)$A[4]

cor_value <- c()
for (i in 5:dim(data_b)[2]){
  a <- cor(data_b[,i], data_b$hsa.mir.214, method = 'spearman')
  cor_value <- rbind(cor_value, a)
}

summary(p_value)
summary(cor_value)

genes <- colnames(data)[5:dim(data)[2]]
gene_mir214_cor <- genes[abs(cor_value) > 0.4 & is.na(cor_value) == FALSE]
gene_mir214_inverse_cor <- genes[cor_value < -0.4 & is.na(cor_value) == FALSE]

cor_value <- c()
for (i in 5:dim(data_la)[2]){
  a <- cor(data_la[,i], data_la$hsa.mir.214, method = 'spearman')
  cor_value <- rbind(cor_value, a)
}

summary(p_value)
summary(cor_value)

genes <- colnames(data)[5:dim(data)[2]]
gene_mir214_cor_la <- genes[abs(cor_value) > 0.4 & is.na(cor_value) == FALSE]
gene_mir214_inverse_cor <- genes[cor_value < -0.4 & is.na(cor_value) == FALSE]

cor_value <- c()
for (i in 5:dim(data_h)[2]){
  a <- cor(data_h[,i], data_h$hsa.mir.214, method = 'spearman')
  cor_value <- rbind(cor_value, a)
}

summary(p_value)
summary(cor_value)

genes <- colnames(data)[5:dim(data)[2]]
gene_mir214_cor_h <- genes[abs(cor_value) > 0.4 & is.na(cor_value) == FALSE]
gene_mir214_inverse_cor <- genes[cor_value < -0.4 & is.na(cor_value) == FALSE]


write.table(gene_mir214_inverse_cor,"genes_mir214_inverse_cor.txt",quote=F,col.names=F,row.names=F)

gene_mir214_positive_cor <- genes[cor_value > 0.4 & is.na(cor_value) == FALSE]

write.table(gene_mir214_positive_cor,"genes_mir214_positive_cor.txt",quote=F,col.names=F,row.names=F)

cairo_pdf(filename = 'mir214_cor.pdf',
          width = 7, height = 7, pointsize = 12,
          family = "sans", bg = "white")
hist(cor_value)
dev.off()

cairo_pdf(filename = 'mir214_cor_gene.pdf',
          width = 14, height = 10, pointsize = 20,
          family = "sans", bg = "white")
par(mfrow = c(2,3))

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$CXCL12, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 CXCL12')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$PDGFRB, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 PDGFRB')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$MMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 MMP2')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$TIMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 TIMP2')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$RECK, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 RECK')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$ERCC3, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 ERCC3')

dev.off()

cairo_pdf(filename = 'mir214_cor_negative_gene.pdf',
          width = 14, height = 10, pointsize = 20,
          family = "sans", bg = "white")
par(mfrow = c(2,3))



plot(log(data_b$hsa.mir.214, base = 10), log(data_b$PDGFRB, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 PDGFRB')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$MMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 MMP2')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$TIMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 TIMP2')

plot(log(data_b$hsa.mir.214, base = 10), log(data_b$RECK, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 RECK')
dev.off()

cairo_pdf(filename = 'mir214_PTEN.pdf',
          width = 7, height = 7, pointsize = 12,
          family = "sans", bg = "white")
plot(log(data_b$hsa.mir.214, base = 10), log(data_b$PTEN, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 PTEN')
text (2.5, 2.5, 'cor = 0.07')
dev.off()

cor_value [genes == 'PTEN']

############## Luminal A ###############################################
p_value <- c()

for (i in 5:dim(data_la)[2]){
  a <- cor.test(data_la[,i], data_la$hsa.mir.214, method = 'spearman')
  b <- a$p.value
  p_value <- rbind(p_value, b)
}


#tapply(data[,3], data$group_2, quantile)$A[4]

cor_value <- c()
for (i in 5:dim(data_la)[2]){
  a <- cor(data_la[,i], data_la$hsa.mir.214, method = 'spearman')
  cor_value <- rbind(cor_value, a)
}

summary(p_value)
summary(cor_value)

genes <- colnames(data)[5:dim(data)[2]]

gene_mir214_inverse_cor <- genes[cor_value < -0.2 & is.na(cor_value) == FALSE]

write.table(gene_mir214_inverse_cor,"genes_mir214_inverse_cor_la.txt",quote=F,col.names=F,row.names=F)

gene_mir214_positive_cor <- genes[cor_value > 0.4 & is.na(cor_value) == FALSE]

write.table(gene_mir214_positive_cor,"genes_mir214_positive_cor_la.txt",quote=F,col.names=F,row.names=F)

cairo_pdf(filename = 'mir214_cor_la.pdf',
          width = 7, height = 7, pointsize = 12,
          family = "sans", bg = "white")
hist(cor_value)
dev.off()

cairo_pdf(filename = 'mir214_cor_gene_la.pdf',
          width = 14, height = 10, pointsize = 20,
          family = "sans", bg = "white")
par(mfrow = c(2,3))

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$CXCL12, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 CXCL12')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$PDGFRB, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 PDGFRB')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$MMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 MMP2')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$TIMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 TIMP2')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$RECK, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 RECK')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$ERCC3, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 ERCC3')

dev.off()

cairo_pdf(filename = 'mir214_cor_negative_gene_la.pdf',
          width = 14, height = 10, pointsize = 20,
          family = "sans", bg = "white")
par(mfrow = c(2,3))



plot(log(data_la$hsa.mir.214, base = 10), log(data_la$PDGFRB, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 PDGFRB')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$MMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 MMP2')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$TIMP2, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 TIMP2')

plot(log(data_la$hsa.mir.214, base = 10), log(data_la$RECK, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 RECK')
dev.off()

cairo_pdf(filename = 'mir214_PTEN_la.pdf',
          width = 7, height = 7, pointsize = 12,
          family = "sans", bg = "white")
plot(log(data_la$hsa.mir.214, base = 10), log(data_la$PTEN, base = 10),
     xlab = 'log10 Mir-214',
     ylab = 'log10 PTEN')
text (2.5, 2.5, 'cor = 0.07')
dev.off()


########################################################
########### Mir-30a, DLL3 ##############################

CairoPNG(filename = "mir30_DLL3_cor.png", width = 1200, height = 1200,
         pointsize = 24, bg = "white")

par(mfrow = c(2,2))

cor(log(data_b$hsa.mir.30a, base = 10) , log(data_b$DLL4, base = 10))
cor.test(log(data_b$hsa.mir.30a, base = 10) , log(data_b$DLL4, base = 10))
plot(log(data_b$hsa.mir.30a, base = 10) , log(data_b$DLL4, base = 10),
     xlab = 'Mir-30a expression (log-scale)',
     ylab = 'DLL4 expression (log-scale)',
     main = 'Basal-like')
legend('bottomleft', 'cor = -0.24, P = 0.018', bty = 'n')
abline(lm(log(data_b$DLL4, base = 10) ~ log(data_b$hsa.mir.30a, base = 10)))


cor(log(data_la$hsa.mir.30a, base = 10) , log(data_la$DLL4, base = 10))
cor.test(log(data_la$hsa.mir.30a, base = 10) , log(data_la$DLL4, base = 10))
plot(log(data_la$hsa.mir.30a, base = 10) , log(data_la$DLL4, base = 10),
     xlab = 'Mir-30a expression (log-scale)',
     ylab = 'DLL4 expression (log-scale)',
     main = 'Luminal A')
legend('bottomleft', 'cor = -0.007, P = 0.927', bty = 'n')
abline(lm(log(data_la$DLL4, base = 10) ~ log(data_la$hsa.mir.30a, base = 10)))

cor(log(data_lb$hsa.mir.30a, base = 10) , log(data_lb$DLL4, base = 10))
cor.test(log(data_lb$hsa.mir.30a, base = 10) , log(data_lb$DLL4, base = 10))
plot(log(data_lb$hsa.mir.30a, base = 10) , log(data_lb$DLL4, base = 10),
     xlab = 'Mir-30a expression (log-scale)',
     ylab = 'DLL4 expression (log-scale)',
     main = 'Luminal B')
legend('bottomleft', 'cor = -0.116, P = 0.235', bty = 'n')
abline(lm(log(data_lb$DLL4, base = 10) ~ log(data_lb$hsa.mir.30a, base = 10)))

cor(log(data_h$hsa.mir.30a, base = 10) , log(data_h$DLL4, base = 10))
cor.test(log(data_h$hsa.mir.30a, base = 10) , log(data_h$DLL4, base = 10))
plot(log(data_h$hsa.mir.30a, base = 10) , log(data_h$DLL4, base = 10),
     xlab = 'Mir-30a expression (log-scale)',
     ylab = 'DLL4 expression (log-scale)',
     main = 'HER2')
legend('bottomleft', 'cor = -0.004, P = 0.978', bty = 'n')
abline(lm(log(data_h$DLL4, base = 10) ~ log(data_h$hsa.mir.30a, base = 10)))

dev.off()

cor(log(data_la$hsa.mir.30a, base = 10) , log(data_la$DLL4, base = 10))
plot(log(data_la$hsa.mir.30a, base = 10) , log(data_la$DLL4, base = 10))


cor_value [genes == 'PTEN']

loc <- matrix( c(cor_value[abs(cor_value) > 0.5 & p_value < 0.05], 
                 -log(p_value, base =10)[abs(cor_value) > 0.5 & p_value <0.05]), 
               nrow = 2,
               byrow = TRUE)

colnames(loc)
svg(file = 'volcano_214.svg',
    width =7.5, height = 7.5, pointsize = 16)                
vp <- plot(cor_value, -log(p_value, base = 10),
           xlab = 'Difference of mRNA expression',
           ylab = '-Log10(p-value)')
text(loc[1,], loc[2,], colnames(loc), adj = c(0,0.5),
     cex =0.5)
#arrows(0.3, 4, 0.4, 3.25, angle = 30, length = 0.05)
dev.off()



cairo_pdf(filename = 'DLL4_miR30a_expression.pdf',
          width = 7, height = 7, pointsize = 12,
          family = "sans", bg = "white")
par (mfrow = c(2,1))
boxplot(log(data$DLL4, base = 10) ~ data$subtype,
        ylab = 'DLL4 expression (RSEM)')
boxplot(log(data$miR_30a, base = 10) ~ data$subtype,
        ylab = 'miR30a expression')
dev.off()

data_luminalA <- data[data$ID %in% ID_luminalA,]
data_basal <- data[data$ID %in% ID_basal,]
data_her2 <- data[data$ID %in% ID_her2,]
data_luminalB <- data[data$ID %in% ID_luminalB,]

par (mfrow=c(2,2))

cairo_pdf(filename = 'DLL4_miR30a.pdf',
          width = 7, height = 7, pointsize = 12,
          family = "sans", bg = "white")

par (mfrow=c(2,2),
     mar=c(4,6,4,4))
plot(log(data_luminalA$miR_30a, base = 10), log(data_luminalA$DLL4, base = 10),
     xlab='miR30a expression',
     ylab='DLL4 expression \n(RSEM log-scale)',
     main = 'Luminal A')
legend('bottomleft', 'cor = 0.003 P = 1', bty = 'n')
abline(lm(log(data_luminalA$DLL4, base = 10) ~ log(data_luminalA$miR_30a, base = 10)))

plot(log(data_basal$miR_30a), log(data_basal$DLL4),
     xlab='',
     ylab='',
     main = 'Basal-like')
legend('bottomleft', 'cor = -0.19 P = 0.100', bty = 'n')
abline(lm(log(data_basal$DLL4) ~ log(data_basal$miR_30a)))

plot(log(data_her2$miR_30a), log(data_her2$DLL4),
     xlab='',
     ylab='',
     main = 'HER2')
legend('bottomleft', 'cor = 0.27 P = 0.088', bty = 'n')
abline(lm(log(data_her2$DLL4) ~ log(data_her2$miR_30a)))

plot(log(data_luminalB$miR_30a), log(data_luminalB$DLL4),
     xlab='',
     ylab='',
     main = 'Luminal B')
legend('bottomleft', 'cor = 0.11 P = 0.330', bty = 'n')
abline(lm(log(data_luminalB$DLL4) ~ log(data_luminalB$miR_30a)))

dev.off()

sink(file = 'statistics.txt')

summary(data$subtype)
cor.test(log(data_luminalA$DLL4), log(data_luminalA$miR_30a))
cor.test(log(data_basal$DLL4), log(data_basal$miR_30a))
cor.test(log(data_her2$DLL4), log(data_her2$miR_30a))
cor.test(log(data_luminalB$DLL4), log(data_luminalB$miR_30a))
sink()
file.show('statistics.txt')

clin <- survivalTCGA(clin)
clin$ID <- gsub('\\.', '-', clin$ID)
data_clinic <- merge(data, clin, by = 'ID')
data_clinic$surv_months [data_clinic$surv_months < 0 ] <- NA
library(survival)
library(rms)

summary(data_clinic)
data_clinic$pathologic_stage <- 
  factor(data_clinic$pathologic_stage)
summary(data_clinic$patholoic_stage)
stage1 <- levels(data_clinic$pathologic_stage)[1:3]
stage2 <- levels(data_clinic$pathologic_stage)[4:6]
stage3 <- levels(data_clinic$pathologic_stage)[7:10]
stage4 <- levels(data_clinic$pathologic_stage)[11]

data_clinic$stage[as.character(data_clinic$pathologic_stage) %in% stage1] <- 'Stage1'
data_clinic$stage[as.character(data_clinic$pathologic_stage) %in% stage2] <- 'Stage2'
data_clinic$stage[as.character(data_clinic$pathologic_stage) %in% stage3] <- 'Stage3'
data_clinic$stage[as.character(data_clinic$pathologic_stage) %in% stage4] <- 'Stage4'
summary(factor(data_clinic$stage))
data_clinic$stage <- factor(data_clinic$stage)
data_clinic$miR_30a_G <- factor(data_clinic$miR_30a > 200000,
                                levels = c(FALSE, TRUE),
                                labels = c('Low miR30A', 'High miR30A'))
data_clinic$log_miR30a <- log(data_clinic$miR_30a, base = 10)
strata = levels(data_clinic$miR_30a_G)

library(Cairo)

CairoPDF(file = "Survival.pdf",  width = 5, height = 5, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
par(
  #mfrow = c(1,2), 
  mar=c(6,4,2,6), mgp = c(2, 1, 0))


fit = npsurv(Surv(surv_months, vital_status == 1)~
               miR_30a_G, data = data_clinic,
             subset = (data_clinic$subtype == 'Basal-like'))

fit

diff = survdiff(Surv(surv_months, vital_status == 1)~
                  miR_30a_G, data = data_clinic, 
                subset = (data_clinic$subtype == 'Basal-like'))
diff

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:2),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.3, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 0.6, strata, lty = c(1:2), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(55, 0.47, 'P-value = 0.404', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()

CairoPDF(file = "Histogram.pdf",  width = 5, height = 5, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
hist(data_laasal$miR_30a, breaks = 100)
dev.off()

CairoPDF(file = "miR30a_stage.pdf",  width = 5, height = 5, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
boxplot(data_clinic$log_miR30a ~ stage, 
        data = data_clinic, 
        subset = (subtype == 'Basal-like'),
        ylab = 'miR30a expression, log-scale')
dev.off()

################### Multivariate ######################

cox <- coxph(Surv(TCGA_M, CURATED_VITAL_STATUS == 'Dead')~
                    CURATED_AGE_AT_TCGA_SPECIMEN +
                    #stage (Stage at diagnosis is not compatible with post accession survival)
                    UV.signature, 
                  data = data)

survfit(cox_TCGA)

sink('result.txt')
summary(cox_TCGA)
sink()

library(rjson)
metadata <- fromJSON(file = 'metadata.cart.2017-03-16T03-48-19.410382.json')
metadata <- do.call(rbind.data.frame,metadata)
length(metadata[[2]])
metadata <- data.frame(matrix(unlist(metadata), length(metadata[[2]]), byrow=T))

df_metadata <- c()
for (i in 1:length(metadata)){
  file_ID <- metadata[[i]]$file_id
  submitter_ID<- metadata[[i]]$associated_entities[[1]]$entity_submitter_id
  a <- c(file_ID, submitter_ID)
  df_metadata <- rbind(df_metadata, a)
}
colnames(df_metadata) <- c('file_ID', 'ID')


df_metadata[grep('^0123*', df_metadata[,1]),]

###https://www.biostars.org/p/215175/

manifest= "gdc_manifest.2017-03-29T01-21-07.652710.tsv" #Manifest name 
x=read.table(manifest,header = T)

id= toString(sprintf('"%s"', x$id))

Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '


Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
    "size":"500"} '

Sentence= paste(Part1,id,Part2, collapse=" ") #This creates the search sentence for the command line



write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
file_metadata <- read.delim('File_metadata.txt')


gene_list <- c('EZH2', 'PTEN', 'MAP2K3', 'MAPK8', 'PLXNB1', 'POU4F2', 'SRGAP1', 'XBP1',
               'TWIST1', 'CTNNB1', 'BCL2L2', 'BCL2L11', 'JAG1', 'FGFR1', 'NRAS', 'UBE2I', 
               'HDGF', 'CDK6', 'ALPK2', 'PAPPA', 'LZTS1', 'CPEB4', 'TNFSF9', 'TFAP2A', 
               'ING4', 'LTF', 'TP53', 'MEF2C', 'STAT6')

gene_list [gene_list %in% colnames(data_b)]
par(new=F)
CairoPDF(file = 'gene_list.pdf',
         width =20, height = 16, pointsize = 10)
par(mfrow = c(5,6))
for (i in gene_list){
  hist(data_b[,i],
       xlab = i)
} 

dev.off()

hist(log(data_b$hsa.mir.214))
boxplot(log(data_b$hsa.mir.214) ~ data_b$TP53>1000)

cor_value_gene_list <- c()
for (i in gene_list){
  a <- cor(log(data_b[,i]), log(data_b$hsa.mir.214), method = 'pearson')
  cor_value_gene_list <- rbind(cor_value_gene_list, a)
}
hist(cor_value_gene_list,
     main = 'Mir214 and selected gene expression',
     xlab = 'Correlation coefficient')
gene_list [abs(cor_value_gene_list) > 0.3 ]

CairoPDF(file = 'gene_list_cor.pdf',
         width =10, height = 6, pointsize = 16)
par(mfrow = c(1,2))
plot(log(data_b$hsa.mir.214),log(data_b$ALPK2))
plot(log(data_b$hsa.mir.214), log(data_b$MEF2C))

dev.off()


gs <- gs_gmt('h.all.v6.0.symbols.gmt')

mRNA <- t(data_b[,5:dim(data_b)[2]])
es <- gsva(mRNA, gs, method = 'ssgsea', rnaseq = TRUE )
 
data_t <- scale(t(es), center = TRUE, scale = TRUE)
colnames(data_t) <- gsub('http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_', '', colnames(data_t))

results_col = ConsensusClusterPlus(data_t,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_col',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

results_row = ConsensusClusterPlus(t(data_t),maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_row',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

ann_col <- data.frame(immune_class = factor(results_col[[4]]$consensusClass))
ann <- data.frame(class = as.factor( results_row[[2]]$consensusClass))
rownames(ann) <- gsub('-', '\\.', gsub('-...-...-....-..','',rownames(ann)))


ann$mir214 <- log(data_b$hsa.mir.214)

CairoPDF(file = 'cluster.pdf',
         width =7.5, height = 7.5, pointsize = 16)
aheatmap(data_t,
         hclustfun=function(d) hclust(dist(d, method = 'euclidean'), method = "ward.D2"),
         annRow = ann,
         annCol = ann_col,
         Colv = results_col[[4]]$consensusTree,
         Rowv = results_row[[2]]$consensusTree,
         labRow = rep('',dim(data_t)[1]))


dev.off()

cor_value_gene_set <- c()
for (i in 1:dim(data_t)[2]){
  a <- cor(data_t[,i], data_b$hsa.mir.214, method = 'spearman')
  cor_value_gene_set <- rbind(cor_value_gene_set, a)
}

hist(cor_value_gene_set)
colnames(data_t)[abs(cor_value_gene_set)>.4]
