library(miRNATCGA)
library(RTCGAToolbox)
library(cgdsr)
library(Cairo)
library(ssGSEA4TCGA)

#http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000088
rdate <- getFirehoseRunningDates(last = NULL)
dset <- 'BRCA'

readData <- getFirehoseData (dataset=dset, runDate=rdate[1], RNAseq2_Gene_Norm=TRUE)
mRNA <- getData(readData, 'RNASeq2GeneNorm')

readData = getFirehoseData (dataset=dset, runDate=rdate[1], Clinic = FALSE, miRNASeq_Gene=TRUE)
miRNA <- getData(readData, 'miRNASeqGene')

readData = getFirehoseData (dataset=dset, runDate=rdate[1], Clinic = TRUE)
clin = getData(readData, "Clinical")

mRNA <- t(mRNA)
miRNA <- t(miRNA)

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

############ miRNA ###############################
ID <- gsub('-...-...-....-..', '', rownames(miRNA))
index_normal <- grepl('TCGA-..-....-1..-...-....-..', rownames(miRNA))
index_tumor <- grepl('TCGA-..-....-0..-...-....-..', rownames(miRNA))
miRNA_tumor <- miRNA [index_tumor, ] 
miRNA_normal <- miRNA [index_normal,]

rownames (miRNA_tumor) <- ID[index_tumor]
rownames(miRNA_tumor) 

rownames (miRNA_normal) <- ID[index_normal]
rownames(miRNA_normal) 

sum(duplicated(rownames(miRNA_tumor)))

miRNA_tumor <- miRNA_tumor[duplicated(rownames(miRNA_tumor)) == FALSE,]

sum(duplicated(rownames(miRNA_normal)))

miRNA_normal <- miRNA_normal[duplicated(rownames(miRNA_normal)) == FALSE,]

################################################################################

dll4 <- data.frame(ID = rownames(mRNA_tumor), 
                   DLL4 = mRNA_tumor[,colnames(mRNA_tumor) == 'DLL4'])
miR30a <- data.frame(ID = rownames(miRNA_tumor), 
                     miR_30a = miRNA_tumor[,colnames(miRNA_tumor) == 'hsa-mir-30a'])
data <- merge(dll4, miR30a, by = 'ID')

##############################################################
plot(log(data$miR_30a), log(data$DLL4))
cor.test(log(data$miR_30a), log(data$DLL4))


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

subtype <- rep(NA, dim(data)[1])
subtype [data$ID %in% ID_luminalA] <- 'LuminalA'
subtype [data$ID %in% ID_basal] <- 'Basal-like'
subtype [data$ID %in% ID_her2] <- 'HER2'
subtype [data$ID %in% ID_luminalB] <- 'LuminalB'
subtype [data$ID %in% ID_lobular] <- 'Lobular'
subtype [data$ID %in% ID_others] <- 'others'

data$subtype <- factor(subtype)

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
hist(data_basal$miR_30a, breaks = 100)
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
