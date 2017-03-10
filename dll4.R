library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
studies <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
cs <- getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[23,1]

caselist <- getCaseLists(mycgds,mycancerstudy)

mycaselist = getCaseLists(mycgds,mycancerstudy)[4,1]

# Get available genetic profiles
geneticprofile = getGeneticProfiles(mycgds,mycancerstudy)

mirna = getGeneticProfiles(mycgds,mycancerstudy)[4,1]
mrna = getGeneticProfiles(mycgds,mycancerstudy)[11,1]

# Get data slices for a specified list of genes, genetic profile and case list
dll4 <- getProfileData(mycgds,'DLL4',mrna,mycaselist)

mir_30a <- getProfileData(mycgds,'mir30a',mirna,mycaselist)

hist(uchl1$UCHL1, 
     breaks = 1000, 
     xlim = c(0,1000),
     plot = TRUE
)

log(200, base = 10)
hist(log(uchl1$UCHL1))

##### http://www.gtexportal.org/home/gene/UCHL1 ##########

ID <- gsub('\\.', '-', (gsub('\\.0.', '', rownames(uchl1))))
df_uchl1 <- data.frame(ID = ID, UCHL1 = uchl1$UCHL1, 
                       log_UCHL1 = log(uchl1$UCHL1+1, base = 10),
                       UCHL1_G = factor(log(uchl1$UCHL1, base = 10) > 2.5, 
                                        levels = c(FALSE, TRUE),
                                        labels = c('Low UCHL1 expression',
                                                   'High UCHL1 expression'))
)
colnames(data)
data$ID <- gsub('-0[0-9]{1}', '', data$Name)
data_uchl1 <- merge(df_uchl1, data, by = 'ID', all.x = TRUE)


colnames(data_uchl1)
boxplot(log_UCHL1 ~ UV.signature, data_uchl1)
t.test(log_UCHL1 ~ UV.signature, data_uchl1)

boxplot(log_UCHL1 ~ MUTATIONSUBTYPES, data_uchl1)
summary(aov(log_UCHL1 ~ MUTATIONSUBTYPES, data_uchl1))

boxplot(log_UCHL1 ~ ALL_PRIMARY_VS_METASTATIC, data_uchl1)
boxplot(log_UCHL1 ~ REGIONAL_VS_PRIMARY, data_uchl1)
boxplot(log_UCHL1 ~ CURATED_TCGA_SPECIMEN_SITE, data_uchl1)

a <- aov(log_UCHL1 ~ CURATED_TCGA_SPECIMEN_SITE, data_uchl1)
summary(a)
summary(data_uchl1)

########################################################

library(compareGroups)

compareGroups(UCHL1_G~ 
                CURATED_AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS +
                CURATED_AGE_AT_TCGA_SPECIMEN +
                GENDER+
                CURATED_SITE_OF_PRIMARY_TUMOR_KNOWN_PRIMARY_ONLY +
                REGIONAL_VS_PRIMARY+
                CURATED_BRESLOW +
                CURATED_ULCERATION +
                PIGMENT.SCORE+
                CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE +
                MUTATIONSUBTYPES+
                UV.signature,
              data = data_uchl1) -> table

table(data$Pathologic_stage, data$uv_g)

createTable(table) -> table
export2csv(table, 'table.csv')


##### Overall survival ######### 
library(survival)
library(rms)


data_uchl1$OS_M <- as.numeric(as.character(data_uchl1$CURATED_DAYS_TO_DEATH_OR_LAST_FU))/30.4
fit = npsurv(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
               UCHL1_G, data = data_uchl1)
fit


strata = levels(data_uchl1$UCHL1_G)

library(Cairo)
CairoSVG(file = "OS.svg",  width = 6, height = 6, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
par(mar=c(6,4,2,8), mgp = c(2, 1, 0))
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
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:2), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.445', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()

diff = survdiff(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
                  UCHL1_G, data = data_uchl1)
diff


data_uchl1$stage <- 
  factor(data_uchl1$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE %in% 
           c('Stage III' ,'Stage IV'),
         levels = c(FALSE, TRUE),
         labels = c('Stage 0-II', 'Stage III-IV'))
data_uchl1$stage[is.na(data_uchl1$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE)] <- NA
summary (data_uchl1$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE)

cox <- coxph(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
               stage + 
               CURATED_AGE_AT_TCGA_SPECIMEN+
               UCHL1_G, data = data_uchl1)

summary(cox)

sink('cox_analysis_output.txt')

summary(cox)

sink()