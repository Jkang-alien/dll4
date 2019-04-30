library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
studies <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
cs <- getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[34,1]

caselist <- getCaseLists(mycgds,mycancerstudy)

mycaselist = getCaseLists(mycgds,mycancerstudy)[7,1]

# Get available genetic profiles
geneticprofile = getGeneticProfiles(mycgds,mycancerstudy)

mirna = getGeneticProfiles(mycgds,mycancerstudy)[4,1]
mrna = getGeneticProfiles(mycgds,mycancerstudy)[11,1]

# Get data slices for a specified list of genes, genetic profile and case list
mir30a <- getProfileData(mycgds,'MIR-30A/30A*', mirna, mycaselist)
dll4 <- getProfileData(mycgds,'DLL4',mrna,mycaselist)

ID_mRNA <- gsub('\\.', '-', rownames(dll4))
ID_mRNA <- gsub('-0.', '', ID_mRNA)

dll4 <- data.frame(ID = ID_mRNA, dll4)
miR30a <- data.frame(ID = rownames(data_tumor), 
                     miR_30a = data_tumor[,colnames(data_tumor) == 'hsa-mir-30a'])

data <- merge(dll4, miR30a, by = 'ID')


clinical <- getClinicalData(mycgds,mycaselist)
summary(clinical)
head(clinical)
head(df_mir34a)

clinical$ID <- gsub('-0[0-9]{1}', '', gsub('\\.', '-', rownames(clinical)))

library(dplyr)
clinical <- clinical %>% 
  select(ID, AGE, AJCC_NODES_PATHOLOGIC_PN, AJCC_PATHOLOGIC_TUMOR_STAGE, 
         DFS_MONTHS, DFS_STATUS, OS_MONTHS, OS_STATUS)

data_clinical_mir_34a <- merge(df_mir34a, clinical, by = 'ID')

data_clinical_mir_34a <- data_clinical_mir_34a %>% 
  mutate (mir_34a_median = hsa.mir.34a > median(hsa.mir.34a) )
data_clinical_mir_34a$DFS_STATUS <- factor(data_clinical_mir_34a$DFS_STATUS, exclude = '')
data_clinical_mir_34a <- data_clinical_mir_34a %>%
  mutate(DFS_STATUS_no = DFS_STATUS == 'Recurred/Progressed')
library(survival)
surv_object <- Surv(time = data_clinical_mir_34a$DFS_MONTHS, 
                    event = data_clinical_mir_34a$DFS_STATUS_no)
fit <- survfit(surv_object ~ mir_34a_median, data = data_clinical_mir_34a)
summary(fit)

library(ggplot2)
library(survminer)
ggsurvplot(fit, pval = TRUE)

write.csv(data_clinical_mir_34a, 'data_TCGA.csv', quote = FALSE, row.names = FALSE)
