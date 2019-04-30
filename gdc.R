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



