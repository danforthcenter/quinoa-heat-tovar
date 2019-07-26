#GO analysis
setwd("/Users/jtovar/quinoa_heat/")

#Load quinoa genome annotation, downloaded from https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Cquinoa_er
Cquinoa_392_v1.0.annotation_info <- read.delim("~/Documents/kallisto_test/Cquinoa/v1.0/annotation/Cquinoa_392_v1.0.annotation_info.txt", comment.char="#")
Cquinoa_392_v1.0.annotation_info_sub <- Cquinoa_392_v1.0.annotation_info[,c("transcriptName", "Panther")]
colnames(Cquinoa_392_v1.0.annotation_info_sub)[1]<- "target_id"

#Panther needs only one panther_id per gene, so delete everything after the comma on the column Panther
Cquinoa_392_v1.0.annotation_info_sub$Panther <- gsub("(.*),.*", "\\1", Cquinoa_392_v1.0.annotation_info_sub$Panther)

#Panther requires that for genes with no Panther ID, we use "NOHIT"
Cquinoa_392_v1.0.annotation_info_sub$Panther <- sub("^$", "NOHIT", Cquinoa_392_v1.0.annotation_info_sub$Panther)

#Load lists of differentially expressed genes by day and treatment
hr_day1_DE_genes.qval_0.05 <- read.delim("~/Documents/kallisto_test/hr_day1_DE_genes.qval_0.05.txt")
hr_day11_DE_genes.qval_0.05 <- read.delim("~/Documents/kallisto_test/hr_day11_DE_genes.qval_0.05.txt")
hrs_day1_DE_genes.qval_0.05 <- read.delim("~/Documents/kallisto_test/hrs_day1_DE_genes.qval_0.05.txt")
hrs_day11_DE_genes.qval_0.05 <- read.delim("~/Documents/kallisto_test/hrs_day11_DE_genes.qval_0.05.txt")
hs_day1_DE_genes.qval_0.05 <- read.delim("~/Documents/kallisto_test/hs_day1_DE_genes.qval_0.05.txt")
hs_day11_DE_genes.qval_0.05 <- read.delim("~/Documents/kallisto_test/hs_day11_DE_genes.qval_0.05.txt")

#Subset only relevant info from DE gene lists
hr_day1_DE_genes.qval_0.05_sub <- hr_day1_DE_genes.qval_0.05[,c("target_id","b")]
hr_day11_DE_genes.qval_0.05_sub <- hr_day11_DE_genes.qval_0.05[,c("target_id","b")]
hrs_day1_DE_genes.qval_0.05_sub <- hrs_day1_DE_genes.qval_0.05[,c("target_id","b")]
hrs_day11_DE_genes.qval_0.05_sub <- hrs_day11_DE_genes.qval_0.05[,c("target_id","b")]
hs_day1_DE_genes.qval_0.05_sub <- hs_day1_DE_genes.qval_0.05[,c("target_id","b")]
hs_day11_DE_genes.qval_0.05_sub <- hs_day11_DE_genes.qval_0.05[,c("target_id","b")]

#Create lists for uploading into http://www.pantherdb.org/
hr_day1_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, hr_day1_DE_genes.qval_0.05_sub, by='target_id')
hr_day11_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, hr_day11_DE_genes.qval_0.05_sub, by='target_id')
hrs_day1_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, hrs_day1_DE_genes.qval_0.05_sub, by='target_id')
hrs_day11_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, hrs_day11_DE_genes.qval_0.05_sub, by='target_id')
hs_day1_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, hs_day1_DE_genes.qval_0.05_sub, by='target_id')
hs_day11_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, hs_day11_DE_genes.qval_0.05_sub, by='target_id')

#Save text files to upload to http://www.pantherdb.org/ These can be used for both statistical enrichment 
#or statistical overrepresentation tests
write.table(hr_day1_DEGs_wpantherid, file='hr_day1_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)
write.table(hr_day11_DEGs_wpantherid, file='hr_day11_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)
write.table(hrs_day1_DEGs_wpantherid, file='hrs_day1_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)
write.table(hrs_day11_DEGs_wpantherid, file='hrs_day11_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)
write.table(hs_day1_DEGs_wpantherid, file='hs_day1_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)
write.table(hs_day11_DEGs_wpantherid, file='hs_day11_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)

#To run a statistical overrepresentation test, we need the genome file. Save it
write.table(Cquinoa_392_v1.0.annotation_info_sub, file='Cquinoa_392_v1.0.panther.annotation.txt', sep="\t",row.names=F, quote=F)

#Obtain the lists of genes differentially expressed in both treatments with heated shoots
#Day 1 of heat treatment
diff_hrs_hs_day1_exclusive <- read.delim("~/Documents/kallisto_test/hrs_hs_day1_exclusive_DE_genes.txt")
hrs_hs_day1_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, diff_hrs_hs_day1_exclusive, by='target_id')
hrs_hs_day1_DEGs_wpantherid1 <- hrs_hs_day1_DEGs_wpantherid[, c("target_id","Panther")]
write.table(hrs_hs_day1_DEGs_wpantherid1, file='hrs_hs_day1_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)

#Day 11 of heat treatment
diff_hrs_hs_day11_exclusive_DE_genes <- read.delim("~/Documents/kallisto_test/hrs_hs_day11_exclusive_DE_genes.txt")
hrs_hs_day11_DEGs_wpantherid <- merge(Cquinoa_392_v1.0.annotation_info_sub, diff_hrs_hs_day11_exclusive_DE_genes, by='target_id')
hrs_hs_day11_DEGs_wpantherid1 <- hrs_hs_day11_DEGs_wpantherid[, c("target_id","Panther")]
write.table(hrs_hs_day11_DEGs_wpantherid1, file='hrs_hs_day11_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)

#Obatin the list of the 394 genes expressed on both days, 
#in both yield impacted treatments, with respective b values
hrs_hs_all_days_DEGs_wpantherid <- merge(hrs_hs_day1_DEGs_wpantherid, hrs_hs_day11_DEGs_wpantherid, by='target_id')
hrs_hs_all_days_DEGs_wpantherid1 <- hrs_hs_all_days_DEGs_wpantherid[, c("target_id","Panther.x", "b.hrs.x", "b.hs.x", "b.hrs.y", "b.hs.y")]
colnames(hrs_hs_all_days_DEGs_wpantherid1)[2:6]<- c("Panther","b.hrs.day.1", "b.hs.day.1", "b.hrs.day.11", "b.hs.day.11")
#Add a column with average b values accross days and treatments, for the statistical enrichment test
hrs_hs_all_days_DEGs_wpantherid1$b.avg <- (hrs_hs_all_days_DEGs_wpantherid1$b.hrs.day.1 + hrs_hs_all_days_DEGs_wpantherid1$b.hs.day.1 + hrs_hs_all_days_DEGs_wpantherid1$b.hrs.day.11 + hrs_hs_all_days_DEGs_wpantherid1$b.hs.day.11)/4
hrs_hs_all_days_DEGs_wpantherid1 <- hrs_hs_all_days_DEGs_wpantherid1[, c("target_id","Panther", "b.avg")]
write.table(hrs_hs_all_days_DEGs_wpantherid1, file='hrs_hs_all_days_DEGs_wpantherid.txt', sep="\t",row.names=F, quote=F)


#Lists ready to run in http://www.pantherdb.org/