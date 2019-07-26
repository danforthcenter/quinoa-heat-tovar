library(sleuth)
library(dplyr)
library(VennDiagram)

#Set working directory and base dir
setwd("/Users/jtovar/Documents/kallisto_test/")
base_dir<-"/Users/jtovar/Documents/kallisto_test/"


#########################################################
########## Samples collected on day 11 of heat ##########
#########################################################

#Get sample ids and base directories
sample_ids11 <- dir(file.path(base_dir,"day_11"))
kallisto_dirs11 <- sapply(sample_ids11, function(id) file.path(base_dir, "day_11", id))

# Read in data
s2c11 <- read.table(file.path(base_dir, "day_11_sample_ids.txt"), header = TRUE, stringsAsFactors=FALSE)

# Read in annotation file
annotation<-read.table(file.path(base_dir, "annotation.txt"), header = TRUE,strip.white = TRUE, stringsAsFactors=FALSE)

# Make annotation with Gene Onthology information from Panther's Arabidopsis thanliana gene list
pantherGeneList <- read.delim("~/Documents/kallisto_test/pantherGeneList.txt")
gene_list<- pantherGeneList[,c("Mapped_IDs","gene_name_gene_symbol","PANTHER_family.subfamily","PANTHER_protein_class","PANTHER_GO.slim_biological_process","PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component")]

#Eliminate duplicate genes in "gene_list"
nonredundant = gene_list[!duplicated(gene_list$Mapped_IDs),]
write.table(nonredundant, file='gene_list_nonredundant.txt', sep="\t",row.names=F, quote=F)
#####Manual step######
###Manually separated all genes into individual rows and deleted duplicated genes.
#Upload manually edited gene list
gene_list_nonredundant <- read.delim("~/Documents/kallisto_test/gene_list_nonredundant.txt")
#Check that there are no more duplicated genes
gene_list_nonredundant[duplicated(gene_list_nonredundant$Mapped_IDs),]
annotation_go<-merge(annotation,gene_list_nonredundant,by.x='Mapped_IDs',all.x=TRUE)

#average absolute tpm for control samples
control_a11<-read.table(file.path(base_dir,"day_11","C16m11","abundance.tsv"),header = TRUE, stringsAsFactors=FALSE)
control_b11<-read.table(file.path(base_dir,"day_11","C28m11","abundance.tsv"),header = TRUE, stringsAsFactors=FALSE)
control_c11<-read.table(file.path(base_dir,"day_11","C30m11","abundance.tsv"),header = TRUE, stringsAsFactors=FALSE)

control_bind11 <- rbind(control_a11,control_b11,control_c11)
control_avg11 <- summarise(group_by(control_bind11,target_id),
                           mean=mean(tpm), sd=sd(tpm))

# Differential expression between control and heated_roots
s2c_hr11<-s2c11[s2c11$condition=='heated_roots'|s2c11$condition=='control',]
sohr <- sleuth_prep(s2c_hr11, ~ condition)
sohr <- sleuth_fit(sohr)
sohr <- sleuth_wt(sohr, 'conditionheated_roots')
results_tablehr <- sleuth_results(sohr, 'conditionheated_roots')
results_orderedhr <- results_tablehr[order(results_tablehr$qval),]
results_orderedhr_ara_hom <- merge(results_orderedhr,annotation_go,by.x='target_id',all.x=TRUE)
write.table(results_orderedhr_ara_hom, file='hr_day11_all_genes.txt', sep="\t",row.names=F, quote=F)
table(results_orderedhr$qval <= 0.05)
write.table(subset(results_orderedhr_ara_hom, qval <= 0.05), file='hr_day11_DE_genes.qval_0.05.txt', sep="\t",row.names=F, quote=F)

diffcontrolhr11<-subset(results_orderedhr_ara_hom, qval <= 0.05)

# Differential expression between control and heated_shoots
s2c_hs11<-s2c11[s2c11$condition=='heated_shoots'|s2c11$condition=='control',]
sohs <- sleuth_prep(s2c_hs11, ~ condition)
sohs <- sleuth_fit(sohs)
sohs <- sleuth_wt(sohs, 'conditionheated_shoots')
results_tablehs <- sleuth_results(sohs, 'conditionheated_shoots')
results_orderedhs <- results_tablehs[order(results_tablehs$qval),]
results_orderedhs_ara_hom <- merge(results_orderedhs,annotation_go,by.x='target_id',all.x=TRUE)
write.table(results_orderedhs_ara_hom, file='hs_day11_all_genes.txt', sep="\t",row.names=F, quote=F)
table(results_orderedhs$qval <= 0.05)
write.table( subset(results_orderedhs_ara_hom, qval <= 0.05), file='hs_day11_DE_genes.qval_0.05.txt', sep="\t",row.names=F, quote=F)

diffcontrolhs11<-subset(results_orderedhs_ara_hom, qval <= 0.05)

# Differential expression between control and heated_roots_and_shoots
s2c_hrs11<-s2c11[s2c11$condition=='heated_roots_and_shoots'|s2c11$condition=='control',]
sohrs <- sleuth_prep(s2c_hrs11, ~ condition)
sohrs <- sleuth_fit(sohrs)
sohrs <- sleuth_wt(sohrs, 'conditionheated_roots_and_shoots')
results_tablehrs <- sleuth_results(sohrs, 'conditionheated_roots_and_shoots')
results_orderedhrs <- results_tablehrs[order(results_tablehrs$qval),]
results_orderedhrs_ara_hom <- merge(results_orderedhrs,annotation_go,by.x='target_id',all.x=TRUE)
write.table(results_orderedhrs_ara_hom, file='hrs_day11_all_genes.txt', sep="\t",row.names=F, quote=F)
table(results_orderedhrs$qval <= 0.05)
write.table( subset(results_orderedhrs_ara_hom, qval <= 0.05), file='hrs_day11_DE_genes.qval_0.05.txt', sep="\t",row.names=F, quote=F)

diffcontrolhrs11<-subset(results_orderedhrs_ara_hom, qval <= 0.05)

#Make a Venn diagram of the genes differentially expressed on day 11 of heat
pdf("day11_DE_genes.pdf")
venn.plot <- venn.diagram(list(diffcontrolhr11$target_id, diffcontrolhs11$target_id, diffcontrolhrs11$target_id), NULL, fill=c("green", "purple", "blue"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Heated Roots", "Heated Shoots", "Heated Roots and Shoots"))
grid.draw(venn.plot)
dev.off()

#Find genes that are differentially expressed in all 3 treatments
diff_day11<-merge(diffcontrolhr11,diffcontrolhrs11,by='target_id')
diff_day11<-merge(diff_day11,diffcontrolhs11,by='target_id')
diff_day11<-merge(diff_day11,control_avg11,by='target_id',all.x=TRUE)
diff_day11 <- subset(diff_day11, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, gene_name_gene_symbol.y, PANTHER_family.subfamily.x, PANTHER_family.subfamily.y, PANTHER_protein_class.x, PANTHER_protein_class.y, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_biological_process.y, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_molecular_function.y, PANTHER_GO.slim_cellular_component.x, PANTHER_GO.slim_cellular_component.y))
colnames(diff_day11)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(diff_day11)[12:21]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(diff_day11)[22:31]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")

#make a file with the differentially expressed genes in all 3 treatments. 
write.table(diff_day11, file='day11_DE_genes_all_treatments.txt', sep="\t",row.names=F, quote=F)


#Find the expression levels for every gene in all 3 treatments
all_day11<-merge(results_orderedhr_ara_hom,results_orderedhrs_ara_hom,by='target_id')
all_day11<-merge(all_day11,results_orderedhs_ara_hom,by='target_id')
all_day11<-merge(all_day11,control_avg11,by='target_id',all.x=TRUE)
all_day11 <- subset(all_day11, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, gene_name_gene_symbol.y, PANTHER_family.subfamily.x, PANTHER_family.subfamily.y, PANTHER_protein_class.x, PANTHER_protein_class.y, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_biological_process.y, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_molecular_function.y, PANTHER_GO.slim_cellular_component.x, PANTHER_GO.slim_cellular_component.y))
colnames(all_day11)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(all_day11)[12:21]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(all_day11)[22:31]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")

#Make a file with every gene's expression data for all 3 treatments 
write.table(all_day11, file='day11_all_genes_all_treatments.txt', sep="\t",row.names=F, quote=F)

#########################################################
########## Samples collected on day 1 of heat ###########
#########################################################

#Get sample ids and base directories
sample_ids1 <- dir(file.path(base_dir,"day_1"))
kallisto_dirs1 <- sapply(sample_ids1, function(id) file.path(base_dir, "day_1", id))

# Read in data
s2c1 <- read.table(file.path(base_dir, "day_1_sample_ids.txt"), header = TRUE, stringsAsFactors=FALSE)

#average absolute tpm for control samples
control_a1<-read.table(file.path(base_dir,"day_1","C20m1","abundance.tsv"),header = TRUE, stringsAsFactors=FALSE)
control_b1<-read.table(file.path(base_dir,"day_1","C22m1","abundance.tsv"),header = TRUE, stringsAsFactors=FALSE)
control_c1<-read.table(file.path(base_dir,"day_1","C27m1","abundance.tsv"),header = TRUE, stringsAsFactors=FALSE)

control_bind1<-rbind(control_a1,control_b1,control_c1)
control_avg1 <- summarise(group_by(control_bind1,target_id),
                          mean=mean(tpm), sd=sd(tpm))

# Differential expression between control and heated_roots
s2c_hr1<-s2c1[s2c1$condition=='heated_roots'|s2c1$condition=='control',]
sohr1 <- sleuth_prep(s2c_hr1, ~ condition)
sohr1 <- sleuth_fit(sohr1)
sohr1 <- sleuth_wt(sohr1, 'conditionheated_roots')
results_tablehr1 <- sleuth_results(sohr1, 'conditionheated_roots')
results_orderedhr1 <- results_tablehr1[order(results_tablehr1$qval),]
results_orderedhr1_ara_hom <- merge(results_orderedhr1,annotation_go,by.x='target_id',all.x=TRUE)
write.table(results_orderedhr1_ara_hom, file='hr_day1_all_genes.txt', sep="\t",row.names=F, quote=F)
table(results_orderedhr1$qval <= 0.05)
write.table( subset(results_orderedhr1_ara_hom, qval <= 0.05), file='hr_day1_DE_genes.qval_0.05.txt', sep="\t",row.names=F, quote=F)

diffcontrolhr1<-subset(results_orderedhr1_ara_hom, qval <= 0.05)

# Differential expression between control and heated_shoots
s2c_hs1<-s2c1[s2c1$condition=='heated_shoots'|s2c1$condition=='control',]
sohs1 <- sleuth_prep(s2c_hs1, ~ condition)
sohs1 <- sleuth_fit(sohs1)
sohs1 <- sleuth_wt(sohs1, 'conditionheated_shoots')
results_tablehs1 <- sleuth_results(sohs1, 'conditionheated_shoots')
results_orderedhs1 <- results_tablehs1[order(results_tablehs1$qval),]
results_orderedhs1_ara_hom <- merge(results_orderedhs1,annotation_go,by.x='target_id',all.x=TRUE)
write.table(results_orderedhs1_ara_hom, file='hs_day1_all_genes.txt', sep="\t",row.names=F, quote=F)
table(results_orderedhs1$qval <= 0.05)
write.table(subset(results_orderedhs1_ara_hom, qval <= 0.05), file='hs_day1_DE_genes.qval_0.05.txt', sep="\t",row.names=F, quote=F)

diffcontrolhs1<-subset(results_orderedhs1_ara_hom, qval <= 0.05)

# Differential expression between control and heated_roots_and_shoots
s2c_hrs1<-s2c1[s2c1$condition=='heated_roots_and_shoots'|s2c1$condition=='control',]
sohrs1 <- sleuth_prep(s2c_hrs1, ~ condition)
sohrs1 <- sleuth_fit(sohrs1)
sohrs1 <- sleuth_wt(sohrs1, 'conditionheated_roots_and_shoots')
results_tablehrs1 <- sleuth_results(sohrs1, 'conditionheated_roots_and_shoots')
results_orderedhrs1 <- results_tablehrs1[order(results_tablehrs1$qval),]
results_orderedhrs1_ara_hom <- merge(results_orderedhrs1,annotation_go,by.x='target_id',all.x=TRUE)
write.table(results_orderedhrs1_ara_hom, file='hrs_day1_all_genes.txt', sep="\t",row.names=F, quote=F)
table(results_orderedhrs1$qval <= 0.05)
write.table(subset(results_orderedhrs1_ara_hom, qval <= 0.05), file='hrs_day1_DE_genes.qval_0.05.txt', sep="\t",row.names=F, quote=F)

diffcontrolhrs1<-subset(results_orderedhrs1_ara_hom, qval <= 0.05)

#Make a Venn diagram of the genes differentially expressed on day 1 of heat treatment
pdf("day1_DE_genes.pdf")
venn.plot1 <- venn.diagram(list(diffcontrolhr1$target_id, diffcontrolhs1$target_id, diffcontrolhrs1$target_id), NULL, fill=c("green", "purple", "blue"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Heated Roots", "Heated Shoots", "Heated Roots and Shoots"))
grid.draw(venn.plot1)
dev.off()

#Find how all genes are expressed in all 3 treatments
all_day1<-merge(results_orderedhr1_ara_hom,results_orderedhrs1_ara_hom,by='target_id')
all_day1<-merge(all_day1,results_orderedhs1_ara_hom,by='target_id')
all_day1<-merge(all_day1,control_avg1,by='target_id',all.x=TRUE)
all_day1 <- subset(all_day1, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, gene_name_gene_symbol.y, PANTHER_family.subfamily.x, PANTHER_family.subfamily.y, PANTHER_protein_class.x, PANTHER_protein_class.y, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_biological_process.y, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_molecular_function.y, PANTHER_GO.slim_cellular_component.x, PANTHER_GO.slim_cellular_component.y))
colnames(all_day1)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(all_day1)[12:21]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(all_day1)[22:31]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")

#Make a file with every gene's expression data for all 3 treatments
write.table(all_day1, file='day1_all_genes_all_treatments.txt', sep="\t",row.names=F, quote=F)

#Find genes that are differentially expressed in all 3 treatments
diff_hr_hrs_hs_day1<-merge(diffcontrolhr1,diffcontrolhrs1,by='target_id')
diff_hr_hrs_hs_day1<-merge(diff_hr_hrs_hs_day1,diffcontrolhs1,by='target_id')
diff_hr_hrs_hs_day1<-merge(diff_hr_hrs_hs_day1,control_avg1,by='target_id',all.x=TRUE)
diff_hr_hrs_hs_day1 <- subset(diff_hr_hrs_hs_day1, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, gene_name_gene_symbol.y, PANTHER_family.subfamily.x, PANTHER_family.subfamily.y, PANTHER_protein_class.x, PANTHER_protein_class.y, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_biological_process.y, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_molecular_function.y, PANTHER_GO.slim_cellular_component.x, PANTHER_GO.slim_cellular_component.y))
colnames(diff_hr_hrs_hs_day1)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(diff_hr_hrs_hs_day1)[12:21]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(diff_hr_hrs_hs_day1)[22:31]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")

#Make a file with the differentially expressed genes in all 3 treatments. 
write.table(diff_hr_hrs_hs_day1, file='day1_DE_genes_all_treatments.txt', sep="\t",row.names=F, quote=F)

###########################################################################
###Make a file with every gene's expression information for every treatment 
###on both days sampled. This file contains all the RNA-seq results
###########################################################################

all <-merge(all_day1,all_day11,by='target_id')
all <- subset(all, select = -c(Best.hit.arabi.name.x, Mapped_IDs.x, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x, mean.x, sd.x, mean.y, sd.y))
colnames(all)[2:11]<- c("pval.hr.day1","qval.hr.day1", "b.hr.day1", "se_b.hr.day1", "mean_obs.hr.day1", "var_obs.hr.day1", "tech_var.hr.day1", "sigma_sq.hr.day1", "smooth_sigma_sq.hr.day1", "final_sigma_sq.hr.day1")
colnames(all)[12:21]<- c("pval.hrs.day1","qval.hrs.day1", "b.hrs.day1", "se_b.hrs.day1", "mean_obs.hrs.day1", "var_obs.hrs.day1", "tech_var.hrs.day1", "sigma_sq.hrs.day1", "smooth_sigma_sq.hrs.day1", "final_sigma_sq.hrs.day1")
colnames(all)[22:31]<- c("pval.hs.day1","qval.hs.day1", "b.hs.day1", "se_b.hs.day1", "mean_obs.hs.day1", "var_obs.hs.day1", "tech_var.hs.day1", "sigma_sq.hs.day1", "smooth_sigma_sq.hs.day1", "final_sigma_sq.hs.day1")
colnames(all)[32:41]<- c("pval.hr.day11","qval.hr.day11", "b.hr.day11", "se_b.hr.day11", "mean_obs.hr.day11", "var_obs.hr.day11", "tech_var.hr.day11", "sigma_sq.hr.day11", "smooth_sigma_sq.hr.day11", "final_sigma_sq.hr.day11")
colnames(all)[42:51]<- c("pval.hrs.day11","qval.hrs.day11", "b.hrs.day11", "se_b.hrs.day11", "mean_obs.hrs.day11", "var_obs.hrs.day11", "tech_var.hrs.day11", "sigma_sq.hrs.day11", "smooth_sigma_sq.hrs.day11", "final_sigma_sq.hrs.day11")
colnames(all)[52:61]<- c("pval.hs.day11","qval.hs.day11", "b.hs.day11", "se_b.hs.day11", "mean_obs.hs.day11", "var_obs.hs.day11", "tech_var.hs.day11", "sigma_sq.hs.day11", "smooth_sigma_sq.hs.day11", "final_sigma_sq.hs.day11")
colnames(all)[62:69]<- c("Mapped_IDs","Best.hit.arabi.name", "gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")

#Make a file with every gene's expression data for all 3 treatments, on both days sampled
write.table(all, file='all_genes_all_days_all_treatments.txt', sep="\t",row.names=F, quote=F)


########################################################################################
#Find genes differentially expressed in all treatments, during both days of heat sampled
########################################################################################

#Find which genes are expressed both in day 1 and day 11 of heating in all treatments
diffall<-merge(diff_hr_hrs_hs_day1, diff_day11, by='target_id')
diffall <- subset(diffall, select = -c(Best.hit.arabi.name.x, mean.x, sd.x, mean.y, sd.y, Mapped_IDs.x, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diffall)[2:11]<- c("pval.hr.day1","qval.hr.day1", "b.hr.day1", "se_b.hr.day1", "mean_obs.hr.day1", "var_obs.hr.day1", "tech_var.hr.day1", "sigma_sq.hr.day1", "smooth_sigma_sq.hr.day1", "final_sigma_sq.hr.day1")
colnames(diffall)[12:21]<- c("pval.hrs.day1","qval.hrs.day1", "b.hrs.day1", "se_b.hrs.day1", "mean_obs.hrs.day1", "var_obs.hrs.day1", "tech_var.hrs.day1", "sigma_sq.hrs.day1", "smooth_sigma_sq.hrs.day1", "final_sigma_sq.hrs.day1")
colnames(diffall)[22:31]<- c("pval.hs.day1","qval.hs.day1", "b.hs.day1", "se_b.hs.day1", "mean_obs.hs.day1", "var_obs.hs.day1", "tech_var.hs.day1", "sigma_sq.hs.day1", "smooth_sigma_sq.hs.day1", "final_sigma_sq.hs.day1")
colnames(diffall)[32:41]<- c("pval.hr.day11","qval.hr.day11", "b.hr.day11", "se_b.hr.day11", "mean_obs.hr.day11", "var_obs.hr.day11", "tech_var.hr.day11", "sigma_sq.hr.day11", "smooth_sigma_sq.hr.day11", "final_sigma_sq.hr.day11")
colnames(diffall)[42:51]<- c("pval.hrs.day11","qval.hrs.day11", "b.hrs.day11", "se_b.hrs.day11", "mean_obs.hrs.day11", "var_obs.hrs.day11", "tech_var.hrs.day11", "sigma_sq.hrs.day11", "smooth_sigma_sq.hrs.day11", "final_sigma_sq.hrs.day11")
colnames(diffall)[52:61]<- c("pval.hs.day11","qval.hs.day11", "b.hs.day11", "se_b.hs.day11", "mean_obs.hs.day11", "var_obs.hs.day11", "tech_var.hs.day11", "sigma_sq.hs.day11", "smooth_sigma_sq.hs.day11", "final_sigma_sq.hs.day11")
colnames(diffall)[62:69] <- c("Mapped_IDs","Best.hit.arabi.name", "gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")

write.table(diffall, file='DE_genes_all_treatments_all_days.txt', sep="\t",row.names=F, quote=F)

##################################################################
#Find genes differentially expressed in only one or two treatments
##################################################################
#############
#For day 1
#############
#Find genes that are expressed in heated_roots and heated_roots_and_shoots treatments
diff_hr_hrs_day1<-merge(diffcontrolhr1,diffcontrolhrs1,by='target_id')
diff_hr_hrs_day1<-merge(diff_hr_hrs_day1,control_avg1,by='target_id',all.x=TRUE)
diff_hr_hrs_day1<-merge(diff_hr_hrs_day1,annotation,by.x='target_id',all.x=TRUE)
diff_hr_hrs_day1 <- subset(diff_hr_hrs_day1, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, mean, sd, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diff_hr_hrs_day1)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(diff_hr_hrs_day1)[12:21]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(diff_hr_hrs_day1)[22:27] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")

write.table(diff_hr_hrs_day1, file='hr_hrs_day1_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes that are expressed in heated_roots and heated_shoots treatments
diff_hr_hs_day1<-merge(diffcontrolhr1,diffcontrolhs1,by='target_id')
diff_hr_hs_day1<-merge(diff_hr_hs_day1,control_avg1,by='target_id',all.x=TRUE)
diff_hr_hs_day1<-merge(diff_hr_hs_day1,annotation,by.x='target_id',all.x=TRUE)
diff_hr_hs_day1 <- subset(diff_hr_hs_day1, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, mean, sd, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diff_hr_hs_day1)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(diff_hr_hs_day1)[12:21]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")
colnames(diff_hr_hs_day1)[22:27] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")

write.table(diff_hr_hs_day1, file='hr_hs_day1_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes that are expressed in heated_roots_and_shoots and heated_shoots treatments
diff_hrs_hs_day1<-merge(diffcontrolhrs1,diffcontrolhs1,by='target_id')
diff_hrs_hs_day1<-merge(diff_hrs_hs_day1,control_avg1,by='target_id',all.x=TRUE)
diff_hrs_hs_day1<-merge(diff_hrs_hs_day1,annotation,by.x='target_id',all.x=TRUE)
diff_hrs_hs_day1 <- subset(diff_hrs_hs_day1, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, mean, sd, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diff_hrs_hs_day1)[2:11]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(diff_hrs_hs_day1)[12:21]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")
colnames(diff_hrs_hs_day1)[22:27] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")

write.table(diff_hrs_hs_day1, file='hrs_hs_day1_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely between heated_roots and heated_roots_and_shoots treatments, 
#excluding genes differentially expressed in all 3 treatments
diff_hr_hrs_day1_exclusive <- diff_hr_hrs_day1[!(diff_hr_hrs_day1$target_id %in% diff_hr_hrs_hs_day1$target_id), ]
write.table(diff_hr_hrs_day1_exclusive, file='hr_hrs_day1_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely between heated_roots and heated_shoots treatments, 
#excluding genes differentially expressed in all 3 treatments
diff_hr_hs_day1_exclusive <- diff_hr_hs_day1[!(diff_hr_hs_day1$target_id %in% diff_hr_hrs_hs_day1$target_id), ]
write.table(diff_hr_hs_day1_exclusive, file='hr_hs_day1_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely between heated_roots_and_shoots and heated_shoots treatments, 
#excluding genes differentially expressed in all 3 treatments
diff_hrs_hs_day1_exclusive <- diff_hrs_hs_day1[!(diff_hrs_hs_day1$target_id %in% diff_hr_hrs_hs_day1$target_id), ]
write.table(diff_hrs_hs_day1_exclusive, file='hrs_hs_day1_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely in the heated_roots treatment, 
#excluding genes differentially expressed in any other treatments
diff_hr_day1<-merge(diffcontrolhr1,control_avg1,by='target_id',all.x=TRUE)
colnames(diff_hr_day1)[3:4] <- c("qval.hr","b.hr") 
diff_hr_day1_exclusive <- diff_hr_day1[ !(diff_hr_day1$target_id %in% diff_hr_hrs_hs_day1$target_id), ]
diff_hr_day1_exclusive <- diff_hr_day1_exclusive[ !(diff_hr_day1_exclusive$target_id %in% diff_hr_hrs_day1_exclusive$target_id), ]
diff_hr_day1_exclusive <- diff_hr_day1_exclusive[ !(diff_hr_day1_exclusive$target_id %in% diff_hr_hs_day1_exclusive$target_id), ]
write.table(diff_hr_day1_exclusive, file='hr_day1_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely in the heated_roots_and_shoots treatment, 
#excluding genes differentially expressed in any other treatments
diff_hrs_day1<-merge(diffcontrolhrs1,control_avg1,by='target_id',all.x=TRUE)
colnames(diff_hrs_day1)[3:4] <- c("qval.hrs","b.hrs")
diff_hrs_day1_exclusive <- diff_hrs_day1[ !(diff_hrs_day1$target_id %in% diff_hr_hrs_hs_day1$target_id), ]
diff_hrs_day1_exclusive <- diff_hrs_day1_exclusive[ !(diff_hrs_day1_exclusive$target_id %in% diff_hr_hrs_day1_exclusive$target_id), ]
diff_hrs_day1_exclusive <- diff_hrs_day1_exclusive[ !(diff_hrs_day1_exclusive$target_id %in% diff_hrs_hs_day1_exclusive$target_id), ]
write.table(diff_hrs_day1_exclusive, file='hrs_day1_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely in the heated_shoots treatment, 
#excluding genes differentially expressed in any other treatments
diff_hs_day1<-merge(diffcontrolhs1,control_avg1,by='target_id',all.x=TRUE)
colnames(diff_hs_day1)[3:4] <- c("qval.hs","b.hs")
diff_hs_day1_exclusive <- diff_hs_day1[ !(diff_hs_day1$target_id %in% diff_hr_hrs_hs_day1$target_id), ]
diff_hs_day1_exclusive <- diff_hs_day1_exclusive[ !(diff_hs_day1_exclusive$target_id %in% diff_hrs_hs_day1_exclusive$target_id), ]
diff_hs_day1_exclusive <- diff_hs_day1_exclusive[ !(diff_hs_day1_exclusive$target_id %in% diff_hr_hs_day1_exclusive$target_id), ]
write.table(diff_hs_day1_exclusive, file='hs_day1_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#############
#For day 11
#############
#Find genes that are expressed in heated_roots and heated_roots_and_shoots treatments
diff_hr_hrs_day11<-merge(diffcontrolhr11,diffcontrolhrs11,by='target_id')
diff_hr_hrs_day11<-merge(diff_hr_hrs_day11,control_avg11,by='target_id',all.x=TRUE)
diff_hr_hrs_day11<-merge(diff_hr_hrs_day11,annotation,by.x='target_id',all.x=TRUE)
diff_hr_hrs_day11 <- subset(diff_hr_hrs_day11, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, mean, sd, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diff_hr_hrs_day11)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(diff_hr_hrs_day11)[12:21]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(diff_hr_hrs_day11)[22:27] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")
write.table(diff_hr_hrs_day11, file='hr_hrs_day11_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes that are expressed in heated_roots and heated_shoots treatments
diff_hr_hs_day11<-merge(diffcontrolhr11,diffcontrolhs11,by='target_id')
diff_hr_hs_day11<-merge(diff_hr_hs_day11,control_avg11,by='target_id',all.x=TRUE)
diff_hr_hs_day11<-merge(diff_hr_hs_day11,annotation,by.x='target_id',all.x=TRUE)
diff_hr_hs_day11 <- subset(diff_hr_hs_day11, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, mean, sd, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diff_hr_hs_day11)[2:11]<- c("pval.hr","qval.hr", "b.hr", "se_b.hr", "mean_obs.hr", "var_obs.hr", "tech_var.hr", "sigma_sq.hr", "smooth_sigma_sq.hr", "final_sigma_sq.hr")
colnames(diff_hr_hs_day11)[12:21]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")
colnames(diff_hr_hs_day11)[22:27] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")
write.table(diff_hr_hs_day11, file='hr_hs_day11_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes that are expressed in heated_roots_and_shoots and heated_shoots treatments
diff_hrs_hs_day11<-merge(diffcontrolhrs11,diffcontrolhs11,by='target_id')
diff_hrs_hs_day11<-merge(diff_hrs_hs_day11,control_avg11,by='target_id',all.x=TRUE)
diff_hrs_hs_day11<-merge(diff_hrs_hs_day11,annotation,by.x='target_id',all.x=TRUE)
diff_hrs_hs_day11 <- subset(diff_hrs_hs_day11, select = -c(Best.hit.arabi.name.x, Best.hit.arabi.name.y, mean, sd, Mapped_IDs.x, Mapped_IDs.y, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diff_hrs_hs_day11)[2:11]<- c("pval.hrs","qval.hrs", "b.hrs", "se_b.hrs", "mean_obs.hrs", "var_obs.hrs", "tech_var.hrs", "sigma_sq.hrs", "smooth_sigma_sq.hrs", "final_sigma_sq.hrs")
colnames(diff_hrs_hs_day11)[12:21]<- c("pval.hs","qval.hs", "b.hs", "se_b.hs", "mean_obs.hs", "var_obs.hs", "tech_var.hs", "sigma_sq.hs", "smooth_sigma_sq.hs", "final_sigma_sq.hs")
colnames(diff_hrs_hs_day11)[22:27] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component")
write.table(diff_hrs_hs_day11, file='hrs_hs_day11_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely between heated_roots and heated_roots_and_shoots treatments, 
#excluding genes differentially expressed in all 3 treatments
diff_hr_hrs_day11_exclusive <- diff_hr_hrs_day11[ !(diff_hr_hrs_day11$target_id %in% diff_day11$target_id), ]
write.table(diff_hr_hrs_day11_exclusive, file='hr_hrs_day11_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely between heated_roots and heated_shoots treatments, 
#excluding genes differentially expressed in all 3 treatments
diff_hr_hs_day11_exclusive <- diff_hr_hs_day11[ !(diff_hr_hs_day11$target_id %in% diff_day11$target_id), ]
write.table(diff_hr_hs_day11_exclusive, file='hr_hs_day11_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely between heated_roots_and_shoots and heated_shoots treatments, 
#excluding genes differentially expressed in all 3 treatments
diff_hrs_hs_day11_exclusive <- diff_hrs_hs_day11[ !(diff_hrs_hs_day11$target_id %in% diff_day11$target_id), ]
write.table(diff_hrs_hs_day11_exclusive, file='hrs_hs_day11_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely in the heated_roots treatment, 
#excluding genes differentially expressed in any other treatments
diff_hr_day11<-merge(diffcontrolhr11,control_avg11,by='target_id',all.x=TRUE)
colnames(diff_hr_day11)[3:4] <- c("qval.hr","b.hr") 
diff_hr_day11_exclusive <- diff_hr_day11[!(diff_hr_day11$target_id %in% diff_day11$target_id), ]
diff_hr_day11_exclusive <- diff_hr_day11_exclusive[!(diff_hr_day11_exclusive$target_id %in% diff_hr_hrs_day11_exclusive$target_id), ]
diff_hr_day11_exclusive <- diff_hr_day11_exclusive[!(diff_hr_day11_exclusive$target_id %in% diff_hr_hs_day11_exclusive$target_id), ]
write.table(diff_hr_day11_exclusive, file='hr_day11_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely in the heated_roots_and_shoots treatment, 
#excluding genes differentially expressed in any other treatments
diff_hrs_day11<-merge(diffcontrolhrs11,control_avg11,by='target_id',all.x=TRUE)
colnames(diff_hrs_day11)[3:4] <- c("qval.hrs", "b.hrs") 
diff_hrs_day11_exclusive <- diff_hrs_day11[!(diff_hrs_day11$target_id %in% diff_day11$target_id), ]
diff_hrs_day11_exclusive <- diff_hrs_day11_exclusive[!(diff_hrs_day11_exclusive$target_id %in% diff_hr_hrs_day11_exclusive$target_id), ]
diff_hrs_day11_exclusive <- diff_hrs_day11_exclusive[!(diff_hrs_day11_exclusive$target_id %in% diff_hrs_hs_day11_exclusive$target_id), ]
write.table(diff_hrs_day11_exclusive, file='hrs_day11_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Find genes differentially expressed uniquely in the heated_shoots treatment, 
#excluding genes differentially expressed in any other treatments
diff_hs_day11<-merge(diffcontrolhs11,control_avg11,by='target_id',all.x=TRUE)
colnames(diff_hs_day11)[3:4] <- c("qval.hs", "b.hs") 
diff_hs_day11_exclusive <- diff_hs_day11[!(diff_hs_day11$target_id %in% diff_day11$target_id), ]
diff_hs_day11_exclusive <- diff_hs_day11_exclusive[!(diff_hs_day11_exclusive$target_id %in% diff_hrs_hs_day11_exclusive$target_id), ]
diff_hs_day11_exclusive <- diff_hs_day11_exclusive[!(diff_hs_day11_exclusive$target_id %in% diff_hr_hs_day11_exclusive$target_id), ]
write.table(diff_hs_day11_exclusive, file='hs_day11_exclusive_DE_genes.txt', sep="\t",row.names=F, quote=F)

##########################################################################
#Find genes differentially expressed in both treatments with reduced yield
#during both days of heat treatment sampled
##########################################################################

#Genes differentially expressed in both heated_roots_and_shoots and in heated_shoots treatments
#on both days 1 and 11 of heat treatment
diff_hrs_hs_both_days <- merge(diff_hrs_hs_day1_exclusive, diff_hrs_hs_day11_exclusive, by='target_id')
diff_hrs_hs_both_days <- subset(diff_hrs_hs_both_days, select = -c(Best.hit.arabi.name.x, Mapped_IDs.x, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(diff_hrs_hs_both_days)[2:11] <- c("pval.hrs.day1","qval.hrs.day1", "b.hrs.day1", "se_b.hrs.day1", "mean_obs.hrs.day1", "var_obs.hrs.day1", "tech_var.hrs.day1", "sigma_sq.hrs.day1", "smooth_sigma_sq.hrs.day1", "final_sigma_sq.hrs.day1")
colnames(diff_hrs_hs_both_days)[12:21] <- c("pval.hs.day1","qval.hs.day1", "b.hs.day1", "se_b.hs.day1", "mean_obs.hs.day1", "var_obs.hs.day1", "tech_var.hs.day1", "sigma_sq.hs.day1", "smooth_sigma_sq.hs.day1", "final_sigma_sq.hs.day1")
colnames(diff_hrs_hs_both_days)[22:31] <- c("pval.hrs.day11","qval.hrs.day11", "b.hrs.day11", "se_b.hrs.day11", "mean_obs.hrs.day11", "var_obs.hrs.day11", "tech_var.hrs.day11", "sigma_sq.hrs.day11", "smooth_sigma_sq.hrs.day11", "final_sigma_sq.hrs.day11")
colnames(diff_hrs_hs_both_days)[32:41] <- c("pval.hs.day11","qval.hs.day11", "b.hs.day11", "se_b.hs.day11", "mean_obs.hs.day11", "var_obs.hs.day11", "tech_var.hs.day11", "sigma_sq.hs.day11", "smooth_sigma_sq.hs.day11", "final_sigma_sq.hs.day11")
colnames(diff_hrs_hs_both_days)[42:49] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component", "Best.hit.arabi.name", "Mapped_IDs")
write.table(diff_hrs_hs_both_days, file='hrs_hs_all_days_DE_genes.txt', sep="\t",row.names=F, quote=F)

#Make a Venn diagram of the differentially expressed genes in 
#both heated_roots_and_shoots and heated_shoots treatments in both days of heat
pdf("hrs_hs_DE_genes_all_days.pdf")
venn.plot <- venn.diagram(list(diff_hrs_hs_day1_exclusive$target_id, diff_hrs_hs_day11_exclusive$target_id), NULL, fill=c("orchid", "deepskyblue"), alpha=c(0.5,0.5), cex = 2, cat.cex = 1.5, cat.fontface=1, cat.dist = 0.05, category.names=c("Day 1", "Day 11"), margin = 0.35, main = "Differentially Expressed Genes in \nHeated Roots and Shoots and Heated Shoots Treatments", main.cex = 1.5, main.pos = c(0.5,0.8))
grid.draw(venn.plot)
dev.off()


##########################################################################
#Find transcription factors homologs in treatments heated_roots_and_shoots
#and heated_shoots, the 2 treatments with reduced yield
##########################################################################

#Find transcription factors differentially expressed in both heated_roots_and_shoots and heated_shoots
#Arabidopsis transcription factors list (Ath_TF_list) downloaded from http://planttfdb.cbi.pku.edu.cn/index.php?sp=Ath
Ath_TF_list <- read.delim("~/Documents/kallisto_test/Ath_TF_list")
Ath_TF_list$Best.hit.arabi.name <- Ath_TF_list$TF_ID
ara_tf_homologs_hrs_hs_day1<-merge(diff_hrs_hs_day1_exclusive,Ath_TF_list,by='Best.hit.arabi.name')
write.table(ara_tf_homologs_hrs_hs_day1, file='ara_tf_homologs_hrs_hs_day1.txt', sep="\t",row.names=F, quote=F)
ara_tf_homologs_hrs_hs_day11<-merge(diff_hrs_hs_day11_exclusive,Ath_TF_list,by='Best.hit.arabi.name')
write.table(ara_tf_homologs_hrs_hs_day11, file='ara_tf_homologs_hrs_hs_day11.txt', sep="\t",row.names=F, quote=F)

#Find transcription factors differentially expressed in both heated_roots_and_shoots and heated_shoots
#in both day_1 and day_11 of heating
ara_tf_homologs_hrs_hs_day1_and_day11<-merge(ara_tf_homologs_hrs_hs_day1,ara_tf_homologs_hrs_hs_day11,by='target_id')
ara_tf_homologs_hrs_hs_day1_and_day11 <- subset(ara_tf_homologs_hrs_hs_day1_and_day11, select = -c(Best.hit.arabi.name.y, TF_ID.x, Gene_ID.x, Family.x, Mapped_IDs.x, gene_name_gene_symbol.x, PANTHER_family.subfamily.x, PANTHER_protein_class.x, PANTHER_GO.slim_biological_process.x, PANTHER_GO.slim_molecular_function.x, PANTHER_GO.slim_cellular_component.x))
colnames(ara_tf_homologs_hrs_hs_day1_and_day11)[2] <- "Best.hit.arabi.name"
colnames(ara_tf_homologs_hrs_hs_day1_and_day11)[3:12] <- c("pval.hrs.day1", "qval.hrs.day1", "b.hrs.day1", "se_b.hrs.day1", "mean_obs.hrs.day1", "var_obs.hrs.day1", "tech_var.hrs.day1", "sigma_sq.hrs.day1", "smooth_sigma_sq.hrs.day1", "final_sigma_sq.hrs.day1")
colnames(ara_tf_homologs_hrs_hs_day1_and_day11)[13:22] <- c("pval.hs.day1", "qval.hs.day1", "b.hs.day1", "se_b.hs.day1", "mean_obs.hs.day1", "var_obs.hs.day1", "tech_var.hs.day1", "sigma_sq.hs.day1", "smooth_sigma_sq.hs.day1", "final_sigma_sq.hs.day1")
colnames(ara_tf_homologs_hrs_hs_day1_and_day11)[23:32] <- c("pval.hrs.day11", "qval.hrs.day11", "b.hrs.day11", "se_b.hrs.day11", "mean_obs.hrs.day11", "var_obs.hrs.day11", "tech_var.hrs.day11", "sigma_sq.hrs.day11", "smooth_sigma_sq.hrs.day11", "final_sigma_sq.hrs.day11")
colnames(ara_tf_homologs_hrs_hs_day1_and_day11)[33:42] <- c("pval.hs.day11", "qval.hs.day11", "b.hs.day11", "se_b.hs.day11", "mean_obs.hs.day11", "var_obs.hs.day11", "tech_var.hs.day11", "sigma_sq.hs.day11", "smooth_sigma_sq.hs.day11", "final_sigma_sq.hs.day11")
colnames(ara_tf_homologs_hrs_hs_day1_and_day11)[43:52] <- c("gene_name_gene_symbol", "PANTHER_family.subfamily", "PANTHER_protein_class", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function", "PANTHER_GO.slim_cellular_component", "Mapped_IDs", "TF_ID", "Gene_ID", "Family")
write.table(ara_tf_homologs_hrs_hs_day1_and_day11, file='ara_tf_homologs_hrs_hs_all_days.txt', sep="\t",row.names=F, quote=F)

dev.off()
