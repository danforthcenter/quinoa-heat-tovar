library(tm)
library(ggplot2)
library(gridExtra)
library(plyr)
library(shiny)
library(reshape2)
library(reshape)
library(stringr)

setwd("/Users/mgehan/Documents/github/quinoa-heat-tovar/data/")

experiment<-read.table(file='experiment.information.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
saveRDS(experiment,"experiment.information.rds")

hr1<-read.table(file='all_genes_lists/hr_day1_all_genes.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=FALSE, na.strings=c(""," ","NA"))
hs1<-read.table(file='all_genes_lists/hs_day1_all_genes.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=FALSE, na.strings=c(""," ","NA"))
hrs1<-read.table(file='all_genes_lists/hrs_day1_all_genes.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=FALSE, na.strings=c(""," ","NA"))
hr11<-read.table(file='all_genes_lists/hr_day11_all_genes.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=FALSE, na.strings=c(""," ","NA"))
hs11<-read.table(file='all_genes_lists/hs_day11_all_genes.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=FALSE, na.strings=c(""," ","NA"))
hrs11<-read.table(file='all_genes_lists/hrs_day11_all_genes.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE, quote="", fill=FALSE, na.strings=c(""," ","NA"))

hr1s<-subset(hr1, select=c(target_id,qval,b,Best.hit.arabi.name, PANTHER_GO.slim_biological_process, PANTHER_GO.slim_molecular_function,PANTHER_GO.slim_cellular_component))
colnames(hr1s)<-c('target_id','hr1_qval','hr1_de','best_hit_arabidopsis', "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component" )

hs1s<-subset(hs1, select=c(target_id,qval,b,Best.hit.arabi.name, PANTHER_GO.slim_biological_process, PANTHER_GO.slim_molecular_function,PANTHER_GO.slim_cellular_component))
colnames(hs1s)<-c('target_id','hs1_qval','hs1_de','best_hit_arabidopsis', "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component")

hrs1s<-subset(hrs1, select=c(target_id,qval,b,Best.hit.arabi.name,PANTHER_GO.slim_biological_process, PANTHER_GO.slim_molecular_function,PANTHER_GO.slim_cellular_component))
colnames(hrs1s)<-c('target_id','hrs1_qval','hrs1_de','best_hit_arabidopsis', "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component")

hr11s<-subset(hr11, select=c(target_id,qval,b,Best.hit.arabi.name, PANTHER_GO.slim_biological_process, PANTHER_GO.slim_molecular_function,PANTHER_GO.slim_cellular_component))
colnames(hr11s)<-c('target_id','hr11_qval','hr11_de','best_hit_arabidopsis', "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component")

hs11s<-subset(hs11, select=c(target_id,qval,b,Best.hit.arabi.name, PANTHER_GO.slim_biological_process, PANTHER_GO.slim_molecular_function,PANTHER_GO.slim_cellular_component))
colnames(hs11s)<-c('target_id','hs11_qval','hs11_de','best_hit_arabidopsis', "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component")

hrs11s<-subset(hrs11, select=c(target_id,qval,b,Best.hit.arabi.name, PANTHER_GO.slim_biological_process, PANTHER_GO.slim_molecular_function,PANTHER_GO.slim_cellular_component))
colnames(hrs11s)<-c('target_id','hrs11_qval','hrs11_de','best_hit_arabidopsis', "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component")

hr1.hs1<-merge(hr1s,hs1s,by=c("target_id","best_hit_arabidopsis", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component"), all=TRUE)
hr1.hs1.hrs1<-merge(hr1.hs1,hrs1s,by=c("target_id","best_hit_arabidopsis", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component"), all=TRUE)
hr1.hs1.hrs1.hr11<-merge(hr1.hs1.hrs1,hr11s,by=c("target_id","best_hit_arabidopsis", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component"), all=TRUE)
hr1.hs1.hrs1.hr11.hs11<-merge(hr1.hs1.hrs1.hr11,hs11s,by=c("target_id","best_hit_arabidopsis", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component"), all=TRUE)
hr1.hs1.hrs1.hr11.hs11.hrs11<-merge(hr1.hs1.hrs1.hr11.hs11,hrs11s,by=c("target_id","best_hit_arabidopsis", "PANTHER_GO.slim_biological_process", "PANTHER_GO.slim_molecular_function","PANTHER_GO.slim_cellular_component"), all=TRUE)
hr1.hs1.hrs1.hr11.hs11.hrs11.sub<-subset(hr1.hs1.hrs1.hr11.hs11.hrs11, select=-c(target_id,best_hit_arabidopsis,PANTHER_GO.slim_biological_process, PANTHER_GO.slim_molecular_function,PANTHER_GO.slim_cellular_component))

rnaseqdata<-hr1.hs1.hrs1.hr11.hs11.hrs11[rowSums(is.na(hr1.hs1.hrs1.hr11.hs11.hrs11.sub)) != ncol(hr1.hs1.hrs1.hr11.hs11.hrs11.sub), ]
saveRDS(rnaseqdata,"quinoa-all-data.rds")
