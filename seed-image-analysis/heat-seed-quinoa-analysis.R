library(ggplot2)
library(lsmeans)

this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

############################################
# Read data and format for analysis
############################################

# Read weight data
weight.data = read.table(file="imaged-seed-weight.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)

# Data extracted from images
filenames = list.files("./seed-results/")

seed.avg<-data.frame(id=character(),avg.area=numeric(), num.seeds=numeric(), stringsAsFactors=FALSE)

# normalize area to marker area, get average seed area, find number of seeds in image (number of rows of data)
for (f in filenames){
  name = paste0("./seed-results/", f)
  seed.data=read.table(file=name,sep="\t",header=TRUE, stringsAsFactors = FALSE)
  
  seed.data$area.normalized<-NA
  seed.data$area.normalized<-seed.data$area/seed.data$marker_area
  
  geno=seed.data$genotype[1]
  avg= mean(seed.data$area.normalized)
  rown=nrow(seed.data)
  data1=c(geno, avg, rown)
  
  print(data1)
  seed.avg=rbind(seed.avg,data1,stringsAsFactors=FALSE)
}

colnames(seed.avg)<-c("id", "avg.area","num.seeds")
seed.merge=merge(weight.data,seed.avg, by.x="result_name",by.y="id", all.x=TRUE)

seed.merge$est.seed.weight.image<-NA
seed.merge$est.seed.weight.image<-seed.merge$sample_weight/as.numeric(seed.merge$num.seeds)

seed.merge$est.total.seeds<-NA
seed.merge$est.total.seeds<-round((seed.merge$yield_g/seed.merge$est.seed.weight.image),0)

# Graph Seed Area
pdf(file="seed.area_round2.pdf",width = 8,height = 10,pointsize = 5,useDingbats = FALSE)
ggplot(seed.merge, aes(x = treatment, y = as.numeric(avg.area), fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Normalized Area (seed pixels/size marker area)", limits = c(0, 0.1)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Quinoa Seed Area after Heat Exposure") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#ANOVA comparisson
area_mod <- glm(data = seed.merge, as.numeric(avg.area) ~ 0+treatment)
pairs(lsmeans(area_mod, specs = "treatment"))

# Graph est seed number
pdf(file="seed.number_round2.pdf",width = 8,height = 10,pointsize = 5,useDingbats = FALSE)
ggplot(seed.merge, aes(x = treatment, y = as.numeric(est.total.seeds), fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Estimated Seed Number (yield(g)/estimated seed weight)", limits = c(0, 8000)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Quinoa Seed Number after Heat Exposure") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#ANOVA comparisson
area_mod <- glm(data = seed.merge, as.numeric(est.total.seeds) ~ 0+treatment)
pairs(lsmeans(area_mod, specs = "treatment"))

# Graph est seed number
pdf(file="seed.weight_round2.pdf",width = 8,height = 10,pointsize = 5,useDingbats = FALSE)
ggplot(seed.merge, aes(x = treatment, y = as.numeric(est.seed.weight.image), fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Estimated Seed Weight (seed weight(g)/number seeds)", limits = c(0, 0.015)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Quinoa Seed Weight after Heat Exposure") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#ANOVA comparisson
area_mod <- glm(data = seed.merge, as.numeric(est.seed.weight.image) ~ 0+treatment)
pairs(lsmeans(area_mod, specs = "treatment"))

