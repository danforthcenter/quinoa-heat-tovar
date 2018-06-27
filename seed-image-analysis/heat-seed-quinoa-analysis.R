library(ggplot2)

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

# Graph


# Statistics


