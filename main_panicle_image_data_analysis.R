library(WRS2)
library(ggplot2)

setwd("/Users/jtovar/quinoa_heat/")

main.panicles.full.maturity <- read.csv("~/quinoa_heat/main-panicles-full-maturity.csv")

##Analyze area
#Summary statistics
area.sub <- main.panicles.full.maturity[, c("treatment","area")]
area.sub_melted <- melt(area.sub, id.vars = c("treatment"))
ddply(area.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), median(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=area))+geom_histogram(binwidth = 60000)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(area ~ root_heating*shoot_heating, data = main.panicles.full.maturity)
mcp2a(area ~ root_heating*shoot_heating, data = main.panicles.full.maturity, 
      est = "median")

#Test with Kruskal-Wallis to see which treatment is different from which
kruskal.test(area ~ treatment, data = main.panicles.full.maturity)
pairwise.wilcox.test(main.panicles.full.maturity$area, main.panicles.full.maturity$treatment, p.adjust.method = "BY")

#Graph area
pdf(file="main.panicle.area1.pdf",width = 8,height = 6.5,pointsize = 5,useDingbats = FALSE)
ggplot(main.panicles.full.maturity, aes(x = treatment, y = area, fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Main Panicle Area (pixels)", limits = c(0, 1200000)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Main Panicle Area after Heat Exposure") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


##Analyze width

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=width))+geom_histogram(binwidth = 50)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(width ~ root_heating*shoot_heating, data = main.panicles.full.maturity)

#Test with Kruskal-Wallis to see which treatment is different from which
kruskal.test(width ~ treatment, data = main.panicles.full.maturity)
pairwise.wilcox.test(main.panicles.full.maturity$width, main.panicles.full.maturity$treatment, p.adjust.method = "BY")


##Analyze height

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=height))+geom_histogram(binwidth = 50)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(height ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze hull.area

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=hull.area))+geom_histogram(binwidth = 100000)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(hull.area ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze solidity
#Summary statistics
solidity.sub <- main.panicles.full.maturity[, c("treatment","solidity")]
solidity.sub_melted <- melt(solidity.sub, id.vars = c("treatment"))
ddply(solidity.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), median(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=solidity))+geom_histogram(binwidth = 0.05)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(solidity ~ root_heating*shoot_heating, data = main.panicles.full.maturity)
mcp2a(solidity ~ root_heating*shoot_heating, data = main.panicles.full.maturity, 
      est = "median")

#Test with Kruskal-Wallis to see which treatment is different from which
kruskal.test(solidity ~ treatment, data = main.panicles.full.maturity)
pairwise.wilcox.test(main.panicles.full.maturity$solidity, main.panicles.full.maturity$treatment, p.adjust.method = "BY")

#Graph solidity
pdf(file="main.panicle.solidity.pdf",width = 8,height = 6.5,pointsize = 5,useDingbats = FALSE)
ggplot(main.panicles.full.maturity, aes(x = treatment, y = solidity, fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Main Panicle Solidity (panicle area/hull area)", limits = c(0, 0.8)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Main Panicle Solidity after Heat Exposure") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

##Analyze perimeter

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=perimeter))+geom_histogram(binwidth = 2000)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(perimeter ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze longest_axis

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=longest_axis))+geom_histogram(binwidth = 1000)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(longest_axis ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze center.of.mass.x

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=center.of.mass.x))+geom_histogram(binwidth = 100)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(center.of.mass.x ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze center.of.mass.y

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=center.of.mass.y))+geom_histogram(binwidth = 50)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(center.of.mass.y ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze hull_vertices

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=hull_vertices))+geom_histogram(binwidth = 3)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(hull_vertices ~ root_heating*shoot_heating, data = main.panicles.full.maturity)

#Test with Kruskal-Wallis to see which treatment is different from which
kruskal.test(hull_vertices ~ treatment, data = main.panicles.full.maturity)
pairwise.wilcox.test(main.panicles.full.maturity$hull_vertices, main.panicles.full.maturity$treatment, p.adjust.method = "BY")


##Analyze ellipse_center_x

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=ellipse_center_x))+geom_histogram(binwidth = 50)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(ellipse_center_x ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze ellipse_center_y

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=ellipse_center_y))+geom_histogram(binwidth = 50)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(ellipse_center_y ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze ellipse_major_axis

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=ellipse_major_axis))+geom_histogram(binwidth = 100)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(ellipse_major_axis ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze ellipse_minor_axis

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=ellipse_minor_axis))+geom_histogram(binwidth = 100)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(ellipse_minor_axis ~ root_heating*shoot_heating, data = main.panicles.full.maturity)


##Analyze ellipse_angle

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=ellipse_angle))+geom_histogram(binwidth = 20)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(ellipse_angle ~ root_heating*shoot_heating, data = main.panicles.full.maturity)

#Test with Kruskal-Wallis to see which treatment is different from which
kruskal.test(ellipse_angle ~ treatment, data = main.panicles.full.maturity)
pairwise.wilcox.test(main.panicles.full.maturity$ellipse_angle, main.panicles.full.maturity$treatment, p.adjust.method = "BY")


##Analyze ellipse_eccentricity

#Check if data is parametric
ggplot(main.panicles.full.maturity,aes(x=ellipse_eccentricity))+geom_histogram(binwidth = 0.05)+facet_grid(~treatment)+theme_bw()
#Data is not normally distributed, so use non-parametrc test

med2way(ellipse_eccentricity ~ root_heating*shoot_heating, data = main.panicles.full.maturity)

#Test with Kruskal-Wallis to see which treatment is different from which
kruskal.test(ellipse_eccentricity ~ treatment, data = main.panicles.full.maturity)
pairwise.wilcox.test(main.panicles.full.maturity$ellipse_eccentricity, main.panicles.full.maturity$treatment, p.adjust.method = "BY")

