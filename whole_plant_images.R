library(ggplot2)
library(WRS2)
library(reshape2)
library(plyr)

setwd("/Users/jtovar/quinoa_heat/")
image_shapes_mature_whole_plants <- read.delim("~/quinoa_heat/image_shapes_mature_whole_plants.txt")

#Normalized area statistics
wp.narea.sub <- image_shapes_mature_whole_plants[, c("treatment","normalized_area_pot_side_squared")]
wp.narea.sub_melted <- melt(wp.narea.sub, id.vars = c("treatment"))
ddply(wp.narea.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), median(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(normalized_area_pot_side_squared ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)

kruskal.test(normalized_area_pot_side_squared ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$normalized_area_pot_side_squared, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Graph normalized_area
pdf(file="whole_plants_narea_pot2.pdf",width = 8,height = 6.5,pointsize = 5,useDingbats = FALSE)
ggplot(image_shapes_mature_whole_plants, aes(x = treatment, y = normalized_area_pot_side_squared, fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Normalized Area", limits = c(0, 35)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Quinoa Plant Size") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Normalize area based on side of pot (linear measurement)
wp.nareal.sub <- image_shapes_mature_whole_plants[, c("treatment","normalized_area_pot_side")]
wp.nareal.sub_melted <- melt(wp.nareal.sub, id.vars = c("treatment"))
ddply(wp.nareal.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(normalized_area_pot_side ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)

kruskal.test(normalized_area_pot_side ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$normalized_area_pot_side, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Graph normalized_area
pdf(file="whole_plants_narea_pot2.pdf",width = 8,height = 6.5,pointsize = 5,useDingbats = FALSE)
ggplot(image_shapes_mature_whole_plants, aes(x = treatment, y = normalized_area_pot_side_squared, fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Normalized Area", limits = c(0, 35)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Quinoa Plant Size") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Normalized hull-area statistics
wp.nharea.sub <- image_shapes_mature_whole_plants[, c("treatment","normalized_hull.area")]
wp.nharea.sub_melted <- melt(wp.nharea.sub, id.vars = c("treatment"))
ddply(wp.nharea.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(normalized_hull.area ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)

kruskal.test(normalized_hull.area ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$normalized_hull.area, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Solidity statistics
wp.solidity.sub <- image_shapes_mature_whole_plants[, c("treatment","solidity")]
wp.solidity.sub_melted <- melt(wp.solidity.sub, id.vars = c("treatment"))
ddply(wp.solidity.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(solidity ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)
kruskal.test(solidity ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$solidity, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Graph solidity
pdf(file="whole_plants_solidity.pdf",width = 8,height = 6.5,pointsize = 5,useDingbats = FALSE)
ggplot(image_shapes_mature_whole_plants, aes(x = treatment, y = solidity, fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Solidity", limits = c(0, 0.6)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Quinoa Plant Solidity") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Normalized Perimeter statistics
wp.perimeter.sub <- image_shapes_mature_whole_plants[, c("treatment","normalized_perimeter")]
wp.perimeter.sub_melted <- melt(wp.perimeter.sub, id.vars = c("treatment"))
ddply(wp.perimeter.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(normalized_perimeter ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)
kruskal.test(normalized_perimeter ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$normalized_perimeter, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Normalized width statistics
wp.width.sub <- image_shapes_mature_whole_plants[, c("treatment","normalized_width")]
wp.width.sub_melted <- melt(wp.width.sub, id.vars = c("treatment"))
ddply(wp.width.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(normalized_width ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)
kruskal.test(normalized_width ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$normalized_width, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Normalized height statistics
wp.height.sub <- image_shapes_mature_whole_plants[, c("treatment","normalized_height")]
wp.height.sub_melted <- melt(wp.height.sub, id.vars = c("treatment"))
ddply(wp.height.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(normalized_height ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)
kruskal.test(normalized_height ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$normalized_height, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Normalized longest axis statistics
wp.longest.axis.sub <- image_shapes_mature_whole_plants[, c("treatment","normalized_longest_axis")]
wp.longest.axis.sub_melted <- melt(wp.longest.axis.sub, id.vars = c("treatment"))
ddply(wp.longest.axis.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(normalized_longest_axis ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)
kruskal.test(normalized_longest_axis ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$normalized_longest_axis, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

#Yellow vs. green plant areas (maturity)
wp.gy.sub <- image_shapes_mature_whole_plants[, c("treatment","percent_yellow")]
wp.gy.sub_melted <- melt(wp.gy.sub, id.vars = c("treatment"))
ddply(wp.gy.sub_melted, c("treatment","variable"), summarise, 
      mean = mean(value), meadian = median(value), sd = sd(value),sem = sd(value)/sqrt(length(value)))

med2way(percent_yellow ~ root_heating*shoot_heating, data = image_shapes_mature_whole_plants)
kruskal.test(percent_yellow ~ treatment, data = image_shapes_mature_whole_plants)
pairwise.wilcox.test(image_shapes_mature_whole_plants$percent_yellow, image_shapes_mature_whole_plants$treatment, p.adjust.method = "BY")

pdf(file="whole_plants_yellowness.pdf",width = 8,height = 6.5,pointsize = 5,useDingbats = FALSE)
ggplot(image_shapes_mature_whole_plants, aes(x = treatment, y = percent_yellow, fill = treatment)) + geom_boxplot(aes(fill = treatment)) + 
  geom_point(color="black", pch=21, size = 2, position=position_dodge(width = 0.6)) + 
  theme_bw() + scale_x_discrete("Treatment") + 
  scale_y_continuous("Percent Yellow", limits = c(0, 100)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)) +
  ggtitle("Quinoa Plant Maturity") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
