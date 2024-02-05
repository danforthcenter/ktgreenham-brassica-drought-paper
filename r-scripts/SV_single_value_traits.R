library(tidyverse)
library(ggridges)
library(ggpubr)

setwd("~/Documents/brassica-image-analysis/shape_data/")

############################
# read in data and currate #
############################

z500 <- read.csv("SV_z500_2-single-value-traits.csv")
z500 <- z500 %>% rename("camera" = "frame", "frame" = "camera")
z1500 <- read.csv("VIS_SV_z1500_v2-single-value-traits.csv")
z2500 <- read.csv("VIS_SV_z2500_v2-single-value-traits.csv")

# read in barcode file with genotype and treatment information
barcodes <- read.csv("~/Documents/brassica-image-analysis/barcodes.csv") #read in the barcodes and genotype info
barcodes = select(barcodes, Barcodes, Genotype.ID, Treatment.2, Timepoint, Replicate)

# filter the data files to get just the columns you need
filtered500 <- select(z500, frame, zoom, timestamp, image, plantbarcode, area, convex_hull_area, solidity, perimeter,
                      width, height, longest_path, ellipse_major_axis, ellipse_minor_axis, ellipse_angle, ellipse_eccentricity,
                      hue_circular_mean, hue_circular_std, hue_median)
filtered1500 <- select(z1500, frame, zoom, timestamp, image, plantbarcode, area, convex_hull_area, solidity, perimeter,
                       width, height, longest_path, ellipse_major_axis, ellipse_minor_axis, ellipse_angle, ellipse_eccentricity,
                       hue_circular_mean, hue_circular_std, hue_median)
filtered2500 <- select(z2500, frame, zoom, timestamp, image, plantbarcode, area, convex_hull_area, solidity, perimeter,
                       width, height, longest_path, ellipse_major_axis, ellipse_minor_axis, ellipse_angle, ellipse_eccentricity,
                       hue_circular_mean, hue_circular_std, hue_median)

TS <- c("date","time") #make a vector for the column names of each part of the timestamp
filtered500 <- separate(filtered500, col = timestamp, into = c("date","time"), sep = " ")
filtered1500 <- separate(filtered1500, col = timestamp, into = c("date","time"), sep = " ")
filtered2500 <- separate(filtered2500, col = timestamp, into = c("date","time"), sep = " ") #separate each component of timestamp


#calculate days after first image
startdate <- as.POSIXct("2016-01-13")
filtered500$days=NA
filtered1500$days = NA
filtered2500$days = NA
filtered500$days = as.integer(difftime(filtered500$date, startdate, units = "days"))
filtered1500$days = as.integer(difftime(filtered1500$date, startdate, units = "days"))
filtered2500$days = as.integer(difftime(filtered2500$date, startdate, units = "days"))

# zoom correction
# the zoom correction factors were cacluated from another script and is based on a size marker
filtered500 <- filtered500 %>% mutate_at(c("area", "convex_hull_area"), list(calib = ~ ./1229.355))
filtered500 <- filtered500 %>% mutate_at(c("perimeter", "width", "height", "longest_path", "ellipse_major_axis", "ellipse_minor_axis"), 
                                         list(calib = ~ ./33.95782))
filtered1500 <- filtered1500 %>% mutate_at(c("area", "convex_hull_area"), list(calib = ~ ./2610.580))
filtered1500 <- filtered1500 %>% mutate_at(c("perimeter", "width", "height", "longest_path", "ellipse_major_axis", "ellipse_minor_axis"), 
                                         list(calib = ~ ./49.41139))
filtered2500 <- filtered2500 %>% mutate_at(c("area", "convex_hull_area"), list(calib = ~ ./5543.662))
filtered2500 <- filtered2500 %>% mutate_at(c("perimeter", "width", "height", "longest_path", "ellipse_major_axis", "ellipse_minor_axis"), 
                                           list(calib = ~ ./71.89759))

# separate each filtered DF by imaging day and match barcodes 
DF0 <- filter(filtered2500, filtered2500$days == "0")
DF0 = merge(DF0, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF0 = drop_na(DF0)
DF1 <- filter(filtered2500, filtered2500$days == "1")
DF1 = merge(DF1, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF1 = drop_na(DF1)
DF2 <- filter(filtered2500, filtered2500$days == "2")
DF2 = merge(DF2, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF2 = drop_na(DF2)
DF3 <- filter(filtered2500, filtered2500$days == "3")
DF3 = merge(DF3, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF3 = drop_na(DF3)
DF4 <- filter(filtered2500, filtered2500$days == "4")
DF4 = merge(DF4, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF4 = drop_na(DF4)
DF5 <- filter(filtered2500, filtered2500$days == "5")
DF5 = merge(DF5, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF5 = drop_na(DF5)
DF6 <- filter(filtered2500, filtered2500$days == "6")
DF6 = merge(DF6, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF6 = drop_na(DF6)
DF7 <- filter(filtered2500, filtered2500$days == "7")
DF7 = merge(DF7, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF7 = drop_na(DF7)
DF8 <- filter(filtered2500, filtered2500$days == "8")
DF8 = merge(DF8, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF8 = drop_na(DF8)
DF9 <- filter(filtered2500, filtered2500$days == "9")
DF9 = merge(DF9, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF9 = drop_na(DF9)
DF10 <- filter(filtered2500, filtered2500$days == "10")
DF10 = merge(DF10, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF10 = drop_na(DF10)
DF10.1 <- filter(filtered1500, filtered1500$days == "10")
DF10.1 = merge(DF10, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF10.1 = drop_na(DF10)
DF11 <- filter(filtered1500, filtered1500$days == "11")
DF11 = merge(DF11, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF11 = drop_na(DF11)
DF12 <- filter(filtered1500, filtered1500$days == "12")
DF12 = merge(DF12, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF12 = drop_na(DF12)
DF13 <- filter(filtered1500, filtered1500$days == "13")
DF13 = merge(DF13, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF13 = drop_na(DF13)
DF13.1 <- filter(filtered500, filtered500$days == "13")
DF13.1 = merge(DF13.1, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF13.1 = drop_na(DF13.1)
DF14 <- filter(filtered500, filtered500$days == "14")
DF14 = merge(DF14, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF14 = drop_na(DF14)
DF15 <- filter(filtered500, filtered500$days == "15")
DF15 = merge(DF15, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF15 = drop_na(DF15)
DF16 <- filter(filtered500, filtered500$days == "16")
DF16 = merge(DF16, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF16 = drop_na(DF16)
DF17 <- filter(filtered500, filtered500$days == "17")
DF17 = merge(DF17, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF17 = drop_na(DF17)
DF18 <- filter(filtered500, filtered500$days == "18")
DF18 = merge(DF18, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF18 = drop_na(DF18)
DF19 <- filter(filtered500, filtered500$days == "19")
DF19 = merge(DF19, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF19 = drop_na(DF19)
DF20 <- filter(filtered500, filtered500$days == "20")
DF20 = merge(DF20, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF20 = drop_na(DF20)
DF21 <- filter(filtered500, filtered500$days == "21")
DF21 = merge(DF21, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF21 = drop_na(DF21)
DF22 <- filter(filtered500, filtered500$days == "22")
DF22 = merge(DF22, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF22 = drop_na(DF22)
DF23 <- filter(filtered500, filtered500$days == "23")
DF23 = merge(DF23, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF23 = drop_na(DF23)
DF24 <- filter(filtered500, filtered500$days == "24")
DF24 = merge(DF24, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF24 = drop_na(DF24)
DF25 <- filter(filtered500, filtered500$days == "25")
DF25 = merge(DF25, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF25 = drop_na(DF25)
DF26 <- filter(filtered500, filtered500$days == "26")
DF26 = merge(DF26, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF26 = drop_na(DF26)
DF27 <- filter(filtered500, filtered500$days == "27")
DF27 = merge(DF27, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF27 = drop_na(DF27)
DF28 <- filter(filtered500, filtered500$days == "28")
DF28 = merge(DF28, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF28 = drop_na(DF28)
DF <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF10.1,DF11,DF12,DF13,
                DF13.1,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28)

rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF10.1,DF11,DF12,DF13,
   DF13.1,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28)

write.csv(DF, "~/Documents/brassica-image-analysis/shape_data/VIS_SV_data_compiled_currated_zoomcorrected.csv", row.names = FALSE)


###################################################
# read in compiled and currated data for plotting #
###################################################
setwd("~/Documents/brassica-image-analysis/shape_data/")

DF<- read.csv("VIS_SV_data_compiled_currated_zoomcorrected.csv") # read in data

DF$Treatment.2 <- as.factor(DF$Treatment.2) # make sure the treatments are factors

avgDF <- DF %>% group_by(Genotype.ID, Treatment.2, days, Replicate, Timepoint) %>% 
  summarize(area = sum(area_calib), hue_mean = mean(hue_circular_mean)) #add 0 and 90 deg frame values together
avgDF <- ungroup(avgDF)

# filter DF to have a DF of only the parent genotypes
parents <- filter(DF, Genotype.ID == "Ab4-28"|Genotype.ID == "Al1-29"| Genotype.ID == "Av1-30"|
                    Genotype.ID == "Br2-31"| Genotype.ID == "Ca4-32"| Genotype.ID == "Da2-44"|
                    Genotype.ID == "DH12.6-33"| Genotype.ID == "DH20"| Genotype.ID == "Gr2-35"|
                    Genotype.ID == "Mu4-37"| Genotype.ID == "Ne3-38"| Genotype.ID == "Qu1-39"|
                    Genotype.ID == "Se2-40"| Genotype.ID == "st1-45"|
                    Genotype.ID == "Yu1-42"| Genotype.ID == "Ze2")

###################################################################
# This list of parents is for making the parent/hybrid list in the for loop below
parents_list <- c("Ab4-28", "Al1-29", "Av1-30", "Br2-31", "Ca4-32", "Da2-44", "DH12.6-33", "DH20",
             "Gr2-35", "Mu4-37", "Ne3-38", "Qu1-39", "Se2-40", "st1-45", "Yu1-42", "Ze2")

avgDF$PHlist <- avgDF$Genotype.ID %in% parents_list #make a column that compares the genotype column to the parents list
avgDF <- as.data.frame(avgDF)

avgDF$PH <- NA
for (i in seq_along(avgDF$PHlist)){
    if (avgDF[i,'PHlist'] == TRUE){
      avgDF[i, 'PH'] <- "parent"
    } else if (avgDF[i,'PHlist'] == FALSE){
      avgDF[i,'PH'] <- "hybrid"}} #label genotypes as hybrid or parent

###################################################
#plot hue circular mean
#cols <- c("20"="#99CCFF", "100" = "#3366CC")
ggplot(data = parents, aes(x = days, y = hue_mean , group = Treatment.2)) +
  geom_smooth(method="loess", aes(color = Treatment.2))+
  facet_wrap("Genotype.ID")+
  labs(title = "hue_circular mean SV", x = "days", y = "mean hue angle")+
  scale_x_continuous(breaks = seq(0,28, by = 4))+
  scale_color_manual(name= "Treatment", breaks = c("20","100"), labels= c("Drought","Control"), values = cal_palette(name = "sierra2", type = "discrete"))+
  theme(axis.title.x=element_text(hjust = 0.5,face="bold"),
        axis.title.y=element_text(face="bold"),
        strip.text.x = element_text(size = 10, margin = margin(t=1, b=1), color = "black", hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "grey"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        axis.ticks = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_text(size = 9))
dev.copy(png, "./SV_nirmean.png", width=1000, height=700, res=120)
dev.off()

plotDF <- filter(parents, days < 27)
p <- ggline(data = plotDF, x = 'days', y = 'hue_circular_mean', group = 'Treatment.2', color = 'Treatment.2', 
       add = c("mean_se"), palette = c("#E7B800","#00AFBB"), facet.by = 'Genotype.ID',
       xlab = "Days", ylab = "Hue circular mean (degrees)")
ggpar(p, legend.title = "Treatment")

###################################################
#plot NIR means
NIR <- read.csv('~/Documents/brassica-image-analysis/NIR_data/side_view/NIR_SV_singlevalue_data_compiled_currated.csv')

NIR$Treatment.2 <- as.factor(NIR$Treatment.2) # make sure the treatments are factors

#avgNIR <- NIR %>% group_by(Genotype.ID, Treatment.2, days, Replicate, Timepoint) %>% 
#  summarize(nir_mean = mean(nir_mean)) #average 0 and 90 deg frame values together
#avgNIR <- ungroup(avgNIR)

# filter DF to have a DF of only the parent genotypes
NIRparents <- filter(NIR, Genotype.ID == "Ab4-28"|Genotype.ID == "Al1-29"| Genotype.ID == "Av1-30"|
                    Genotype.ID == "Br2-31"| Genotype.ID == "Ca4-32"| Genotype.ID == "Da2-44"|
                    Genotype.ID == "DH12.6-33"| Genotype.ID == "DH20"| Genotype.ID == "Gr2-35"|
                    Genotype.ID == "Mu4-37"| Genotype.ID == "Ne3-38"| Genotype.ID == "Qu1-39"|
                    Genotype.ID == "Se2-40"| Genotype.ID == "st1-45"|
                    Genotype.ID == "Yu1-42"| Genotype.ID == "Ze2")

# we are only using NIR data from the last zoom setting
NIRparents <- filter(NIRparents, days > 14 & days < 27)
NIRparents$Genotype.ID <- gsub("Ab4-28", "Ab", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Al1-29", "Al", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Av1-30", "Av", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Br2-31", "Br", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Ca4-32", "Ca", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Da2-44", "Da", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("DH12.6-33", "DH12", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Gr2-35", "Gr", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Mu4-37", "Mu", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Ne3-38", "Ne", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Qu1-39", "Qu", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Se2-40", "Se", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("st1-45", "St", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Yu1-42", "Yu", NIRparents$Genotype.ID)
NIRparents$Genotype.ID <- gsub("Ze2", "Ze", NIRparents$Genotype.ID)

NIRparents$Treatment.2 <- gsub("20", "Water limited", NIRparents$Treatment.2)
NIRparents$Treatment.2 <- gsub("100", "Well watered", NIRparents$Treatment.2)


p <- ggline(data = NIRparents, x = 'days', y = 'nir_mean', group = 'Treatment.2', color = 'Treatment.2', 
       add = c("mean_se"), palette = c("orange1","cornflowerblue"), facet.by = 'Genotype.ID',
       xlab = "Days", ylab = "NIR reflectance intensity")
ggpar(p, legend.title = "Treatment", legend = "right")

#NIR_sm <- filter(avgNIRparents, days ==25 & Genotype.ID == "Ab")
# ggboxplot(NIR_sm, x = 'Treatment.2', y = 'nir_mean', facet.by = 'Genotype.ID',
#           ylim=  c(20000,26500)) +
#   # rotate_x_text(30)+
#   stat_compare_means(method = 't.test', label = 'p.signif', label.y = 25000, label.x =1.5)


# KW test to compare mean NIR signal on day 15, day 20, and day 26

avgNIRparents <- NIRparents %>%group_by(Genotype.ID, Treatment.2, days, Timepoint, Replicate) %>% summarise(mean_NIR = mean(nir_mean))
avgNIRparents <- ungroup(avgNIRparents)

NIR_sm <- filter(avgNIRparents, days == 26 & Genotype.ID == "Ze")
kruskal.test(mean_NIR ~ Treatment.2, data = NIR_sm)



####################################################
#plot all genotypes pixel area faceted by treatment

plotDF <- unite(avgDF, Genotype.ID, Timepoint, col = label, sep = "_", remove = FALSE)

#cols <- c("ZT1"="blue", "ZT7" = "green", "ZT13" = "pink", "ZT19" = "yellow")
ggplot(data = plotDF, aes(x = days, y = area, group = label)) +
  geom_smooth(aes(color = PH))+
  #geom_point(aes(color = Treatment.2))+
  facet_wrap("Treatment.2")+
  labs(title = "Area SV", x = "days", y = "area")+
  scale_color_brewer(palette="Set2")+
  scale_x_continuous(breaks = seq(0,28, by = 4))+
  #scale_color_manual(name= "Timepoint", breaks = c("ZT1", "ZT7", "ZT13", "ZT19"), labels= c("ZT1", "ZT7", "ZT13", "ZT19"), values = cols)+
  theme(axis.title.x=element_text(hjust = 0.5, face="bold"),
        axis.title.y=element_text(face="bold"),
        strip.text.x = element_text(size = 10, margin = margin(t=1, b=1), color = "black", hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "grey"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        axis.ticks = element_line(colour = "grey", size = 0.2))


######################################################################
# color data with yellow leaves removed using naive bayes classifier #
######################################################################

# read in raw plantcv data and combine

z500 <- read.csv("color_results-SV_z500_single-value-traits.csv")
z500 <- z500 %>% rename("camera" = "frame", "frame" = "camera")
z1500 <- read.csv("color_results_SV_z1500-single-value-traits.csv")
z2500 <- read.csv("color_results_SV_z2500-single-value-traits.csv")

# read in barcode file with genotype and treatment information
barcodes <- read.csv("~/Documents/brassica-image-analysis/barcodes.csv") #read in the barcodes and genotype info
barcodes = select(barcodes, Barcodes, Genotype.ID, Treatment.2, Timepoint, Replicate)

# filter the data files to get just the columns you need
filtered500 <- select(z500, frame, zoom, timestamp, image, plantbarcode,
                      hue_circular_mean, hue_circular_std, hue_median, percent_scenesed)
filtered1500 <- select(z1500, frame, zoom, timestamp, image, plantbarcode,
                       hue_circular_mean, hue_circular_std, hue_median, percent_scenesed)
filtered2500 <- select(z2500, frame, zoom, timestamp, image, plantbarcode,
                       hue_circular_mean, hue_circular_std, hue_median, percent_scenesed)

TS <- c("date","time") #make a vector for the column names of each part of the timestamp
filtered500 <- separate(filtered500, col = timestamp, into = c("date","time"), sep = " ")
filtered1500 <- separate(filtered1500, col = timestamp, into = c("date","time"), sep = " ")
filtered2500 <- separate(filtered2500, col = timestamp, into = c("date","time"), sep = " ") #separate each component of timestamp


#calculate days after first image
startdate <- as.POSIXct("2016-01-13")
filtered500$days=NA
filtered1500$days = NA
filtered2500$days = NA
filtered500$days = as.integer(difftime(filtered500$date, startdate, units = "days"))
filtered1500$days = as.integer(difftime(filtered1500$date, startdate, units = "days"))
filtered2500$days = as.integer(difftime(filtered2500$date, startdate, units = "days"))

# zoom correction
# the zoom correction factors were cacluated from another script and is based on a size marker
filtered500 <- filtered500 %>% mutate_at(c("area", "convex_hull_area"), list(calib = ~ ./1229.355))
filtered500 <- filtered500 %>% mutate_at(c("perimeter", "width", "height", "longest_path", "ellipse_major_axis", "ellipse_minor_axis"), 
                                         list(calib = ~ ./33.95782))
filtered1500 <- filtered1500 %>% mutate_at(c("area", "convex_hull_area"), list(calib = ~ ./2610.580))
filtered1500 <- filtered1500 %>% mutate_at(c("perimeter", "width", "height", "longest_path", "ellipse_major_axis", "ellipse_minor_axis"), 
                                           list(calib = ~ ./49.41139))
filtered2500 <- filtered2500 %>% mutate_at(c("area", "convex_hull_area"), list(calib = ~ ./5543.662))
filtered2500 <- filtered2500 %>% mutate_at(c("perimeter", "width", "height", "longest_path", "ellipse_major_axis", "ellipse_minor_axis"), 
                                           list(calib = ~ ./71.89759))

# separate each filtered DF by imaging day and match barcodes 
DF0 <- filter(filtered2500, filtered2500$days == "0")
DF0 = merge(DF0, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF0 = drop_na(DF0)
DF1 <- filter(filtered2500, filtered2500$days == "1")
DF1 = merge(DF1, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF1 = drop_na(DF1)
DF2 <- filter(filtered2500, filtered2500$days == "2")
DF2 = merge(DF2, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF2 = drop_na(DF2)
DF3 <- filter(filtered2500, filtered2500$days == "3")
DF3 = merge(DF3, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF3 = drop_na(DF3)
DF4 <- filter(filtered2500, filtered2500$days == "4")
DF4 = merge(DF4, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF4 = drop_na(DF4)
DF5 <- filter(filtered2500, filtered2500$days == "5")
DF5 = merge(DF5, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF5 = drop_na(DF5)
DF6 <- filter(filtered2500, filtered2500$days == "6")
DF6 = merge(DF6, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF6 = drop_na(DF6)
DF7 <- filter(filtered2500, filtered2500$days == "7")
DF7 = merge(DF7, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF7 = drop_na(DF7)
DF8 <- filter(filtered2500, filtered2500$days == "8")
DF8 = merge(DF8, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF8 = drop_na(DF8)
DF9 <- filter(filtered2500, filtered2500$days == "9")
DF9 = merge(DF9, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF9 = drop_na(DF9)
DF10 <- filter(filtered2500, filtered2500$days == "10")
DF10 = merge(DF10, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF10 = drop_na(DF10)
DF10.1 <- filter(filtered1500, filtered1500$days == "10")
DF10.1 = merge(DF10, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF10.1 = drop_na(DF10)
DF11 <- filter(filtered1500, filtered1500$days == "11")
DF11 = merge(DF11, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF11 = drop_na(DF11)
DF12 <- filter(filtered1500, filtered1500$days == "12")
DF12 = merge(DF12, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF12 = drop_na(DF12)
DF13 <- filter(filtered1500, filtered1500$days == "13")
DF13 = merge(DF13, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF13 = drop_na(DF13)
DF13.1 <- filter(filtered500, filtered500$days == "13")
DF13.1 = merge(DF13.1, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF13.1 = drop_na(DF13.1)
DF14 <- filter(filtered500, filtered500$days == "14")
DF14 = merge(DF14, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF14 = drop_na(DF14)
DF15 <- filter(filtered500, filtered500$days == "15")
DF15 = merge(DF15, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF15 = drop_na(DF15)
DF16 <- filter(filtered500, filtered500$days == "16")
DF16 = merge(DF16, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF16 = drop_na(DF16)
DF17 <- filter(filtered500, filtered500$days == "17")
DF17 = merge(DF17, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF17 = drop_na(DF17)
DF18 <- filter(filtered500, filtered500$days == "18")
DF18 = merge(DF18, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF18 = drop_na(DF18)
DF19 <- filter(filtered500, filtered500$days == "19")
DF19 = merge(DF19, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF19 = drop_na(DF19)
DF20 <- filter(filtered500, filtered500$days == "20")
DF20 = merge(DF20, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF20 = drop_na(DF20)
DF21 <- filter(filtered500, filtered500$days == "21")
DF21 = merge(DF21, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF21 = drop_na(DF21)
DF22 <- filter(filtered500, filtered500$days == "22")
DF22 = merge(DF22, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF22 = drop_na(DF22)
DF23 <- filter(filtered500, filtered500$days == "23")
DF23 = merge(DF23, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF23 = drop_na(DF23)
DF24 <- filter(filtered500, filtered500$days == "24")
DF24 = merge(DF24, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF24 = drop_na(DF24)
DF25 <- filter(filtered500, filtered500$days == "25")
DF25 = merge(DF25, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF25 = drop_na(DF25)
DF26 <- filter(filtered500, filtered500$days == "26")
DF26 = merge(DF26, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF26 = drop_na(DF26)
DF27 <- filter(filtered500, filtered500$days == "27")
DF27 = merge(DF27, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF27 = drop_na(DF27)
DF28 <- filter(filtered500, filtered500$days == "28")
DF28 = merge(DF28, barcodes, by.x = "plantbarcode", by.y = "Barcodes", all = TRUE)
DF28 = drop_na(DF28)
DF <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF10.1,DF11,DF12,DF13,
                DF13.1,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28)

rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF10.1,DF11,DF12,DF13,
   DF13.1,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28)

# toDelete <- seq(1, nrow(DF), 2)
# DF <- DF[-toDelete,]

write.csv(DF, "~/Documents/brassica-image-analysis/yellow_leaf_naive_bayes/VIS_SV_color_results_without_yellow_leaves.csv", row.names = FALSE)


###################################################
# read in compiled and currated data for plotting #
###################################################
setwd("~/Documents/brassica-image-analysis/yellow_leaf_naive_bayes/")
library(tidyverse)
library(ggpubr)


DF<- read.csv("VIS_SV_color_results_without_yellow_leaves.csv") # read in data

DF$Treatment.2 <- as.factor(DF$Treatment.2) # make sure the treatments are factors

# avgDF <- DF %>% group_by(Genotype.ID, Treatment.2, days, Replicate, Timepoint) %>% 
#   dplyr::summarise(hue_mean = mean(hue_circular_mean), percent_scenesed = mean(percent_scenesed*100)) #add 0 and 90 deg frame values together
# avgDF <- ungroup(avgDF)

# filter DF to have a DF of only the parent genotypes
parents <- filter(DF, Genotype.ID == "Ab4-28"|Genotype.ID == "Al1-29"| Genotype.ID == "Av1-30"|
                    Genotype.ID == "Br2-31"| Genotype.ID == "Ca4-32"| Genotype.ID == "Da2-44"|
                    Genotype.ID == "DH12.6-33"| Genotype.ID == "DH20"| Genotype.ID == "Gr2-35"|
                    Genotype.ID == "Mu4-37"| Genotype.ID == "Ne3-38"| Genotype.ID == "Qu1-39"|
                    Genotype.ID == "Se2-40"| Genotype.ID == "st1-45"|
                    Genotype.ID == "Yu1-42"| Genotype.ID == "Ze2")


####################################################
#plot hue circular mean
#cols <- c("20"="#99CCFF", "100" = "#3366CC")
ggplot(data = parents, aes(x = days, y = percent_scenesed , group = Treatment.2)) +
  geom_smooth(method="loess", aes(color = Treatment.2))+
  facet_wrap("Genotype.ID")+
  labs(title = "hue_circular mean SV", x = "days", y = "mean hue angle")+
  scale_x_continuous(breaks = seq(0,28, by = 4))+
  scale_color_manual(name= "Treatment", breaks = c("20","100"), labels= c("Drought","Control"), values = cal_palette(name = "sierra2", type = "discrete"))+
  theme(axis.title.x=element_text(hjust = 0.5,face="bold"),
        axis.title.y=element_text(face="bold"),
        strip.text.x = element_text(size = 10, margin = margin(t=1, b=1), color = "black", hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "grey"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        axis.ticks = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_text(size = 9))
dev.copy(png, "./SV_nirmean.png", width=1000, height=700, res=120)
dev.off()

plotDF <- filter(parents, days == "25")
ggplot(data = plotDF, aes(x = Treatment.2, y = hue_mean)) +
  geom_boxplot()+
  facet_wrap("Genotype.ID")

ggplot(data = plotDF, aes(x = Treatment.2, y = percent_scenesed)) +
  geom_boxplot()+
  facet_wrap("Genotype.ID")

plotDF <- filter(parents, days >=11 & days <= 25)
p <- ggline(data = plotDF, x = 'days', y = 'hue_circular_mean', group = 'Treatment.2', color = 'Treatment.2', 
            add = c("mean_se"), palette = c("#E7B800","#00AFBB"), facet.by = 'Genotype.ID',
            xlab = "Days", ylab = "Hue circular mean (degrees)")
ggpar(p, legend.title = "Treatment")

small.par <- filter(parents, days >= 20 & days <= 25)
ggline(data =small.par, x = 'days', y = 'percent_scenesed', group = 'Treatment.2', color = 'Treatment.2', 
       add = c("mean_se"), palette = c("#E7B800", "#00AFBB"), facet.by = 'Genotype.ID', 
       xlab = "Days", ylab = "Percentage of leaf area scenesced (%)")

par25 <- filter(parents, days == 25)
ggboxplot(par25, x = 'Treatment.2', y = 'hue_circular_mean', facet.by = 'Genotype.ID',
          ylim=  c(50,150)) +
  # rotate_x_text(30)+
  stat_compare_means(method = 't.test', label = 'p.signif', label.y = 130, label.x = 1.5)

#########
# T-tests of last day

sm_par <- filter(parents, Genotype.ID == "Ab4-28" & days == 25)
x <- sm_par %>% filter(Treatment.2 == "100") %>% select(hue_circular_mean)
y <- sm_par %>% filter(Treatment.2 == "20") %>% select(hue_circular_mean)

t.test(x, y, alternative = "two.sided", var.equal = FALSE)
