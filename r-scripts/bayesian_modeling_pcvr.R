library(tidyverse)
library(pcvr)
library(data.table)
library(patchwork)

setwd("~/Documents/brassica-image-analysis/shape_data/")

DF <- read.csv("./VIS_SV_data_compiled_currated_zoomcorrected.csv")
DF <- select(DF, frame, area_calib, Genotype.ID, Treatment.2, days, Replicate, Timepoint)

DF$Treatment.2 <- as.factor(DF$Treatment.2) # make sure the treatments are factors

avgDF <- dplyr::group_by(DF, Genotype.ID, Treatment.2, days, Replicate, Timepoint) %>% dplyr::summarise(area = sum(area_calib)) #add 0 and 90 deg frame values together
avgDF <- ungroup(avgDF)

# filter DF to have a DF of only the parent genotypes
parents <- filter(avgDF, Genotype.ID == "Ab4-28"|Genotype.ID == "Al1-29"| Genotype.ID == "Av1-30"|
                    Genotype.ID == "Br2-31"| Genotype.ID == "Ca4-32"| Genotype.ID == "Da2-44"|
                    Genotype.ID == "DH12.6-33"| Genotype.ID == "DH20"| Genotype.ID == "Gr2-35"|
                    Genotype.ID == "Mu4-37"| Genotype.ID == "Ne3-38"| Genotype.ID == "Qu1-39"|
                    Genotype.ID == "Se2-40"| Genotype.ID == "st1-45"|
                    Genotype.ID == "Yu1-42"| Genotype.ID == "Ze2")

parents <- parents %>% unite("plant", c(Genotype.ID, Replicate, Timepoint), sep= "_", remove = FALSE)

#parents <- filter(parents, days <= 26)
#colnames(parents)[3] <- "trt"
#parents$trt <- sub("20", "Water limited", parents$trt)
#parents$trt <- sub("100", "Well watered", parents$trt)


sm.parents <- filter(parents, Genotype.ID == "Ze2")
sm.parents <- filter(sm.parents, days <= 26)
colnames(sm.parents)[3] <- "trt"
sm.parents$trt <- sub("20", "Water limited", sm.parents$trt)
sm.parents$trt <- sub("100", "Well watered", sm.parents$trt)

ss<-growthSS(model="gompertz",
             form =  area~days|Genotype.ID/trt,
             sigma="logistic",
             df=sm.parents,
             start = list("A" = 130, "B" = 15, "C" = 0.25, "subA" = 20, "subB" = 10, "subC" = 3),
             type="brms")


Ze_fit <- fitGrowth(ss, iter = 2000,
                      cores = 2, chains = 2, backend = "cmdstanr",
                      control = list(adapt_delta = 0.999, max_treedepth = 20))

save(Ze_fit, file = "./Gompertz_curves/new_pcvr_plots/Ze_fit.RData")
rm(Ze_fit)
#Av1 <-summary(Av1_fit)

###########################

setwd("~/Documents/brassica-image-analysis/shape_data/")

load("./Gompertz_curves/new_pcvr_plots/Ab_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Al_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Av_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Br_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Ca_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Da_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/DH12_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/DH20_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Gr_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Mu_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Ne_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Qu_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Se_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/St_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Yu_fit.RData")
load("./Gompertz_curves/new_pcvr_plots/Ze_fit.RData")

# brmPlot(Ab4_fit, form = area~days|plant/trt, df = ss$df)+ 
#   labs(y=expression("Area"~"(cm"^2~")"))
# 
sum1 <- summary(Ze_fit)
sum1 <- sum1[["fixed"]]
write.csv(sum1, "./Gompertz_curves/new_pcvr_plots/fit_summaries/Ze_fit_summary.csv")

# 
# brmViolin(model = Ab4_fit, params = NULL,
#           hyp="num/denom>1.05", compareX = "Water limited", againstY = "Well Watered",
#           x="trt", facet = "Genotype.ID",
#           returnData=FALSE)

summary(Ab4_fit)

#brms::hypothesis(st1_fit, "A_trtWellwatered/A_trtWaterlimited > 1.05")

#testing that the well watered plants are at least 50% larger than the water limited treatment
testAb <- brms::hypothesis(Ab4_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testAb$Genotype <- "Ab"
testAl <- brms::hypothesis(Al_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testAl$Genotype <- "Al"
testAv <- brms::hypothesis(Av_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testAv$Genotype <- "Av"
testBr <- brms::hypothesis(Br_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testBr$Genotype <- "Br"
testCa <- brms::hypothesis(Ca_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testCa$Genotype <- "Ca"
testDa <- brms::hypothesis(Da_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testDa$Genotype <- "Da"
testDH12 <- brms::hypothesis(DH12_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testDH12$Genotype <- "DH12"
testDH20 <- brms::hypothesis(DH20_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testDH20$Genotype <- "DH20"
testGr <- brms::hypothesis(Gr_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testGr$Genotype <- "Gr"
testMu <- brms::hypothesis(Mu_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testMu$Genotype <- "Mu"
testNe <- brms::hypothesis(Ne_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testNe$Genotype <- "Ne"
testQu <- brms::hypothesis(Qu_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testQu$Genotype <- "Qu"
testSe <- brms::hypothesis(Se_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testSe$Genotype <- "Se"
testSt <- brms::hypothesis(St_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testSt$Genotype <- "St"
testYu <- brms::hypothesis(Yu_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testYu$Genotype <- "Yu"
testZe <- brms::hypothesis(Ze_fit, "A_trtWellwatered > 1.5 * A_trtWaterlimited")$hyp
testZe$Genotype <- "Ze"

testall <- rbind(testAb, testAl, testAv, testBr, testCa, testDa, testDH12, testDH20, 
                   testGr, testMu, testNe, testQu, testSe, testSt, testYu, testZe)

write.csv(testall, "./pcvr_bayesian_test_A_results.csv", row.names = FALSE)

#testing that the well watered plant growth rate is at least 5% larger than the water limited treatment
testAb <- brms::hypothesis(Ab4_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testAb$Genotype <- "Ab"
testAl <- brms::hypothesis(Al_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testAl$Genotype <- "Al"
testAv <- brms::hypothesis(Av_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testAv$Genotype <- "Av"
testBr <- brms::hypothesis(Br_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testBr$Genotype <- "Br"
testCa <- brms::hypothesis(Ca_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testCa$Genotype <- "Ca"
testDa <- brms::hypothesis(Da_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testDa$Genotype <- "Da"
testDH12 <- brms::hypothesis(DH12_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testDH12$Genotype <- "DH12"
testDH20 <- brms::hypothesis(DH20_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testDH20$Genotype <- "DH20"
testGr <- brms::hypothesis(Gr_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testGr$Genotype <- "Gr"
testMu <- brms::hypothesis(Mu_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testMu$Genotype <- "Mu"
testNe <- brms::hypothesis(Ne_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testNe$Genotype <- "Ne"
testQu <- brms::hypothesis(Qu_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testQu$Genotype <- "Qu"
testSe <- brms::hypothesis(Se_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testSe$Genotype <- "Se"
testSt <- brms::hypothesis(St_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testSt$Genotype <- "St"
testYu <- brms::hypothesis(Yu_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testYu$Genotype <- "Yu"
testZe <- brms::hypothesis(Ze_fit, "C_trtWellwatered > 1.05 * C_trtWaterlimited")$hyp
testZe$Genotype <- "Ze"

testall <- rbind(testAb, testAl, testAv, testBr, testCa, testDa, testDH12, testDH20, 
                 testGr, testMu, testNe, testQu, testSe, testSt, testYu, testZe)

write.csv(testall, "./pcvr_bayesian_test_C_results.csv", row.names = FALSE)


#testing that the well watered plant B is at least 5% larger than the water limited treatment
testAb <- brms::hypothesis(Ab4_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testAb$Genotype <- "Ab"
testAl <- brms::hypothesis(Al_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testAl$Genotype <- "Al"
testAv <- brms::hypothesis(Av_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testAv$Genotype <- "Av"
testBr <- brms::hypothesis(Br_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testBr$Genotype <- "Br"
testCa <- brms::hypothesis(Ca_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testCa$Genotype <- "Ca"
testDa <- brms::hypothesis(Da_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testDa$Genotype <- "Da"
testDH12 <- brms::hypothesis(DH12_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testDH12$Genotype <- "DH12"
testDH20 <- brms::hypothesis(DH20_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testDH20$Genotype <- "DH20"
testGr <- brms::hypothesis(Gr_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testGr$Genotype <- "Gr"
testMu <- brms::hypothesis(Mu_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testMu$Genotype <- "Mu"
testNe <- brms::hypothesis(Ne_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testNe$Genotype <- "Ne"
testQu <- brms::hypothesis(Qu_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testQu$Genotype <- "Qu"
testSe <- brms::hypothesis(Se_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testSe$Genotype <- "Se"
testSt <- brms::hypothesis(St_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testSt$Genotype <- "St"
testYu <- brms::hypothesis(Yu_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testYu$Genotype <- "Yu"
testZe <- brms::hypothesis(Ze_fit, "B_trtWellwatered > 1.05 * B_trtWaterlimited")$hyp
testZe$Genotype <- "Ze"

testall <- rbind(testAb, testAl, testAv, testBr, testCa, testDa, testDH12, testDH20, 
                 testGr, testMu, testNe, testQu, testSe, testSt, testYu, testZe)

write.csv(testall, "./pcvr_bayesian_test_B_results.csv", row.names = FALSE)

########
# plots

Ab4 <- growthPlot(Ab4_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Ab", size = 5, fontface = "plain")+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 11)) #axis numbers text size

Al1 <- growthPlot(Al_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Al", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank())

Av1 <- growthPlot(Av_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Av", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

Br2 <- growthPlot(Br_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Br", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

Ca4 <- growthPlot(Ca_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Ca", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank()) #axis numbers text size

Da2 <- growthPlot(Da_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Da", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

DH12 <- growthPlot(DH12_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=10, y=300, label="DH12", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

DH20 <- growthPlot(DH20_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=10, y=300, label="DH20", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

Gr2 <- growthPlot(Gr_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Gr", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank()) #axis numbers text size

Mu4 <- growthPlot(Mu_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Mu", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

Ne3 <- growthPlot(Ne_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Ne", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

Qu1 <- growthPlot(Qu_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Qu", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank()) #axis numbers text size

Se2 <- growthPlot(Se_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Se", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks.x = element_blank(),
                        axis.text = element_text(size = 11)) #axis numbers text size

st1 <- growthPlot(St_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="St", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.text = element_text(size = 11)) #axis numbers text size

Yu1 <- growthPlot(Yu_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Yu", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.text = element_text(size = 11)) #axis numbers text size

Ze2 <- growthPlot(Ze_fit, form = area~days|plant/trt)+
  labs(x = "Days", y = paste("Area", "(cm", "\u00B2", ")", sep = ""))+
  geom_text(x=5, y=300, label="Ze", size = 5)+
  ylim(0, 350)+
  theme_light() + theme(panel.grid = element_line(linewidth = 1, color = "gray92"), #panel background grid line size
                        panel.background = element_rect(fill = "white", color = "black", linetype = 1), #panel background color
                        strip.background = element_rect(fill = "gray92", color = "black", linetype = 1), #panel label background color
                        strip.text = element_text(size = 11, color = "black"), #panel label text size
                        plot.margin = margin(5, 1.5, 5, 1.5),
                        axis.ticks = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.text = element_text(size = 11)) #axis numbers text size

Ab4 + Al1 + Av1 + Br2 + Ca4 + Da2 + DH12 + DH20 + Gr2 + Mu4 + Ne3 + Qu1 + Se2 + st1 + Yu1 + Ze2
