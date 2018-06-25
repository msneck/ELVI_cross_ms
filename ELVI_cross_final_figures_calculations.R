## Date: 6/13/2018
## Author: Tom
## Purpose: prepare final figures and some other miscellany for Michelle's Chapter 2 manuscript submission

# Setup -------------------------------------------------------------------

setwd("C:/Users/tm9/Desktop/git local/ELVI_cross_ms")
library(tidyverse)
library(lme4)
library(bbmle)

invlogit<-function(x){exp(x)/(1+exp(x))} 


# Import data -------------------------------------------------------------

## 1. transmission
trans_cross <- read.csv("Working_Micro_Cross_2_23_17.csv")
trans <- trans_cross %>%
  select(cross_type,gen,mom_E_status,mom_heat,ID,dad_ID,mean_dist,count_num_assay,count_E_plus,cross)%>%
  filter(cross_type=="intra" | cross_type=="inter" | cross_type=="hyb",
         mom_E_status==1,count_num_assay>=1,!is.na(mean_dist)) %>%
  mutate(M_ID = as.factor(ID),
         D_ID = as.factor(dad_ID),
         cross = as.factor(cross),
         cross_ID = interaction(M_ID,D_ID),
         dist_scale = (mean_dist - mean(mean_dist,na.rm=T))/sd(mean_dist,na.rm=T),
         transmission = count_E_plus/count_num_assay,
         trans_color = ifelse(cross_type=="intra","#1b9e77",ifelse(cross_type=="inter","#d95f02","#7570b3")),
         trans_pch = ifelse(cross_type=="intra",15,ifelse(cross_type=="inter",16,17)),
        cross_type = fct_relevel(cross_type,"intra","inter"))
f1_trans <- trans %>% filter(gen=="f1") 
f2_trans <- trans %>% filter(gen=="f2") 

## 2. hyphal density
endo_length2 <- read.csv("endo_length2.csv")
## there is an "na" where there should be "NA" in the seed data
endo_length2$seed_tiller_level[which(endo_length2$seed_tiller_level=="na")]<-NA
hyphae<-endo_length2 %>%
  mutate(dist_scale = (mean_dist - mean(mean_dist,na.rm=T))/sd(mean_dist,na.rm=T),
         M_ID = as.factor(M_ID),
         D_ID = as.factor(D_ID),
         endo = as.factor(endo),
         plant = as.factor(plant),
         cross_ID = interaction(M_ID,D_ID),
         tiller_seeds = as.numeric(seed_tiller_level),
         log_hyphae = log(mean_length),
         plant_level_trans = plant_level_trans,
         hyph_color = ifelse(cross_type=="intra","#1b9e77",ifelse(cross_type=="inter","#d95f02","#7570b3")),
         hyph_pch = ifelse(cross_type=="intra",15,ifelse(cross_type=="inter",16,17)),
         cross_type = fct_relevel(cross_type,"intra","inter"))

## 3. germination
germ <- trans_cross %>%
  select(cross_type,gen,mom_E_status,mom_heat,ID,dad_ID,mean_dist,count_germ,count_num_germ_assay)%>%
  filter(cross_type=="intra" | cross_type=="inter" | cross_type=="hyb",
         mom_E_status==1,count_num_germ_assay>=1,gen=="f1",!is.na(mean_dist)) %>%
  mutate(M_ID = as.factor(ID),
         D_ID = as.factor(dad_ID),
         cross_ID = interaction(M_ID,D_ID),
         dist_scale = (mean_dist - mean(mean_dist,na.rm=T))/sd(mean_dist,na.rm=T),
         mean_germ = count_germ/count_num_germ_assay,
         germ_color = ifelse(cross_type=="intra","#1b9e77",ifelse(cross_type=="inter","#d95f02","#7570b3")),
         germ_pch = ifelse(cross_type=="intra",15,ifelse(cross_type=="inter",16,17)),
         cross_type = fct_relevel(cross_type,"intra","inter"))

## 4. Fertility
cg_data <- read.csv("common_garden_2_22_17.csv")
fertility <- cg_data %>%
  select(e_status,cross_type,indiv_ID,M_ID,D_ID,mean_dist,Rep_Till_Num_16,unbag_count,bag_count,
         cross,infl_collected_16,tot_assay_true)%>%
  filter(cross_type=="intra" | cross_type=="inter" | cross_type=="hyb",
         e_status==1, !is.na(unbag_count)) %>%
  mutate(M_ID = as.factor(M_ID),
         D_ID = as.factor(D_ID),
         cross = as.factor(cross),
         dist_scale = (mean_dist - mean(mean_dist,na.rm=T))/sd(mean_dist,na.rm=T),
         tot_seeds = Rep_Till_Num_16*unbag_count,
         fert_color = ifelse(cross_type=="intra","#1b9e77",ifelse(cross_type=="inter","#d95f02","#7570b3")),
         fert_pch = ifelse(cross_type=="intra",15,ifelse(cross_type=="inter",16,17)),
         cross_type = fct_relevel(cross_type,"intra","inter"))

## 5. Combine greenhouse and common garden data into fitness metric
fit_dat <- read.csv("host_fitness_composite_05_6_18.csv")
f1_dat <- fit_dat %>% 
  filter(gen=="f1",
         mom_E_status==1) %>% 
  mutate(cross_ID = interaction(mom_ID,dad_ID)) %>% 
  select(cross_ID,count_germ,count_num_germ_assay,count_E_plus,count_num_assay,mean_dist,cross_type) %>% 
  mutate(germ = ifelse(count_num_germ_assay>0,count_germ/count_num_germ_assay,NA),
         trans_f1 = ifelse(count_num_assay>0,count_E_plus/count_num_assay,NA),
         dist_scale = (mean_dist - mean(mean_dist,na.rm=T))/sd(mean_dist,na.rm=T))

f2_dat <- fit_dat %>% 
  filter(gen=="f2",
         mom_E_status==1) %>% 
  mutate(cross_ID = interaction(mom_ID,dad_ID)) %>% 
  select(cross_ID,Rep_Till_Num_16,unbag_count,count_E_plus,count_num_assay) %>% 
  mutate(trans_f2 = ifelse(count_num_assay>0,count_E_plus/count_num_assay,NA),
         infs = Rep_Till_Num_16,
         seeds_per_inf = unbag_count)

comp_fit <- f2_dat %>% 
  group_by(cross_ID) %>% 
  summarise(mean_trans_f2 = mean(trans_f2,na.rm=T),
            mean_infs = mean(infs,na.rm=T),
            mean_seeds = mean(seeds_per_inf,na.rm=T)) %>% 
  full_join(f1_dat,.,by="cross_ID") %>% 
  mutate(host_fitness = germ * mean_infs * mean_seeds,
         endo_fitness = ifelse(host_fitness>0, trans_f1 * host_fitness * mean_trans_f2,0),
         fit_color = ifelse(cross_type=="intra","#1b9e77",ifelse(cross_type=="inter","#d95f02","#7570b3")),
         host_pch = ifelse(cross_type=="intra",15,ifelse(cross_type=="inter",16,17)),
         endo_pch = ifelse(cross_type=="intra",0,ifelse(cross_type=="inter",1,2)),
         cross_type = fct_relevel(cross_type,"intra","inter"))

# Check on sample sizes ---------------------------------------------------

## what was the range of seedlings/seeds assyed for transmission estimation?
hist(f1_trans$count_num_assay);range(f1_trans$count_num_assay);mean(f1_trans$count_num_assay);sum(f1_trans$count_num_assay)
hist(f2_trans$count_num_assay);range(f2_trans$count_num_assay);mean(f2_trans$count_num_assay);sum(f2_trans$count_num_assay)

## outliers in fitness data
## these are the three way-out outliers
comp_fit %>% filter(host_fitness > 1000)

# Model selection ---------------------------------------------------------

## F1 model selection
f1_trans_models<-list()
f1_trans_models[[1]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              (1|M_ID) + (1|D_ID), 
                            data=f1_trans, family="binomial")
f1_trans_models[[2]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              dist_scale + (1|M_ID) + (1|D_ID), 
                            data=f1_trans, family="binomial")
f1_trans_models[[3]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              cross_type + (1|M_ID) + (1|D_ID), 
                            data=f1_trans, family="binomial")
f1_trans_models[[4]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              dist_scale + cross_type + (1|M_ID) + (1|D_ID), 
                            data=f1_trans, family="binomial")
f1_trans_models[[5]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              dist_scale * cross_type + (1|M_ID) + (1|D_ID), 
                            data=f1_trans, family="binomial")
f1_trans_AIC<-AICtab(f1_trans_models,weights=T,sort=T)

## F2 model selection
f2_trans_models<-list()
f2_trans_models[[1]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              (1|M_ID) + (1|D_ID) + (1|cross), 
                            data=f2_trans, family="binomial")
f2_trans_models[[2]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              dist_scale + (1|M_ID) + (1|D_ID) + (1|cross), 
                            data=f2_trans, family="binomial")
f2_trans_models[[3]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              cross_type + (1|M_ID) + (1|D_ID) + (1|cross), 
                            data=f2_trans, family="binomial")
f2_trans_models[[4]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              dist_scale + cross_type + (1|M_ID) + (1|D_ID) + (1|cross), 
                            data=f2_trans, family="binomial")
f2_trans_models[[5]]<-glmer(cbind(count_E_plus,count_num_assay-count_E_plus)~ 
                              dist_scale * cross_type + (1|M_ID) + (1|D_ID) + (1|cross), 
                            data=f2_trans, family="binomial")
AICtab(f2_trans_models,weights=T,sort=T)
f2_trans_AICmodavg <- AICtab(f2_trans_models[[1]],f2_trans_models[[2]],base=T,weights=T)
  
## Hyphal density model selection
hyphae_models<-list()
hyphae_models[[1]] <- lmer(log_hyphae ~ (1|M_ID) +(1|D_ID) + (1|plant),REML=F,data=hyphae)
hyphae_models[[2]] <- lmer(log_hyphae ~ dist_scale + (1|M_ID) +(1|D_ID) + (1|plant),REML=F,data=hyphae)
hyphae_models[[3]] <- lmer(log_hyphae ~ cross_type + (1|M_ID) +(1|D_ID) + (1|plant),REML=F,data=hyphae)
hyphae_models[[4]] <- lmer(log_hyphae ~ dist_scale + cross_type + (1|M_ID) +(1|D_ID) + (1|plant),REML=F,data=hyphae)
hyphae_models[[5]] <- lmer(log_hyphae ~ dist_scale * cross_type + (1|M_ID) +(1|D_ID) + (1|plant),REML=F,data=hyphae)
AICtab(hyphae_models,weights=T,sort=T)
hyphae_AICmodavg <- AICtab(hyphae_models[[1]],hyphae_models[[2]],
                             hyphae_models[[3]],hyphae_models[[4]],base=T,sort=F,weights=T)


## Germination model selection
germ_models<-list()
germ_models[[1]] <- glmer(cbind(count_germ,count_num_germ_assay-count_germ)~ (1|M_ID) + (1|D_ID), data=germ, family="binomial")
germ_models[[2]] <- glmer(cbind(count_germ,count_num_germ_assay-count_germ)~ dist_scale + (1|M_ID) + (1|D_ID), data=germ, family="binomial")
germ_models[[3]] <- glmer(cbind(count_germ,count_num_germ_assay-count_germ)~ cross_type + (1|M_ID) + (1|D_ID), data=germ, family="binomial")
germ_models[[4]] <- glmer(cbind(count_germ,count_num_germ_assay-count_germ)~ dist_scale + cross_type + (1|M_ID) + (1|D_ID), data=germ, family="binomial")
germ_models[[5]] <- glmer(cbind(count_germ,count_num_germ_assay-count_germ)~ dist_scale * cross_type + (1|M_ID) + (1|D_ID), data=germ, family="binomial")
AICtab(germ_models,weights=T,sort=T)

## Repro tiller model selection
Rtiller_models<-list()
Rtiller_models[[1]] <- glmer((Rep_Till_Num_16) ~ (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
Rtiller_models[[2]] <- glmer((Rep_Till_Num_16) ~ dist_scale + (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
Rtiller_models[[3]] <- glmer((Rep_Till_Num_16) ~ cross_type + (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
Rtiller_models[[4]] <- glmer((Rep_Till_Num_16) ~ dist_scale + cross_type + (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
Rtiller_models[[5]] <- glmer((Rep_Till_Num_16) ~ dist_scale * cross_type + (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
AICtab(Rtiller_models,weights=T,sort=T)
Rtiller_AICmodavg <- AICtab(Rtiller_models[[2]],Rtiller_models[[4]],
                            Rtiller_models[[5]],base=T,sort=F,weights=T)

  
## Seed production model selection
seedcount_models<-list()
seedcount_models[[1]] <- glmer(unbag_count ~ (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
seedcount_models[[2]] <- glmer(unbag_count ~ dist_scale + (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
seedcount_models[[3]] <- glmer(unbag_count ~ cross_type + (1|M_ID) + (1|D_ID),data=fertility,family="poisson")
seedcount_models[[4]] <- glmer(unbag_count ~ dist_scale + cross_type + (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
seedcount_models[[5]] <- glmer(unbag_count ~ dist_scale * cross_type + (1|M_ID) + (1|D_ID) + (1|cross),data=fertility,family="poisson")
AICtab(seedcount_models,weights=T,sort=T)
seed_AICmodavg <- AICtab(seedcount_models[[2]],seedcount_models[[4]],
                            base=T,sort=F,weights=T)

# Figure 1 ----------------------------------------------------------------
## F1 and F2 transmission
coef_trans<-fixef(f1_trans_models[[5]])
f1trans_intra_pred <- tibble(cross_type = "intra",
                             x = seq(range(trans$dist_scale[trans$cross_type=="intra"],na.rm=T)[1],
                                     range(trans$dist_scale[trans$cross_type=="intra"],na.rm=T)[2],0.1),
                             predicted=invlogit(coef_trans[1]+(coef_trans[2])*x),
                             trans_color = unique(trans$trans_color[trans$cross_type=="intra"]))
f1trans_inter_pred <- tibble(cross_type = "inter",
                             x = seq(range(trans$dist_scale[trans$cross_type=="inter"],na.rm=T)[1],
                                     range(trans$dist_scale[trans$cross_type=="inter"],na.rm=T)[2],0.1),
                             predicted=invlogit(coef_trans[1]+coef_trans[3]+(coef_trans[2]+coef_trans[5])*x),
                             trans_color = unique(trans$trans_color[trans$cross_type=="inter"]))
f1trans_hyb_pred <- tibble(cross_type = "hyb",
                           x = seq(range(trans$dist_scale[trans$cross_type=="hyb"],na.rm=T)[1],
                                   range(trans$dist_scale[trans$cross_type=="hyb"],na.rm=T)[2],0.1),
                           predicted=invlogit(coef_trans[1]+coef_trans[4]+(coef_trans[2]+coef_trans[6])*x),
                           trans_color = unique(trans$trans_color[trans$cross_type=="hyb"]))
f1trans_pred <- bind_rows(f1trans_intra_pred,f1trans_inter_pred,f1trans_hyb_pred)%>%
  mutate(gen="f1")

## model average for F2 transmission
modavg_coef_trans_f2<-c()
modavg_coef_trans_f2[1]<-f2_trans_AICmodavg$weight[1]*fixef(f2_trans_models[[1]])[1] +
  f2_trans_AICmodavg$weight[2]*fixef(f2_trans_models[[2]])[1]
modavg_coef_trans_f2[2]<-f2_trans_AICmodavg$weight[1]*0 +
  f2_trans_AICmodavg$weight[2]*fixef(f2_trans_models[[2]])[2]
f2trans_pred <- tibble(x = seq(range(trans$dist_scale,na.rm=T)[1],
                                     range(trans$dist_scale,na.rm=T)[2],0.1),
                             predicted=invlogit(modavg_coef_trans_f2[1]+modavg_coef_trans_f2[2]*x))

trans_means <- trans %>% 
  group_by(gen,cross_type) %>% 
  summarize(mean = mean(transmission),
            SE = sd(transmission)/sqrt(n())) %>% 
  left_join(.,
            unique(select(trans,gen,cross_type,trans_color)),
            by=c("gen","cross_type"))

max_seeds_scored <- max(trans$count_num_assay)
x_range <- range(trans$dist_scale)

win.graph()
test <- layout(matrix(c(1,2,3,4),2,2,byrow=T),
               heights = c(3,3,3,3),widths = c(2,1.25))
layout.show(test)
par(mar=c(4.5,4.5,1.5,1.5))
plot(f1_trans$dist_scale,f1_trans$transmission,col=alpha(f1_trans$trans_color,0.4),
     xlab="Genetic distance",ylab="Vertical transmission",
     cex=(f1_trans$count_num_assay/max_seeds_scored)*5,
     lwd=2,ylim=c(0,1),pch=f1_trans$trans_pch,xlim=x_range,cex.lab=1.4)
title("A",adj=0)
lines(f1trans_intra_pred$x,
      f1trans_intra_pred$predicted,
      col=f1trans_intra_pred$trans_color,
      lwd=3)
lines(f1trans_inter_pred$x,
      f1trans_inter_pred$predicted,
      col=f1trans_inter_pred$trans_color,
      lwd=3)
lines(f1trans_hyb_pred$x,
      f1trans_hyb_pred$predicted,
      col=f1trans_hyb_pred$trans_color,
      lwd=3)

trans_f1_bar <- barplot(trans_means$mean[trans_means$gen=="f1"],col=trans_means$trans_color[trans_means$gen=="f1"],
                        ylim=c(0,1),names.arg=trans_means$cross_type[1:3],cex.names = 1,
                        xlab="Cross type",cex.lab=1.4,ylab="Vertical transmission");box()
title("B",adj=0)
arrows(trans_f1_bar,trans_means$mean[trans_means$gen=="f1"]-2*trans_means$SE[trans_means$gen=="f1"],
       trans_f1_bar,trans_means$mean[trans_means$gen=="f1"]+2*trans_means$SE[trans_means$gen=="f1"],
       code=0)

plot(f2_trans$dist_scale,f2_trans$transmission,col=alpha(f2_trans$trans_color,0.4),
     xlab="Genetic distance",ylab="Vertical transmission",
     cex=(f2_trans$count_num_assay/max_seeds_scored)*5,
     lwd=2,ylim=c(0,1),pch=f1_trans$trans_pch,xlim=x_range,cex.lab=1.4)
title("C",adj=0)
lines(f2trans_pred$x,
      f2trans_pred$predicted,
      lwd=3)

trans_f2_bar <- barplot(trans_means$mean[trans_means$gen=="f2"],col=trans_means$trans_color[trans_means$gen=="f2"],
                        ylim=c(0,1),names.arg=trans_means$cross_type[1:3],
                        xlab="Cross type",cex.lab=1.4,ylab="Vertical transmission");box()
title("D",adj=0)
arrows(trans_f2_bar,trans_means$mean[trans_means$gen=="f2"]-2*trans_means$SE[trans_means$gen=="f2"],
       trans_f2_bar,trans_means$mean[trans_means$gen=="f2"]+2*trans_means$SE[trans_means$gen=="f2"],
       code=0)

# Figure 2 ----------------------------------------------------------------
## model averaging
modavg_coef_hyphae<-c()
modavg_coef_hyphae[1]<-hyphae_AICmodavg$weight[1]*fixef(hyphae_models[[1]])[1]+
  hyphae_AICmodavg$weight[2]*fixef(hyphae_models[[2]])[1]+
  hyphae_AICmodavg$weight[3]*fixef(hyphae_models[[3]])[1]+
  hyphae_AICmodavg$weight[4]*fixef(hyphae_models[[4]])[1]
modavg_coef_hyphae[2]<-hyphae_AICmodavg$weight[1]*0+
  hyphae_AICmodavg$weight[2]*fixef(hyphae_models[[2]])[2]+
  hyphae_AICmodavg$weight[3]*0+
  hyphae_AICmodavg$weight[4]*fixef(hyphae_models[[4]])[2]
modavg_coef_hyphae[3]<-hyphae_AICmodavg$weight[1]*0+
  hyphae_AICmodavg$weight[2]*0+
  hyphae_AICmodavg$weight[3]*fixef(hyphae_models[[3]])[2]+
  hyphae_AICmodavg$weight[4]*fixef(hyphae_models[[4]])[3]  
  
hyphae_inter_pred <- tibble(cross_type = "inter",
                             x = seq(range(hyphae$dist_scale[hyphae$cross_type=="inter"],na.rm=T)[1],
                                     range(hyphae$dist_scale[hyphae$cross_type=="inter"],na.rm=T)[2],0.1),
                             predicted=(modavg_coef_hyphae[1]+modavg_coef_hyphae[2]*x),
                      hyph_color = unique(hyphae$hyph_color[hyphae$cross_type=="inter"]))
hyphae_hyb_pred <- tibble(cross_type = "hyb",
                           x = seq(range(hyphae$dist_scale[hyphae$cross_type=="hyb"],na.rm=T)[1],
                                   range(hyphae$dist_scale[hyphae$cross_type=="hyb"],na.rm=T)[2],0.1),
                           predicted=(modavg_coef_hyphae[1]+modavg_coef_hyphae[3]+modavg_coef_hyphae[2]*x),
                           hyph_color = unique(hyphae$hyph_color[hyphae$cross_type=="hyb"]))
hyph_pred <- bind_rows(hyphae_inter_pred,hyphae_hyb_pred)%>%
  mutate(gen="f1")

## means table
hyph_means <- hyphae %>% 
  group_by(cross_type) %>% 
  summarize(mean = mean(log_hyphae),
            SE = sd(log_hyphae)/sqrt(n())) %>% 
  left_join(.,
            unique(select(hyphae,cross_type,hyph_color)),
            by=c("cross_type"))

win.graph()
test <- layout(matrix(c(1,2),1,2,byrow=T),
               heights = c(3,3),widths = c(2,1.25))
layout.show(test)
par(mar=c(4.5,4.5,1.5,1.5))
plot(hyphae$dist_scale,hyphae$log_hyphae,col=alpha(hyphae$hyph_color,0.4),cex=1.4,
     xlab="Genetic distance",ylab=expression(paste("log Hyphal density ( ",mu,"m)")),
     ylim=range(hyphae$log_hyphae),lwd=2,pch=hyphae$hyph_pch,xlim=x_range,cex.lab=1.2)
title("A",adj=0)
lines(hyphae_inter_pred$x,
      hyphae_inter_pred$predicted,
      col=hyphae_inter_pred$hyph_color,
      lwd=3)
lines(hyphae_hyb_pred$x,
      hyphae_hyb_pred$predicted,
      col=hyphae_hyb_pred$hyph_color,
      lwd=3)

hyph_bar <- barplot(hyph_means$mean,col=hyph_means$hyph_color,
                        names.arg=hyph_means$cross_type[1:2],cex.names = 1,
                        xlab="Cross type",cex.lab=1.2,ylim=range(hyphae$log_hyphae),
                    ylab=expression(paste("log Hyphal density ( ",mu,"m)")),xpd=F);box()
title("B",adj=0)
arrows(hyph_bar,hyph_means$mean-2*hyph_means$SE,
       hyph_bar,hyph_means$mean+2*hyph_means$SE,
       code=0)


# Figure 3 ----------------------------------------------------------------
germ_means <- germ %>% 
  group_by(cross_type) %>% 
  summarize(mean = mean(mean_germ),
            SE = sd(mean_germ)/sqrt(n())) %>% 
  left_join(.,
            unique(select(germ,cross_type,germ_color)),
            by=c("cross_type"))

## model predictions
coef_germ<-fixef(germ_models[[5]])
germ_intra_pred <- tibble(cross_type = "intra",
                             x = seq(range(germ$dist_scale[germ$cross_type=="intra"],na.rm=T)[1],
                                     range(germ$dist_scale[germ$cross_type=="intra"],na.rm=T)[2],0.1),
                             predicted=invlogit(coef_germ[1]+(coef_germ[2])*x),
                             germ_color = unique(germ$germ_color[germ$cross_type=="intra"]))
germ_inter_pred <- tibble(cross_type = "inter",
                             x = seq(range(germ$dist_scale[germ$cross_type=="inter"],na.rm=T)[1],
                                     range(germ$dist_scale[germ$cross_type=="inter"],na.rm=T)[2],0.1),
                             predicted=invlogit(coef_germ[1]+coef_germ[3]+(coef_germ[2]+coef_germ[5])*x),
                          germ_color = unique(germ$germ_color[germ$cross_type=="inter"]))
germ_hyb_pred <- tibble(cross_type = "hyb",
                           x = seq(range(germ$dist_scale[germ$cross_type=="hyb"],na.rm=T)[1],
                                   range(germ$dist_scale[germ$cross_type=="hyb"],na.rm=T)[2],0.1),
                           predicted=invlogit(coef_germ[1]+coef_germ[4]+(coef_germ[2]+coef_germ[6])*x),
                        germ_color = unique(germ$germ_color[germ$cross_type=="hyb"]))
germ_pred <- bind_rows(germ_intra_pred,germ_inter_pred,germ_hyb_pred)
max_germ_assay <-max(germ$count_num_germ_assay)
x_range_germ <- range(germ$dist_scale)

win.graph()
test <- layout(matrix(c(1,2),1,2,byrow=T),
               heights = c(3,3),widths = c(2,1.25))
layout.show(test)
par(mar=c(4.5,4.5,1.5,1.5))

plot(germ$dist_scale,germ$mean_germ,col=alpha(germ$germ_color,0.4),
     xlab="Genetic distance",ylab="Germination",
     cex=(germ$count_num_germ_assay/max_germ_assay)*5,
     lwd=2,ylim=c(0,1),pch=germ$germ_pch,xlim=x_range_germ,cex.lab=1.4)
title("A",adj=0)
lines(germ_intra_pred$x,
      germ_intra_pred$predicted,
      col=germ_intra_pred$germ_color,
      lwd=3)
lines(germ_inter_pred$x,
      germ_inter_pred$predicted,
      col=germ_inter_pred$germ_color,
      lwd=3)
lines(germ_hyb_pred$x,
      germ_hyb_pred$predicted,
      col=germ_hyb_pred$germ_color,
      lwd=3)

germ_bar <- barplot(germ_means$mean,col=germ_means$germ_color,
                        ylim=c(0,1),names.arg=germ_means$cross_type[1:3],cex.names = 1,
                        xlab="Cross type",cex.lab=1.4,ylab="Germination");box()
title("B",adj=0)
arrows(germ_bar,germ_means$mean-2*germ_means$SE,
       germ_bar,germ_means$mean+2*germ_means$SE,
       code=0)


# Figure 4 ----------------------------------------------------------------

Rtill_means <- fertility %>% 
  filter(!is.na(Rep_Till_Num_16)) %>% 
  group_by(cross_type) %>% 
  summarize(mean = mean(Rep_Till_Num_16),
            SD = sd(Rep_Till_Num_16),
            SE = sd(Rep_Till_Num_16)/sqrt(n()))%>% 
  left_join(.,
            unique(select(fertility,cross_type,fert_color,fert_pch)),
            by=c("cross_type"))

seed_means <- fertility %>% 
  group_by(cross_type) %>% 
  summarize(mean = mean(unbag_count),
            SE = sd(unbag_count)/sqrt(n())) %>% 
  left_join(.,
            unique(select(fertility,cross_type,fert_color,fert_pch)),
            by=c("cross_type"))

## model average predictions for R tillers
modavg_coef_Rtill<-c()
modavg_coef_Rtill[1] <- Rtiller_AICmodavg$weight[1]*fixef(Rtiller_models[[2]])[1]+
  Rtiller_AICmodavg$weight[2]*fixef(Rtiller_models[[4]])[1]+
  Rtiller_AICmodavg$weight[3]*fixef(Rtiller_models[[5]])[1]
modavg_coef_Rtill[2] <- Rtiller_AICmodavg$weight[1]*fixef(Rtiller_models[[2]])[2]+
  Rtiller_AICmodavg$weight[2]*fixef(Rtiller_models[[4]])[2]+
  Rtiller_AICmodavg$weight[3]*fixef(Rtiller_models[[5]])[2]
modavg_coef_Rtill[3] <- Rtiller_AICmodavg$weight[1]*0+
  Rtiller_AICmodavg$weight[2]*fixef(Rtiller_models[[4]])[3]+
  Rtiller_AICmodavg$weight[3]*fixef(Rtiller_models[[5]])[3]
modavg_coef_Rtill[4] <- Rtiller_AICmodavg$weight[1]*0+
  Rtiller_AICmodavg$weight[2]*fixef(Rtiller_models[[4]])[4]+
  Rtiller_AICmodavg$weight[3]*fixef(Rtiller_models[[5]])[4]  
modavg_coef_Rtill[5] <- Rtiller_AICmodavg$weight[1]*0+
  Rtiller_AICmodavg$weight[2]*0+
  Rtiller_AICmodavg$weight[3]*fixef(Rtiller_models[[5]])[5]
modavg_coef_Rtill[6] <- Rtiller_AICmodavg$weight[1]*0+
  Rtiller_AICmodavg$weight[2]*0+
  Rtiller_AICmodavg$weight[3]*fixef(Rtiller_models[[5]])[6]

coef_Rtill<-modavg_coef_Rtill
Rtill_intra_pred <- tibble(cross_type = "intra",
                          x = seq(range(fertility$dist_scale[fertility$cross_type=="intra"],na.rm=T)[1],
                                  range(fertility$dist_scale[fertility$cross_type=="intra"],na.rm=T)[2],0.1),
                          predicted=exp(coef_Rtill[1]+(coef_Rtill[2])*x),
                          fert_color = unique(fertility$fert_color[fertility$cross_type=="intra"]))
Rtill_inter_pred <- tibble(cross_type = "inter",
                          x = seq(range(fertility$dist_scale[fertility$cross_type=="inter"],na.rm=T)[1],
                                  range(fertility$dist_scale[fertility$cross_type=="inter"],na.rm=T)[2],0.1),
                          predicted=exp(coef_Rtill[1]+coef_Rtill[3]+(coef_Rtill[2] + coef_Rtill[5])*x),
                          fert_color = unique(fertility$fert_color[fertility$cross_type=="inter"]))
Rtill_hyb_pred <- tibble(cross_type = "hyb",
                        x = seq(range(fertility$dist_scale[fertility$cross_type=="hyb"],na.rm=T)[1],
                                range(fertility$dist_scale[fertility$cross_type=="hyb"],na.rm=T)[2],0.1),
                        predicted=exp(coef_Rtill[1]+coef_Rtill[4]+(coef_Rtill[2] + coef_Rtill[6])*x),
                        fert_color = unique(fertility$fert_color[fertility$cross_type=="hyb"]))
Rtill_pred <- bind_rows(Rtill_intra_pred,Rtill_inter_pred,Rtill_hyb_pred)

## model average predictions for seeds
modavg_coef_seeds <- c()
modavg_coef_seeds[1] <- seed_AICmodavg$weight[1]*fixef(seedcount_models[[2]])[1]+
  seed_AICmodavg$weight[2]*fixef(seedcount_models[[4]])[1]
modavg_coef_seeds[2] <- seed_AICmodavg$weight[1]*fixef(seedcount_models[[2]])[2]+
  seed_AICmodavg$weight[2]*fixef(seedcount_models[[4]])[2]
modavg_coef_seeds[3] <- seed_AICmodavg$weight[1]*0+
  seed_AICmodavg$weight[2]*fixef(seedcount_models[[4]])[3]
modavg_coef_seeds[4] <- seed_AICmodavg$weight[1]*0+
  seed_AICmodavg$weight[2]*fixef(seedcount_models[[4]])[4]

coef_seeds<-modavg_coef_seeds
seeds_intra_pred <- tibble(cross_type = "intra",
                           x = seq(range(fertility$dist_scale[fertility$cross_type=="intra"],na.rm=T)[1],
                                   range(fertility$dist_scale[fertility$cross_type=="intra"],na.rm=T)[2],0.1),
                           predicted=exp(coef_seeds[1]+(coef_seeds[2])*x),
                           fert_color = unique(fertility$fert_color[fertility$cross_type=="intra"]))
seeds_inter_pred <- tibble(cross_type = "inter",
                           x = seq(range(fertility$dist_scale[fertility$cross_type=="inter"],na.rm=T)[1],
                                   range(fertility$dist_scale[fertility$cross_type=="inter"],na.rm=T)[2],0.1),
                           predicted=exp(coef_seeds[1]+coef_seeds[3]+(coef_seeds[2])*x),
                           fert_color = unique(fertility$fert_color[fertility$cross_type=="inter"]))
seeds_hyb_pred <- tibble(cross_type = "hyb",
                         x = seq(range(fertility$dist_scale[fertility$cross_type=="hyb"],na.rm=T)[1],
                                 range(fertility$dist_scale[fertility$cross_type=="hyb"],na.rm=T)[2],0.1),
                         predicted=exp(coef_seeds[1]+coef_seeds[4]+(coef_seeds[2])*x),
                         fert_color = unique(fertility$fert_color[fertility$cross_type=="hyb"]))


x_range_fert <- range(fertility$dist_scale,na.rm = T)
y_range_fert <- range(fertility$Rep_Till_Num_16,na.rm = T)
y_range_seeds <- range(fertility$unbag_count,na.rm = T)

win.graph()
test <- layout(matrix(c(1,2,3,4),2,2,byrow=T),
               heights = c(3,3,3,3),widths = c(2,1.25))
layout.show(test)
par(mar=c(4.5,4.5,1.5,1.5))

plot(fertility$dist_scale,fertility$Rep_Till_Num_16,col=alpha(fertility$fert_color,0.4),
     xlab="Genetic distance",ylab="Inflorescences",cex=1.4,
     lwd=2,pch=fertility$fert_pch,xlim=x_range_fert,ylim=y_range_fert,cex.lab=1.4)
title("A",adj=0)
lines(Rtill_intra_pred$x,
      Rtill_intra_pred$predicted,
      col=Rtill_intra_pred$fert_color,
      lwd=3)
lines(Rtill_inter_pred$x,
      Rtill_inter_pred$predicted,
      col=Rtill_inter_pred$fert_color,
      lwd=3)
lines(Rtill_hyb_pred$x,
      Rtill_hyb_pred$predicted,
      col=Rtill_hyb_pred$fert_color,
      lwd=3)

Rtill_bar <- barplot(Rtill_means$mean,col=Rtill_means$fert_color,
                    ylim=y_range_fert,names.arg=Rtill_means$cross_type[1:3],cex.names = 1,
                    xlab="Cross type",cex.lab=1.4,ylab="Inflorescences");box()
title("B",adj=0)
arrows(Rtill_bar,Rtill_means$mean-2*Rtill_means$SE,
       Rtill_bar,Rtill_means$mean+2*Rtill_means$SE,
       code=0)


plot(fertility$dist_scale,fertility$unbag_count,col=alpha(fertility$fert_color,0.4),
     xlab="Genetic distance",ylab="Seeds / inflorescence",cex=1.4,
     lwd=2,pch=fertility$fert_pch,xlim=x_range_fert,ylim=y_range_seeds,cex.lab=1.4)
title("C",adj=0)
lines(seeds_intra_pred$x,
      seeds_intra_pred$predicted,
      col=seeds_intra_pred$fert_color,
      lwd=3)
lines(seeds_inter_pred$x,
      seeds_inter_pred$predicted,
      col=seeds_inter_pred$fert_color,
      lwd=3)
lines(seeds_hyb_pred$x,
      seeds_hyb_pred$predicted,
      col=seeds_hyb_pred$fert_color,
      lwd=3)

seeds_bar <- barplot(seed_means$mean,col=seed_means$fert_color,
                     ylim=y_range_seeds,names.arg=seed_means$cross_type[1:3],cex.names = 1,
                     xlab="Cross type",cex.lab=1.4,ylab="Seeds / inflorescence");box()
title("D",adj=0)
arrows(seeds_bar,seed_means$mean-2*seed_means$SE,
       seeds_bar,seed_means$mean+2*seed_means$SE,
       code=0)


# Figure 5 ----------------------------------------------------------------
comp_fit_means <- comp_fit %>% 
  select(host_fitness,endo_fitness,cross_type) %>% 
  group_by(cross_type) %>% 
  summarise(host = mean(host_fitness,na.rm=T),
            endo = mean(endo_fitness,na.rm=T),
            host_se = sd(host_fitness,na.rm=T)/sqrt(n()),
            endo_se = sd(endo_fitness,na.rm=T)/sqrt(n()))%>% 
  left_join(.,
            unique(select(comp_fit,cross_type,fit_color,host_pch,endo_pch)),
            by=c("cross_type"))

x_range_fit <- range(comp_fit$dist_scale,na.rm = T)
y_range_fit <- range(comp_fit$host_fitness,na.rm = T)


win.graph()
test <- layout(matrix(c(1,2),1,2,byrow=T),
               heights = c(3,3),widths = c(2,1.25))
layout.show(test)
par(mar=c(4.5,4.5,1.5,1.5))
plot(comp_fit$dist_scale,comp_fit$host_fitness,col=comp_fit$fit_color,
     xlab="Genetic distance",ylab="Annual fitness",cex=1.4,
     lwd=2,pch=comp_fit$host_pch,xlim=x_range_fit,ylim=y_range_fit,cex.lab=1.4)
points(comp_fit$dist_scale+0.15,comp_fit$endo_fitness,col=comp_fit$fit_color,
     pch=comp_fit$endo_pch,cex=1.3)
title("A",adj=0)

fit_bar <- barplot(as.matrix(t(comp_fit_means[,c("host","endo")])),beside=T,
                   border=rep(comp_fit_means$fit_color,each=2),cex.names = 1,
                   col=c(comp_fit_means$fit_color[1],"white",
                         comp_fit_means$fit_color[2],"white",
                         comp_fit_means$fit_color[3],"white"),
                  ylim=c(0,600),names.arg=comp_fit_means$cross_type[1:3],
                  xlab="Cross type",cex.lab=1.4,ylab="Annual fitness");box()
title("B",adj=0)
arrows(fit_bar[1,],comp_fit_means$host-2*comp_fit_means$host_se,
       fit_bar[1,],comp_fit_means$host+2*comp_fit_means$host_se,
       code=0)
arrows(fit_bar[2,],comp_fit_means$endo-2*comp_fit_means$endo_se,
       fit_bar[2,],comp_fit_means$endo+2*comp_fit_means$endo_se,
       code=0)


# Drop fitness outliers ---------------------------------------------------

drop_outliers <- comp_fit %>% filter(host_fitness < 1000) %>% 
  select(host_fitness,endo_fitness,cross_type) %>% 
  group_by(cross_type) %>% 
  summarise(host = mean(host_fitness,na.rm=T),
            endo = mean(endo_fitness,na.rm=T),
            host_se = sd(host_fitness,na.rm=T)/sqrt(n()),
            endo_se = sd(endo_fitness,na.rm=T)/sqrt(n()))
  

# Zero fitness stats ------------------------------------------------------

zero_fit_stats <- comp_fit %>% 
  mutate(zero_fit = host_fitness==0) %>%
  group_by(cross_type) %>% 
  summarize(zero = mean(zero_fit,na.rm=T))
