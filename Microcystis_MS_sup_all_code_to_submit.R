#######################################################################################
# Setup -------------------------------------------------------------------
setwd("C:/Users/mirte/Documents/UCSD/Projects/Microcystis Warming")

rm(list = ls())

library("dplyr") #For data manipulation
library("reshape2") #for data manipulation

library("nlme") #For creating linear models
library("effectsize") #For calculating the effect sizes of linear mixed effect models
library("emmeans") #For PostHoc tests form linear mixed effect models
library("fBasics") #Needed for: normalTest to test whether expression data follows a normal distribution
library("MASS") #Needed for: Box-Cox Transformations for Linear Models AND postHoc of multivariate analysis
library("mgcv") #Needed for GAMs
library("gratia") #For plotting GAMs

library("ggplot2") #For plotting
library("ggnewscale") #To allow two colour scales in Figure 2
library("ggpubr") #For ggdensity plotting
library("patchwork") #For arranging facets of Figure 4, S3, S4 & S5
library("scales") #For better axis labeling in Figure S4
library("ggrepel") #For labels in Figure 5
library("ggtree") #For maps

options(es.use_symbols = TRUE) # get nice symbols when printing! 

#######################################################################################
# Figure 1 ----------------------------------------------------------------
tree<-read.tree("RAxML_bipartitions.Testing_Oct_4_2024_noCR1901_10kboots.txt")#,node.label = 'support')
#alignment length 490 bp
metadata_full<-read.csv("Microcystis_mapping.csv")
dd_full<-data.frame(Strain=metadata_full$Strain,TreeGroups=metadata_full$Tree_Groups,datasets=metadata_full$dataset)

#tree_rest<-ggtree(tree)+geom_tiplab(align=TRUE)
full_tree_p<-ggtree(tree) + geom_tree() +  theme_tree() +  geom_tiplab(size=7,hjust = -0.25)#geom_nodelab(geom='label', aes(label=support, subset=support > 80))

dd_full$TreeGroups<-as.factor(dd_full$TreeGroups)
dd_full$datasets<-as.factor(dd_full$datasets)
dd_full$newlabel<-label_pad(dd_full$Strain)
my_full_tree <- full_tree_p %<+% dd_full + geom_tippoint(aes(label=newlabel,color=TreeGroups,alpha=0.5),size=4)+ geom_treescale()+ theme(legend.position="right")+scale_color_manual(values = c("nutrient-rich"="green", "oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","none"="gray"))+guides(color=FALSE,alpha=FALSE)+geom_text2(aes(label=label, size=2, hjust=1.5, vjust= -0.75,subset = !is.na(as.numeric(label)) & as.numeric(label) > 0))#+xlim(NA,0.15)
my_full_tree
head(dd_full)
test<-rotate(my_full_tree)

#keep bootstraps above 50 like in Molecular Ecology paper 
Fig_1 <- full_tree_p %<+% dd_full + 
  geom_tippoint(aes(x=x+0.0008,color=TreeGroups,fill=datasets,hjust=-0.75),size=5,shape=21,stroke=4)+ 
  # geom_nodelab(geom='label', aes(label=support, subset=support > 50))+
  geom_text2(aes(label=label,x=branch,hjust=1.3, vjust= -0.95,subset = !is.na(as.numeric(label)) & as.numeric(label) >50),size=7)+
  theme(legend.position="right")+
  geom_treescale(fontsize=7)+
  scale_color_manual(name="Genotype",labels=c("none"="N/A","nutrient-rich"="HN/HG","pseudo-oligotrophic"="HN/LG","oligotrophic"="LN/LG","LN/HG"='LN/HG'),values = c("LN/HG"="black","nutrient-rich"="green", "oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","none"="gray"))+
  guides(alpha=FALSE,color=FALSE,size=FALSE,fill=FALSE)+
  scale_fill_manual(name="Dataset",labels=c("old"="Original","new"="New"),values=c("old"="black","new"="white","none"="gray"),
                    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)) 

Fig_1


#######################################################################################
# Figure 2 ----------------------------------------------------------------
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# Removing CR19-01, which was the sole LN/HG
data_exp<-data[(data$Phase=="Exp"),]#Sorting for exponential phase growth

#Check data for normality
hist(data_exp$Growth)
shapiro.test(data_exp$Growth)
normalTest(data_exp$Growth, method = "da")

#boxcox transformation of data for better normality
b <- boxcox(lm(data_exp$Growth +1 ~ 1))
lambda <- b$x[which.max(b$y)]
lambda

#check effect of boxcox transformation on data 
hist(((data_exp$Growth + 1)^lambda - 1)/lambda)
ggdensity(((data_exp$Growth + 1)^lambda - 1)/lambda, y = "density", fill = "grey") 
shapiro.test(((data_exp$Growth + 1)^lambda - 1)/lambda)
normalTest(((data_exp$Growth + 1)^lambda - 1)/lambda, method = "da")

#normality of data is improved (though not perfect), so use this transformation
#Transform data
data_exp <- data_exp %>%
  mutate(Growth_transformed = ((Growth + 1)^lambda - 1)/lambda) 

full<-nlme::lme(Growth_transformed~factor(combo)*factor(temp), random=~1|strain, data=data_exp, method="REML") 
anova(full)

#Addition of effect size stats
effectsize::eta_squared(full) 

#Remove outlier
data_exp_remO <-data_exp[(data_exp$Outlier=="n"),]# Removing an LN/LG outlier from the main dataset, but it is added back in below 

#check removing outlier not effecting results
full_remO <-nlme::lme(Growth_transformed~combo*factor(temp), random=~1|strain, data=data_exp_remO, method="REML")
anova(full_remO)

#Addition of effect size stats
effectsize::eta_squared(full_remO) 

#plotting
combo_names <- list('HNHG'="HL/HG",'HNLG'="HL/LG",'LNLG'="LL/LG")
combo_labeller <- function(variable,value){
  return(combo_names[value])
}

Fig_2 <- ggplot(data_exp_remO,aes(factor(temp),Growth,colour=combo))+
  geom_boxplot(outlier.shape = NA, position=position_dodge(width=0.75),color="black", size = 1) +
  theme_bw() +
  #ADD LINES
  geom_line(aes(group=strain, colour = strain), size=1.25, 
            alpha = 0.75,  show.legend = FALSE)+
  scale_colour_manual(name = "", values = c('FA19-02' = "#000000", 'FA19-06' ="#191919", 'FA19-05' ="#2b2b2b", 'SM19-01' ="#3e3e3e", 
                                            'F19-02' ="#525252", 'WH19-02' ="#676767", 'WH19-04' ="#7d7d7d",
                                            'SM19-03' ="#000000", 'WH19-01' ="#2b2b2b", 'CH19-02' ="#3e3e3e", 'AR19-01' = "#676767",
                                            'SM19-04' ="#7d7d7d",
                                            'CR19-06' ="#000000", 'G19-01' ="#191919", 'BU19-01' ="#2b2b2b", 'G19-02' ="#3e3e3e", 
                                            'CR19-05' ="#525252", 'BU19-04' ="#676767", 'BU19-02' ="#7d7d7d")) +
  #RESET COLOUR SCALE AND ADD POINTS
  new_scale_colour() +
  scale_color_manual(name="Bacterial Type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),
                     values=c("HNHG"="green","HNLG"="#56B4E9","LNLG"="blue"))+
  geom_point(aes(colour = combo), size = 3, shape = 21, stroke = 1.75) +
  #ADD OUTLIER
  geom_point(data = data.frame(temp = 24, Growth = 0.8184, combo = "LNLG"), colour="blue",size=3, stroke=1, shape=21) + #adding the outlier back in so it is on the graph but not connected by the lines 
  #THEME STUFF
  theme(strip.text.x = element_text(size = 16),plot.title=element_text(hjust=0.5),axis.text=element_text(size=16),axis.title=element_text(size=16))+
  xlab("Temperature of Growth Experiment (°C)")+
  ylab("Rate of Exponential Growth \n(per day)")+
  guides(size="none",alpha="none") +
  #FACET
  facet_grid(~combo, labeller=as_labeller(c('HNHG'="HL/HG",'HNLG'="HL/LG",'LNLG'="LL/LG")))

Fig_2

##Calculate how far out of range outlier is

data_exp_remO %>% 
  dplyr::filter(temp == "24") %>% 
  summarise(.by = combo, sd=sd(Growth), mean = mean(Growth))

data_exp_remO %>% 
  dplyr::filter(temp == "24") %>% 
  summarise(sd=sd(Growth), mean = mean(Growth))

#AT 24
#HNLG: sd = 0.09154848, mean = 0.2677000
#LNLG: sd = 0.04190859, mean = 0.1667833
#HNHG: sd = 0.09777963, mean = 0.2751714
#All: sd = 0.09208864, mean = 0.2369667

#outlier = 0.8184

(0.8184 - 0.2369667)/0.09208864 #sd from all strains mean
(0.8184 - 0.1667833)/0.04190859 #sd from mean of LL/LG

# Exponential phase GAM model selection ----------------------------------------------
#following peerj.com/articles/6876
#see also https://dfzljdn9uc3pi.cloudfront.net/2019/6876/1/fig-4-2x.jpg for explanation of models

#make sure data has variables in right formats
data_exp$strain <- as.factor(data_exp$strain)
data_exp$combo <- as.factor(data_exp$combo)

#model with only a global smoother
G <- gam(Growth ~ s(temp, k=4, m=2) + 
           s(strain, k = 19, bs = 're'), #random effects should always have k = to the number of levels of the factor
         method = 'REML',
         data = data_exp)

gratia::draw(G)
summary(G)

#global smoother plus group-level smoothers that have the same wiggliness - Global smoother with individual effects that have a Shared penalty
GS <-  gam(Growth ~ s(temp, k=4, m=2) +
             s(temp, combo, k=4, bs="fs", m=2) +
             s(strain, k = 19, bs = 're'), #random effects should always have k = to the number of levels of the factor
           method = 'REML',
           data = data_exp)

gratia::draw(GS)
summary(GS)
AIC(G, GS)

#A global smoother plus group-level smoothers with differing wiggliness - Global smoother with individual effects that have Individual penalties
GI <- gam(Growth ~ s(temp, k=4, m=2) +
            s(temp, by = combo, k=4, m=1) +
            s(strain, k = 19, bs = 're'), #random effects should always have k = to the number of levels of the factor
          method = 'REML',
          data = data_exp)

gratia::draw(GI)
summary(GI)
AIC(G, GI)

#Group-specific smoothers without a global smoother, but with all smoothers having the same wiggliness
S <- gam(Growth ~ s(temp, combo, k=4, bs="fs", m=2) +
           s(strain, k = 19, bs = 're'), #random effects should always have k = to the number of levels of the factor
         method = 'REML',
         data = data_exp)

gratia::draw(S)
summary(S)
AIC(G, S)
AIC(GS, S)

#Group-specific smoothers with different wiggliness
I <- gam(Growth ~ s(temp, by = combo, k=4, m=1) +
           s(strain, k = 19, bs = 're'), #random effects should always have k = to the number of levels of the factor
         method = 'REML',
         data = data_exp)

gratia::draw(I)
summary(I)
AIC(G, I)

#compare all models
AIC_table <- AIC(G,GS,GI,S,I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0))) %>% 
  arrange(-AIC)

AIC_table

#######################################################################################
# Figure 3 ----------------------------------------------------------------
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# This removes CR19-01, which is the sole LN/HG

#Check data for normality
hist(data$Growth)
shapiro.test(data$Growth)
normalTest(data$Growth, method = "da")

#boxcox transformation of data for better normality
b <- boxcox(lm(data$Growth +1 ~ 1))
lambda <- b$x[which.max(b$y)]
lambda

#check effect of boxcox transformation on data 
hist(((data$Growth + 1)^lambda - 1)/lambda)
ggdensity(((data$Growth + 1)^lambda - 1)/lambda, y = "density", fill = "grey") 
shapiro.test(((data$Growth + 1)^lambda - 1)/lambda)
normalTest(((data$Growth + 1)^lambda - 1)/lambda, method = "da")

#normality of data is improved (though not much), so use this transformation
#Transform data
data <- data %>%
  mutate(Growth_transformed = ((Growth + 1)^lambda - 1)/lambda) 
#while transformation doesn't make a huge difference the histogram and density looks a bit better and it is consistent to use this as used it earlier

#make Phase an ordered factor for when it is an effect
data$Phase <- factor(data$Phase, levels = c("Exp", "wk2", "wk3", "wk4"))
#make sure strain & type are factors too (otherwise GAM will error)
data$strain <- factor(data$strain)
data$combo <- factor(data$combo)

#LMEs ------------
#Growth rate with phase (week) as a random effect
full<-lme(Growth_transformed~factor(combo)*factor(temp), random=~1|Phase/strain, 
          data=data, method="REML")
#full<-lme(Growth_transformed~factor(combo)*factor(temp), random=~1|Phase, 
#          data=data, method="REML")
anova(full)
effectsize::eta_squared(full) 

#Model with phase
no_ints <- lme(Growth_transformed~ factor(combo) + factor(temp) + factor(Phase), 
              random=~1|strain,data=data,method="ML")
lim_ints <- lme(Growth_transformed~ factor(combo) + factor(temp) + factor(Phase)+ 
                  factor(combo):factor(temp) + factor(combo):factor(Phase), 
                random=~1|strain,data=data,method="ML")
anova(no_ints,lim_ints)
two_way<-lme(Growth_transformed~ factor(combo) + factor(temp) + factor(Phase) + 
            factor(combo):factor(temp) + factor(combo):factor(Phase) + factor(Phase):factor(temp), 
          random=~1|strain,data=data,method="ML")
anova(no_ints,two_way)
anova(lim_ints,two_way)
#full<-lme(Growth_transformed~ factor(combo)*factor(temp)*factor(Phase), 
#          random=~1|strain,data=data,method="ML") 
#anova(two_way, full)

#limited interaction model seems the best
#switch back to REML method for final version now you're not comparing models 
lim_ints <- lme(Growth_transformed~ factor(combo) + factor(temp) + factor(Phase)+ 
                  factor(combo):factor(temp) + factor(combo):factor(Phase), 
                random=~1|strain,data=data,method="REML")

anova(lim_ints) 
effectsize::eta_squared(lim_ints) 

#GAMs -----------------
#following peerj.com/articles/6876
#model with only a global smoother
G <- gam(Growth ~ s(temp, k=4) + s(strain, k = 19, bs = 're') + 
           s(Phase, k = 4, bs = 're'), 
         method = 'REML',
         data = data)

gam.check(G)
gratia::draw(G) 
summary(G)


#global smoother plus group-level smoothers that have the same wiggliness - Global smoother with individual effects that have a Shared penalty
GS <-  gam(Growth ~ s(temp, k=4, m=2) + s(strain, k = 19, bs = 're') +
             s(temp, combo, k=4, bs="fs", m=2) +
             s(Phase, k = 4, bs = 're'), 
           method = 'REML',
           data = data)
#ignore warning, this is because you have global and group level trends across temp for growth rate

gam.check(GS)
gratia::draw(GS)
summary(GS)
AIC(G, GS)

#A global smoother plus group-level smoothers with differing wiggliness - Global smoother with individual effects that have Individual penalties
GI <- gam(Growth ~ s(temp, k=4, m=2) + s(strain, k = 19, bs = 're') +
            s(temp, by = combo, k=4, bs="tp", m=1) +
            s(Phase, k = 4, bs = 're'),
          method = 'REML',
          data = data)

gam.check(GI)
gratia::draw(GI)
summary(GI)
AIC(G, GI)

#Group-specific smoothers without a global smoother, but with all smoothers having the same wiggliness
S <- gam(Growth ~ s(temp, combo, k=4, bs="fs", m=2) + s(strain, k = 19, bs = 're') +
           s(Phase, k = 4, bs = 're'), 
         method = 'REML',
         data = data)

gam.check(S)
gratia::draw(S)
summary(S)
AIC(G, S)
AIC(GS, S)

#Group-specific smoothers with different wiggliness
I <- gam(Growth ~ s(temp, by = combo, k=4, m=1) + s(strain, k = 19, bs = 're') +
           s(Phase, k = 4, bs = 're'), 
         method = 'REML',
         data = data)

gam.check(I)
gratia::draw(I)
summary(I)
AIC(G,I)
AIC(GI,I)

#compare all models
AIC_table <- AIC(G,GS,GI,S,I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0))) %>% 
  arrange(-AIC)

AIC_table

#NOW FOR THE WITH PHASE MODEL
#model with only a global smoother
G <- gam(Growth ~ s(temp, k=4) + 
           s(temp, Phase, k = 4, bs = "fs", m = 2) +
           s(strain, k = 19, bs = 're'), #random effects should always have k = to the number of levels of the factor
         method = 'REML',
         data = data)

gratia::draw(G)
summary(G)

#global smoother plus group-level smoothers that have the same wiggliness - Global smoother with individual effects that have a Shared penalty
GS <-  gam(Growth ~ s(temp, k=4, m=2)  + #Phase + 
             s(temp, combo, k=4, bs="fs", m=2) +
             s(temp, Phase, k=4, bs="fs", m=2) +
             s(strain, k = 19, bs = 're'), #random effects should always have k = to the number of levels of the factor
           method = 'REML',
           data = data)

gratia::draw(GS)
summary(GS)
AIC(G, GS)

#A global smoother plus group-level smoothers with differing wiggliness - Global smoother with individual effects that have Individual penalties
GI <- gam(Growth ~ s(temp, k=4, m=2) + s(strain, k = 19, bs = 're') +
            s(temp, by = combo, k=4, m=1) + 
            s(temp, Phase, k=4, bs="fs", m=1),
          method = 'REML',
          data = data)

gam.check(GI)
gratia::draw(GI)
summary(GI)
AIC(G, GI)

#Group-specific smoothers without a global smoother, but with all smoothers having the same wiggliness
S <- gam(Growth ~ s(temp, combo, k=4, bs="fs", m=2) + s(strain, k = 19, bs = 're') +
           s(temp, Phase, k=4, bs="fs", m=1), 
         method = 'REML',
         data = data)

gam.check(S)
gratia::draw(S)
summary(S)
AIC(G, S)
AIC(GS, S)

#Group-specific smoothers with different wiggliness
I <- gam(Growth ~ s(temp, by = combo, k=4, m=1) + s(strain, k = 19, bs = 're') +
           s(temp, Phase, k=4, bs="fs", m=1), 
         method = 'REML',
         data = data)

gam.check(I)
gratia::draw(I)
summary(I)
AIC(G,I)
AIC(GI,I)

#compare all models
AIC_table <- AIC(G,GS,GI,S,I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0))) %>% 
  arrange(-AIC)

AIC_table

#POST HOCs for LMEs ----------------
full<-lme(Growth_transformed~factor(combo)*as.numeric(temp), random=~1|Phase/strain, data=data, method="REML")
emmeans(full, pairwise ~ combo, adjust = "tukey")
regrid(emmeans(full, pairwise ~ combo, adjust = "tukey"))
#to get comparisons across the temperatures use temp as an ordered factor
data$tempF <- factor(data$temp, levels = c("20", "24", "28", "32"))
full<-lme(Growth_transformed~factor(combo)*tempF, random=~1|Phase/strain, data=data, method="REML")
#double check this doesn't have a meaningful change on the model
anova(full) #same terms remain significant
emmeans(full, pairwise ~ combo, adjust = "tukey") #same posthoc comparisons remain significant
regrid(emmeans(full, pairwise ~ combo, adjust = "tukey"))
#so good to go ahead with comparing within temperatures
emmeans(full, pairwise ~ combo|tempF, adjust = "tukey")
regrid(emmeans(full, pairwise ~ combo|tempF, adjust = "tukey"))

#full<-lme(sqrt(Growth+1)~ combo + factor(temp) + factor(Phase) + combo:factor(temp) + combo:factor(Phase), 
#          random=~1|strain,data=data,method="REML")
#emmeans(full, pairwise ~ combo, adjust = "tukey")
#emmeans(full, pairwise ~ Phase, adjust = "tukey")
#emmeans(full, pairwise ~ combo|temp, adjust = "tukey")
#emmeans(full, pairwise ~ combo|Phase, adjust = "tukey")


#Cumulative Link Mixed Model -------------
#these are useful for ordered categorical data which is what the crash data is!
library("ordinal")
#Cumulative Link Mixed Model for ranked categorical/ordinal data
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# This removes CR19-01, which is the sole LN/HG

data$NegativePositiveGrowth<-as.factor(data$NegativePositiveGrowth)
null<-clmm(NegativePositiveGrowth~(1|Phase)+(1|strain),data=data)
mid_1<-clmm(NegativePositiveGrowth~as.factor(combo)+(1|Phase)+(1|strain),data=data)
anova(null,mid_1)#p = 0.007084 (mid_1 wins) (i.e. combination is a significant predictor of crash)

mid_2<-clmm(NegativePositiveGrowth~as.factor(combo)+factor(temp)+(1|Phase)+(1|strain),data=data)
anova(mid_1,mid_2)#p = 0.007804  (mid_2 wins) (i.e. temperature is a significant predictor of crash)

full<-clmm(NegativePositiveGrowth~as.factor(combo)*factor(temp)+(1|Phase)+(1|strain),data=data)
anova(mid_2,full) #p = 0.4501

#significance of fixed effects
library(car)
library(RVAideMemoire)

Anova.clmm(full,
           type = "II")

#Plotting ==================
combo_names <- list('HNHG'="HL/HG",'HNLG'="HL/LG",'LNLG'="LL/LG")
combo_labeller <- function(variable,value){
  return(combo_names[value])
}
week_names <- list('Exp'="Exponential (Week 1)",'wk2'="Week 2",'wk3'="Week 3",'wk4'="Week 4")
week_labeller <- function(variable,value){
  return(week_names[value])
}

Fig_3<- ggplot(data,aes(factor(temp),Growth))+
  geom_boxplot(aes(fill=combo), outlier.shape = NA,
               position=position_dodge(width=0.75), 
               cex = 0.75, alpha = 0.25) +
  scale_fill_manual(name="Bacterial Type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),values=c("HNHG"="green","HNLG"="#56B4E9","LNHG"="black","LNLG"="blue")) +
  geom_hline(yintercept=0, linetype="dashed",color = "red", size=1)+
  #RESET COLOUR SCALE AND ADD POINTS
  new_scale_fill() +
  #geom_jitter(aes(fill=NegativePositiveGrowth, colour=combo),
  #            size=2.5, stroke=2.5, 
  #            shape=21,width=0.36,height=0.0)+
  geom_point(aes(fill=NegativePositiveGrowth, colour=combo),
             size=2.5, stroke=2.5, 
             shape=21, position = position_dodge2(width = 0.6)) +
  theme_bw()+
  theme(strip.text.y=element_text(size=10),strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),axis.text=element_text(size=16),
        axis.title=element_text(size=16), legend.position = "bottom")+
  labs(x=expression(paste("Temperature (", degree~C, ")")),y="Growth Rate") +
  scale_fill_manual(name="Growth Rate (per day)",values=c("negative"="black","positive"="white"))+
  scale_color_manual(name="Bacterial Type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),values=c("HNHG"="green","HNLG"="#56B4E9","LNHG"="black","LNLG"="blue"))+
  guides(size="none",alpha="none")+
  facet_grid(Phase~combo,labeller = labeller(Phase = week_labeller, combo = combo_labeller))

Fig_3


# Separate statistical models for each GE-type ####
#### Within each category of Microcystis, does exponential growth differ by temperature treatment?
data <- read.csv("Sara_Revised_New_FixedSM.csv")

#boxcox transformation of data for better normality
b <- boxcox(lm(data$Growth +1 ~ 1))
lambda <- b$x[which.max(b$y)]
lambda

#Transform data
data <- data %>%
  mutate(Growth_transformed = ((Growth + 1)^lambda - 1)/lambda) 
#while transformation doesn't make a huge difference the histogram and density looks a bit better and it is consistent to use this as used it earlier

#separate into types
HNHG<-data[(data$combo=="HNHG"),]  
HNLG<-data[(data$combo=="HNLG"),]  
LNHG<-data[(data$combo=="LNHG"),] # no stats because just one strain 
LNLG<-data[(data$combo=="LNLG"),]  

#Linear mixed-effects models factoring in strain identity
full<-lme(Growth_transformed~factor(temp),random=~1|Phase/strain,data=HNHG,method="REML")
anova(full)
effectsize::eta_squared(full) 
emmeans(full, list(pairwise ~ temp), adjust = "tukey")

full<-lme(Growth_transformed~factor(temp),random=~1|Phase/strain,data=HNLG,method="ML")
anova(full)
#effectsize::eta_squared(full) 
#emmeans(full, list(pairwise ~ temp), adjust = "tukey")

full<-lme(Growth_transformed~factor(temp),random=~1|Phase/strain,data=LNLG,method="ML")
anova(full)
#effectsize::eta_squared(full) 
#emmeans(full, list(pairwise ~ temp), adjust = "tukey")


#######################################################################################
# Figure 4 ----------------------------------------------------------------
##### Load data #####
dat <- read.csv("cq.data.csv")
met <- read.csv("Microcystis growth and metadata.csv")

#### Process data ####
#process metadata
met <- met %>%
  dplyr::select(!c("lake", "strain", "combo_temp", "quality", "combo", "genotype",
                   "Trophic", "temp"))

sdat <- dat %>%
  #remove blanks
  dplyr::filter(!Sample %in% c("", "Blank", "blank", "BLANK")) %>%
  #the group by factors identifying unique treatments, also include any other factors that need to be kept in the output
  group_by(flask, Target, lake, strain, temp, genotype, Trophic, combo) %>%
  #summarise gene expression within each grouping
  summarise(mean.b20.CT.metric = mean(Final.metric...20.as.control, na.rm = TRUE)) %>%
  #join metadata and gene expression data
  left_join(y = met, by = join_by(flask == flask)) %>%
  #add a column that gives the type of gene
  mutate(Type = case_when(Target == "Hep" ~ "Toxin Gene", .default = "HSP Gene"))  %>%
  #remove unecessary columns
  dplyr::select(!c("X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6", "X.7",))  %>%
  dplyr::filter(temp != "20") %>% 
  #remove the reference condition and gene which were used for calculating the gene expression metric
  dplyr::filter(!Target %in% c("", "rpoA"))

#create labels for the temp facets
temp.labels <- c("24\u00B0C", "28\u00B0C", "32\u00B0C")
names(temp.labels) <- c("24", "28", "32")

#create a vector of gene targets to help with plotting
primer <- unique(sdat$Target)

#make data frames to tally highest genotype across facets
sums <- sdat %>%
  group_by(temp, Target, combo) %>%
  summarise(mean = mean(mean.b20.CT.metric))

#correct genotype label 
sdat <- sdat %>%
  mutate(combo = case_when(
    combo == "HNHG" ~ "HL/HG",
    combo == "HNLG" ~ "HL/LG",
    combo == "LNLG" ~ "LL/LG"
  ))

#### Normalization of data ####
###Box cox
#Check data as is
hist(sdat$mean.b20.CT.metric)
shapiro.test(sdat$mean.b20.CT.metric)
normalTest(sdat$mean.b20.CT.metric, method = "da")

#boxcox transformation
b <- boxcox(lm(sdat$mean.b20.CT.metric ~ 1))
lambda <- b$x[which.max(b$y)]
lambda

#check effect of boxcox transformation on data 
hist((sdat$mean.b20.CT.metric^lambda - 1)/lambda)
ggdensity((sdat$mean.b20.CT.metric^lambda - 1)/lambda, y = "density", fill = "grey") 
shapiro.test((sdat$mean.b20.CT.metric^lambda - 1)/lambda)
normalTest((sdat$mean.b20.CT.metric^lambda - 1)/lambda, method = "da")

#Transform data
tdat <- sdat %>%
  mutate(CT.metric = (mean.b20.CT.metric^lambda - 1)/lambda) 

#make sure various variables are in the correct format
tdat$temp <- as.factor(tdat$temp)
tdat$combo <- as.factor(tdat$combo)
tdat$strain <- as.factor(tdat$strain)
tdat$mean.b20.CT.metric <- as.numeric(tdat$mean.b20.CT.metric)

#remove Hep from the general model
hsps <- tdat %>%
  dplyr::filter(!Target %in% c("Hep"))

ggdensity(hsps$CT.metric, fill = "grey")

#### Model for heat shock protein genes together ####
m1 <- lme(CT.metric ~ factor(temp)*factor(combo), random = list(Target = ~1, strain = ~1), method = "ML", data = hsps)

#Test assumptions
res_aov <- residuals(m1)

# histogram
hist(res_aov)
ggdensity(res_aov, fill = "grey")

# QQ-plot
car::qqPlot(res_aov,
            id = FALSE, # id = FALSE to remove point identification
)

#perform statistics
anova(m1)
effectsize::eta_squared(m1) 

#### Model individually for each heat shock protein gene ####

#make a vector of the target genes to loop through
primer <- unique(tdat$Target)

#use a for loop to simplify creation and testing of the same model for each gene
for(p in 1:length(primer)){
  #create a temporary dataframe holding the information for only a single gene target
  tmp <- tdat %>% 
    dplyr::filter(Target == primer[p])
  
  #make a model for this gene target
  mp <- lme(CT.metric ~ factor(temp)*factor(combo), random = list(strain = ~1), method = "ML", data = tmp)
  
  #Test assumptions within the subset
  res_aov <- residuals(mp)
  
  par(mfrow = c(1, 2)) # combine plots
  
  # histogram
  hist(res_aov, xlab = "ANOVA Residuals", main = primer[p], breaks = 10)
  
  # QQ-plot
  car::qqPlot(res_aov,
              id = FALSE, # id = FALSE to remove point identification
              main = primer[p],
              xlab = "Norm Quantiles",
              ylab = "ANOVA residuals"
  )
  #pring results of staistical testing
  print("") 
  print("")
  print("")
  print(paste("For target of primer", primer[p], "statistical results are:"))
  print(anova(mp))
  print(effectsize::eta_squared(mp))
  print("")
  print("")
  print("")
}

#### Post hoc tests ####
#Post hoc tests for all HSP target genes together
emmeans(m1, pairwise ~ temp, adjust = "tukey")
emmeans(m1, pairwise ~ temp|combo, adjust = "tukey")
emmeans(m1, pairwise ~ combo|temp, adjust = "tukey")

#Post hoc tests for each HSP target gene separately
#make a vector of the target genes to loop through
primer <- unique(tdat$Target)

#loop through each primer
for(p in 1:length(primer)){
  #create a temporary dataframe holding the information for only a single gene target
  tmp <- tdat %>% 
    dplyr::filter(Target == primer[p])
  
  #make a model for this gene target
  mp <- lme(CT.metric ~ factor(temp)*factor(combo), random = list(strain = ~1), method = "ML", data = tmp)
  
  #Test assumptions within smaller dataset
  res_aov <- residuals(mp)
  
  par(mfrow = c(1, 2)) # combine plots
  
  # histogram
  hist(res_aov)
  
  # QQ-plot
  car::qqPlot(res_aov,
              id = FALSE, # id = FALSE to remove point identification
              main = primer[p]
  )
  
  print("")
  print("")
  print(paste("For target of primer", primer[p], "statistical results are:"))
  print(emmeans(mp, pairwise ~ temp|combo, adjust = "tukey"))
  #print(emmeans(mp, pairwise ~ temp, adjust = "tukey"))
  #print(emmeans(mp, pairwise ~ combo|temp, adjust = "tukey"))
  print("")
  print("")
  
}
# Figure 4 plotting  ####
#### Primer 1 = ClpB ####
#set target
p <- 1

#make temporary dataframes of the data for set target into subsets
tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

#plot subsets separately to allow axes to vary
p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.889))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.889))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28  

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.9287))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32  

#assign a name to the three plots joined together
nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)

Plot_ClpB


#### Primer 2 = DnaJ ####
p <- 2

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.93))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.943))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.929))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 3 = DnaK-fp ####
p <- 3

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.947))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.8999))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.9289))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 4 = DnaK1 ####
p <- 4

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.945))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.947))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.93))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 5 = DnaK3 ####
p <- 5

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.95))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.89))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.93))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 6 = GroEL ####
p <- 6

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.94469))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.94))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0.1, 0.855))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 7 = GroES ####
p <- 7

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.8925))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.999))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.925))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 8 = GrpE ####
p <- 8

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.95))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.8945))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.938))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 9 = Hep ####
p <- 9

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.95))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0, 0.944))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.825))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer  10 = HrcA ####
p <- 10

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.898))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0.1, 0.899))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0.1, 0.86))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer 11 = Hsp20 ####
p <- 11

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.9445))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0.1, 0.89))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.93))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer  12 = HspA ####
p <- 12

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0, 0.838))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0.1, 0.915))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0.1, 0.941))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Primer  13 = HtpG ####
p <- 13

tmp24 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp24$mean.b20.CT.metric, c(0.1, 0.8896))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp28$mean.b20.CT.metric, c(0.1, 0.944))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0, 0.8))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


#### Put the smaller plots for all the HSPs together ####
Plot_ClpB | Plot_DnaJ | `Plot_DnaK-fp` | Plot_DnaK1 | Plot_DnaK3 | Plot_GroEL | Plot_GroES | Plot_GrpE | Plot_HrcA | Plot_Hsp20 | Plot_HspA | Plot_HtpG

#plot_spacer() attempt
plot_spacer() + Plot_ClpB + plot_spacer() + Plot_DnaJ + plot_spacer() + `Plot_DnaK-fp` + plot_spacer() + Plot_DnaK1 + plot_spacer() + Plot_DnaK3 + plot_spacer() + Plot_GroEL + plot_spacer() + Plot_GroES + plot_spacer() + Plot_GrpE + plot_spacer() + Plot_HrcA + plot_spacer() + Plot_Hsp20 + plot_spacer() + Plot_HspA + plot_spacer() + Plot_HtpG + plot_spacer() + plot_layout(widths = c(-2.05, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -6)) + plot_layout(nrow = 1, ncol = 27)

#add back legends and facet headings in powerpoint for a smoother final effect



#######################################################################################
# Figure 5 ------------------------------------------------------------------
#expression data
tdat <- read.csv("tdat.csv")

#growth data
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# This removes CR19-01, which is the sole LN/HG

#boxcox transformation of data for better normality (precident and testing done previously)
b <- boxcox(lm(data$Growth +1 ~ 1))
lambda <- b$x[which.max(b$y)]
lambda

#Transform data
data <- data %>%
  #while transformation doesn't make a huge difference the histogram and density looks a bit better and it is consistent to use this as used it earlier
  mutate(Growth_transformed = ((Growth + 1)^lambda - 1)/lambda) %>% 
  #also remove unecessary columns of the data
  dplyr::select(flask, strain, Lake, combo_temp, combo, genotype, Trophic, temp, Crash_Numeric, Crash, Growth_transformed, Phase) %>% 
  #and widen data so each Phase has its own growth data column
  dcast(flask + strain + temp + genotype + Trophic + combo + Crash_Numeric ~ Phase, 
        value.var = "Growth_transformed") %>% 
  #rename these to make it clear they are the transformed values
  rename(Exp_transformed = Exp, wk2_transformed = wk2, wk3_transformed = wk3, wk4_transformed = wk4)
  

#'widen' expression data
tdat_wide <- dcast(data = tdat, flask + combo + temp + genotype + Trophic + lake + strain + Crash_Numeric + ExpPhase + wk2 + wk3 + wk4 ~ Target, 
                   value.var = "CT.metric") %>% 
  left_join(y = data, by = join_by(flask, strain, temp, genotype, Trophic, combo, Crash_Numeric))


full_MANOVA <- manova(cbind(ClpB, DnaJ, `DnaK-fp`, DnaK1, DnaK3, GroEL, GroES, GrpE, Hep, HrcA, Hsp20, HspA, HtpG, 
                            Exp_transformed, wk2_transformed, wk3_transformed, wk4_transformed) 
                      ~ factor(combo)*factor(temp), 
                      data = tdat_wide)

summary(full_MANOVA)
effectsize::eta_squared(full_MANOVA)

expression_MANOVA <- manova(cbind(ClpB, DnaJ, `DnaK-fp`, DnaK1, DnaK3, GroEL, 
                                  GroES, GrpE, Hep, HrcA, Hsp20, HspA, HtpG) 
                            ~ factor(combo)*factor(temp), 
                            data = tdat_wide)

summary(expression_MANOVA)
effectsize::eta_squared(expression_MANOVA)

growth_MANOVA <- manova(cbind(Exp_transformed, wk2_transformed, wk3_transformed, wk4_transformed) 
                        ~ factor(combo)*factor(temp), 
                        data = tdat_wide)

summary(growth_MANOVA)
effectsize::eta_squared(growth_MANOVA)

#post hoc for full model
#test for homogeneity of within-group covariance
t_mat <- as.matrix(tdat_wide[, c("ClpB", "DnaJ", "DnaK-fp", "DnaK1", "DnaK3", "GroEL", 
                                 "GroES", "GrpE", "Hep", "HrcA", "Hsp20", "HspA", "HtpG", "Exp_transformed", 
                                 "wk2_transformed", "wk3_transformed", "wk4_transformed")])

t_mat_dist <- dist(t_mat)

vegan::permutest(vegan::betadisper(t_mat_dist, tdat_wide$combo)) 
#dispersion is okay so can go ahead with lda!


full_post_hoc <- lda(combo~cbind(ClpB, DnaJ, `DnaK-fp`, DnaK1, DnaK3, GroEL, GroES, GrpE, Hep, HrcA, Hsp20, HspA, HtpG, 
                                 Exp_transformed, wk2_transformed, wk3_transformed, wk4_transformed), 
                     data = tdat_wide,
                     CV = F)
full_post_hoc


full_plot_lda <- data.frame(tdat_wide[, "combo"], lda = predict(full_post_hoc)$x)
colnames(full_plot_lda) <- c("combo", "lda.LD1", "lda.LD2" )

#plot
ggplot(full_plot_lda) + 
  geom_point(aes(x = lda.LD1, y = lda.LD2, colour = combo), size = 4) +
  scale_color_manual(name="GE-type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),
                     values=c("HNHG"="green","HNLG"="#56B4E9","LNLG"="blue")) +
  labs(colour = "GE-type",
       x = "Linear Discriminate Analysis Axis 1 \n(Proportion of trace = 0.75)",
       y = "Linear Discriminate Analysis Axis 2 \n(Proportion of trace = 0.25)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16))

#biplot full
full_arrows <- data.frame("x0" = 0, 
                          "y0" = 0,
                          "x1" = coef(full_post_hoc)[,1],
                          "y1" = coef(full_post_hoc)[,2],
                          "labs" = c("ClpB", "DnaJ", "DnaK-fp", "DnaK1", "DnaK3", "GroEL", 
                                     "GroES", "GrpE", "Hep", "HrcA", "Hsp20", "HspA", "HtpG", "Week 1", 
                                     "Week 2", "Week 3", "Week 4"))

ggplot(full_plot_lda) + 
  geom_segment(data = full_arrows, aes(x = x0, xend = x1, y = y0, yend = y1),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, colour = "black") +
  geom_point(aes(x = lda.LD1, y = lda.LD2, fill = combo), size = 5, pch = 21) +
  geom_text_repel(data = full_arrows, aes(x = x1 + (x1/10), y = y1, label = labs), cex = 4, colour = "darkred",
                  fontface = "bold", max.overlaps = 21, min.segment.length = 0.3, #nudge_x = -0.15,
                  segment.curvature = -0.1, point.size = 5) +
  scale_fill_manual(name="GE-type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),
                     values=c("HNHG"="green","HNLG"="#56B4E9","LNLG"="blue")) +
  labs(fill = "GE-type",
       x = "Linear Discriminate Analysis Axis 1 (74.7%)",
       y = "Linear Discriminate Analysis Axis 2 (25.3%)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16))

#unlabeled version for powerpoint improved labels
ggplot(full_plot_lda) + 
  geom_segment(data = full_arrows, aes(x = x0, xend = x1, y = y0, yend = y1),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, colour = "black") +
  geom_point(aes(x = lda.LD1, y = lda.LD2, fill = combo), size = 5, pch = 21) +
  #geom_text_repel(data = full_arrows, aes(x = x1 + (x1/10), y = y1, label = labs), cex = 4, colour = "darkred",
  #                fontface = "bold", max.overlaps = 21, min.segment.length = 0.3, #nudge_x = -0.15,
  #                segment.curvature = -0.1, point.size = 5) +
  scale_fill_manual(name="Bacterial Type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),
                    values=c("HNHG"="green","HNLG"="#56B4E9","LNLG"="blue")) +
  labs(fill = "Bacterial Type",
       x = "Linear Discriminate Analysis Axis 1 (74.7%)",
       y = "Linear Discriminate Analysis Axis 2 (25.3%)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.position = "bottom")

#test how good classifications are
test <- lda(combo~cbind(ClpB, DnaJ, `DnaK-fp`, DnaK1, DnaK3, GroEL, 
                        GroES, GrpE, Hep, HrcA, Hsp20, HspA, HtpG, Exp_transformed, 
                        wk2_transformed, wk3_transformed, wk4_transformed), 
            data = tdat_wide,
            CV = TRUE)

summary(test)

test$class

test.table <- table(tdat_wide$combo, test$class)
#percent correct per category
diag(prop.table(test.table, 1))
#this was far better when using the untransformed growth data...
# total percent correct
sum(diag(prop.table(test.table)))

#######################################################################################
# Figure S1 ---------------------------------------------------------------
full_tree<-read.tree("RAxML_bipartitions.Testing_Aug_23_2024_FigS1A_my_samples_only_10kboots")
#full_tree<-read.tree("final whole tree MLST.txt")
metadata_full<-read.csv("Microcystis_mapping_full_tree.csv")
dd_full<-data.frame(Strain=metadata_full$Strain,TreeGroups=metadata_full$Tree_Groups)
head(dd_full)
full_tree_p<-ggtree(full_tree) + geom_tree() + theme_tree() + geom_tiplab(size=5,hjust = -0.25)
Fig_S1A <-full_tree_p %<+% dd_full + 
  geom_tippoint(aes(color=TreeGroups,alpha=0.5),size=4)+ 
  geom_text2(aes(label=label,x=branch,hjust=1.3, vjust= -0.95,subset = !is.na(as.numeric(label)) & as.numeric(label) >50),size=5)+
  #theme_grey()+
  theme(legend.position="right")+
  scale_color_manual(values = c("nutrient-rich"="green", "oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","none"="gray"))+
  guides(color=FALSE,alpha=FALSE)+
  geom_treescale(x = 0.07, y = 0)+
  xlim(NA,0.15)
Fig_S1A

### Fig S1B

full_its_tree<-read.tree("RAxML_bipartitions.Testing_Aug_17_2024_FigS1B_my_samples_only_10kboots")
ggtree(full_its_tree)+geom_tiplab()
metadata_full<-read.csv("Microcystis_mapping_for_ITSc.csv")
dd_full<-data.frame(Strain=metadata_full$Strain,TreeGroups=metadata_full$Tree_Groups)
full_tree_p<-ggtree(full_its_tree) + geom_tree() + theme_tree() + geom_tiplab(size=5,hjust = -0.015)
my_full_tree <- full_tree_p %<+% dd_full +
  geom_tippoint(aes(color=TreeGroups,alpha=0.5),size=4)+ 
  theme(legend.position="right")+
  geom_treescale()+
  geom_text2(aes(label=label,x=branch,hjust=1.3, vjust= -0.95,subset = !is.na(as.numeric(label)) & as.numeric(label) >50),size=5)+
  scale_color_manual(values = c("nutrient-rich"="green", "oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","none"="gray"))+
  guides(color=FALSE,alpha=FALSE)#+xlim(NA,0.15)
my_full_tree


#######################################################################################
# Figure S3 ---------------------------------------------------------------
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# This removes CR19-01, which is the sole LN/HG

combo_names <- list('HNHG'="HL/HG",'HNLG'="HL/LG",'LNLG'="LL/LG")
combo_labeller <- function(variable,value){
  return(combo_names[value])
}
week_names <- list('Exp'="Exponential (Week 1)",'wk2'="Week 2",'wk3'="Week 3",'wk4'="Week 4")
week_labeller <- function(variable,value){
  return(week_names[value])
}


figS3 <- data %>% 
  dplyr::filter(!strain == "CR19-01") %>% 
  ggplot(aes(x = Phase, y = Growth)) +
  geom_line(aes(group = strain, colour = strain), show.legend = FALSE,
            size = 1.5, alpha = 0.75) +
  scale_colour_manual(name = "", values = c('FA19-02' = "#000000", 'FA19-06' ="#191919", 'FA19-05' ="#2b2b2b", 'SM19-01' ="#3e3e3e", 
                                            'F19-02' ="#525252", 'WH19-02' ="#676767", 'WH19-04' ="#7d7d7d",
                                            'SM19-03' ="#000000", 'WH19-01' ="#2b2b2b", 'CH19-02' ="#3e3e3e", 'AR19-01' = "#676767",
                                            'SM19-04' ="#7d7d7d",
                                            'CR19-06' ="#000000", 'G19-01' ="#191919", 'BU19-01' ="#2b2b2b", 'G19-02' ="#3e3e3e", 
                                            'CR19-05' ="#525252", 'BU19-04' ="#676767", 'BU19-02' ="#7d7d7d")) +
  #RESET COLOUR SCALE AND ADD POINTS
  new_scale_colour() +
  scale_color_manual(name="Bacterial Type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),
                     values=c("HNHG"="green","HNLG"="#56B4E9","LNLG"="blue"))+
  geom_point(aes(colour = combo), pch = 21, stroke = 2.5, size = 3, show.legend = FALSE) +
  geom_hline(yintercept = 0, col = "red", lty = 2) +
  theme_bw() +
  labs(y = "Growth Rate", x = "") +
  theme(strip.text.y=element_text(size=10),strip.text.x = element_text(size = 16),plot.title=element_text(hjust=0.5),
        axis.text=element_text(size=16),axis.title=element_text(size=16), axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3", "Week 4")) +
  facet_grid(temp~combo, labeller = labeller(temp = c("20" = "20°C", "24" = "24°C", "28" = "28°C", "32" = "32°C"), 
                                             combo = combo_labeller))

figS3

tiff("./Final figure files/FigS3.tiff", units="px", width=3000, height=2500, res=300)
figS3
dev.off()

#######################################################################################
# Figure S4 ---------------------------------------------------------------
#Read in and process data
hsps<-read.csv("Heat_Stress_pfam_list.csv",head=F)
hsps_list<-hsps$V2
pfam_table<-read.csv("R_file10r_PFAMs_by_strain.csv",head=F)
pfam_to_label<-pfam_table[,5:6]

#shorten the name of "Zn-dependent protease with chaperone function"
pfam_to_label[12,1]
pfam_to_label[12,1] <- "ZnPs*"

pfam_table<-pfam_table[,1:2]
length(hsps_list)
head(pfam_table)
pfam_table<-subset(pfam_table, subset = V2 %in% c("PF00012.19","PF00011.20","PF01430.18","PF00183.17","PF02518.25","PF13589.5","PF00118.23","PF01025.18","PF03724.15","PF14355.5","PF01435.17","PF08240.11","PF00226.30","PF00684.18","PF01556.17","PF00227.25","PF02190.15","PF00595.23","PF06480.14","PF00140.19","PF04539.15","PF04542.13","PF04545.15","PF08281.11","PF02617.16","PF06689.12","PF02861.19","PF00574.22","PF10431.8"))
pfam_table$V2<-as.factor(pfam_table$V2)
length(levels(pfam_table$V2))
hsp_table<-as.data.frame(table(pfam_table))
head(hsp_table)
hsp_table$Strain<-hsp_table$V1
hsp_table$V1<-NULL
mapping<-read.csv("Microcystis_mapping_full_tree.csv",head=T)
head(mapping)
hsp_table<-dplyr::full_join(hsp_table,mapping,by="Strain")
hsp_table<-hsp_table[!is.na(hsp_table$V2),]
tail(hsp_table)
head(pfam_to_label)
colnames(pfam_to_label)[2] <- "V2"
hsp_table<-dplyr::full_join(hsp_table,pfam_to_label,by="V2")
hsp_table<-hsp_table[!is.na(hsp_table$Strain),]

#plot the results
hsp_table %>%
  dplyr::filter(!V5 == "htpX") %>%
  ggplot(aes(V5, Freq, colour=Tree_Groups)) + 
  geom_boxplot(outlier.shape = NA,position=position_dodge(width = 1), cex = 1)+
  scale_colour_manual(name= expression(paste("Bacterial type of ", italic("Microcystis aeruginosa"))),
                      values=c("nutrient-rich"="green","pseudo-oligotrophic"="#56B4E9","oligotrophic"="blue"),
                      labels=c("nutrient-rich"="HL/HG","pseudo-oligotrophic"="HL/LG", "oligotrophic"="LL/LG"))+
  theme_bw()+
  xlab("\nProtein Family involved in Heat Shock Regulation")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = , vjust = 0.5, hjust=1), axis.text.y = element_text(angle = ))+
  ylab("Frequency of Genes per MAG involved in Heat Shock Regulation\n") +
  coord_flip()


#######################################################################################
# Figure S5 ---------------------------------------------------------------
##### Load data #####
dat <- read.csv("cq.data.csv")
met <- read.csv("Microcystis growth and metadata.csv")

#### Process data ####
#process metadata
met <- met %>%
  dplyr::select(!c("lake", "strain", "combo_temp", "quality", "combo", "genotype",
                   "Trophic", "temp"))

sdat <- dat %>%
  #remove blanks
  dplyr::filter(!Sample %in% c("", "Blank", "blank", "BLANK")) %>%
  #the group by factors identifying unique treatments, also include any other factors that need to be kept in the output
  group_by(flask, Target, lake, strain, temp, genotype, Trophic, combo) %>%
  #summarise gene expression within each grouping
  summarise(mean.b20.CT.metric = mean(Final.metric...20.as.control, na.rm = TRUE)) %>%
  #join metadata and gene expression data
  left_join(y = met, by = join_by(flask == flask)) %>%
  #add a column that gives the type of gene
  mutate(Type = case_when(Target == "Hep" ~ "Toxin Gene", .default = "HSP Gene"))  %>%
  #remove unecessary columns
  dplyr::select(!c("X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6", "X.7",))  %>%
  dplyr::filter(temp != "20") %>% 
  #remove the reference condition and gene which were used for calculating the gene expression metric
  dplyr::filter(!Target %in% c("", "rpoA"))

#create labels for the temp facets
temp.labels <- c("24\u00B0C", "28\u00B0C", "32\u00B0C")
names(temp.labels) <- c("24", "28", "32")

#create a vector of gene targets to help with plotting
primer <- unique(sdat$Target)

#make data frames to tally highest genotype across facets
sums <- sdat %>%
  group_by(temp, Target, combo) %>%
  summarise(mean = mean(mean.b20.CT.metric))

#correct genotype label 
sdat <- sdat %>%
  mutate(combo = case_when(
    combo == "HNHG" ~ "HL/HG",
    combo == "HNLG" ~ "HL/LG",
    combo == "LNLG" ~ "LL/LG"
  ))

# Actual Plotting ####
###### Primer 1 = ClpB #### 
p <- 1

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.945))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.926))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG  

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.84))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG  

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)

Plot_ClpB


###### Primer 2 = DnaJ #### 
p <- 2

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 3)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 2)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 6)) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)

###### Primer 3 = DnaK-fp #### 
p <- 3

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 3)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 4)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 8)) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer 4 = DnaK1 #### 
p <- 4

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.959))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.999))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.945))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer 5 = DnaK3 #### 
p <- 5

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.92))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.93))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.9))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer 6 = GroEL #### 
p <- 6

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.955))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.9))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.83))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer 7 = GroES #### 
p <- 7

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.9))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.91))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.99))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer 8 = GrpE #### 
p <- 8

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.97))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.92))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.91))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer 9 = Hep #### 
p <- 9

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 1000)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 12)) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer  10 = HrcA #### 
p <- 10

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 25)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 4)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 8)) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer 11 = Hsp20 #### 
p <- 11

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.9))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 40)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 125)) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer  12 = HspA #### 
p <- 12

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 10)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.925))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.94))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)


###### Primer  13 = HtpG #### 
p <- 13

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 4)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "") #, title = primer[p])

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 5)) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

pHLLG

plLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  #geom_jitter(aes(fill = combo), width = 0.4, height = 0, pch = 21, cex = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = c(0, 25)) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 12), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 10)) +
  labs(y = "", x = "")

plLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG/pHLLG/plLLG)



###### Put them together ####
Plot_ClpB | Plot_DnaJ | `Plot_DnaK-fp` | Plot_DnaK1 | Plot_DnaK3 | Plot_GroEL | Plot_GroES | Plot_GrpE | Plot_HrcA | Plot_Hsp20 | Plot_HspA | Plot_HtpG

#plot_spacer() attempt
plot_spacer() + Plot_ClpB + plot_spacer() + Plot_DnaJ + plot_spacer() + `Plot_DnaK-fp` + plot_spacer() + Plot_DnaK1 + plot_spacer() + Plot_DnaK3 + plot_spacer() + Plot_GroEL + plot_spacer() + Plot_GroES + plot_spacer() + Plot_GrpE + plot_spacer() + Plot_HrcA + plot_spacer() + Plot_Hsp20 + plot_spacer() + Plot_HspA + plot_spacer() + Plot_HtpG + plot_spacer() + plot_layout(widths = c(-2.05, 10, -1.8, 10, -1.8, 10, -2, 10, -2, 10, -2, 10, -2, 10, -1.8, 10, -2, 10, -2, 10, -2, 10, -1.8, 10, -6)) + plot_layout(nrow = 1, ncol = 27)

#add back legends and facet headings in powerpoint for a smoother final effect





#######################################################################################
# Figure S6 ---------------------------------------------------------------
##### Load data #####
dat <- read.csv("cq.data.csv")
met <- read.csv("Microcystis growth and metadata.csv")

#### Process data ####
#process metadata
met <- met %>%
  dplyr::select(!c("lake", "strain", "combo_temp", "quality", "combo", "genotype",
                   "Trophic", "temp"))

sdat <- dat %>%
  #remove blanks
  dplyr::filter(!Sample %in% c("", "Blank", "blank", "BLANK")) %>%
  #the group by factors identifying unique treatments, also include any other factors that need to be kept in the output
  group_by(flask, Target, lake, strain, temp, genotype, Trophic, combo) %>%
  #summarise gene expression within each grouping
  summarise(mean.b20.CT.metric = mean(Final.metric...20.as.control, na.rm = TRUE)) %>%
  #join metadata and gene expression data
  left_join(y = met, by = join_by(flask == flask)) %>%
  #add a column that gives the type of gene
  mutate(Type = case_when(Target == "Hep" ~ "Toxin Gene", .default = "HSP Gene"))  %>%
  #remove unecessary columns
  dplyr::select(!c("X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6", "X.7",))  %>%
  dplyr::filter(temp != "20") %>% 
  #remove the reference condition and gene which were used for calculating the gene expression metric
  dplyr::filter(!Target %in% c("", "rpoA"))

#create labels for the temp facets
temp.labels <- c("24\u00B0C", "28\u00B0C", "32\u00B0C")
names(temp.labels) <- c("24", "28", "32")

#create a vector of gene targets to help with plotting
primer <- unique(sdat$Target)

#make data frames to tally highest genotype across facets
sums <- sdat %>%
  group_by(temp, Target, combo) %>%
  summarise(mean = mean(mean.b20.CT.metric))

#correct genotype label 
sdat <- sdat %>%
  mutate(combo = case_when(
    combo == "HNHG" ~ "HL/HG",
    combo == "HNLG" ~ "HL/LG",
    combo == "LNLG" ~ "LL/LG"
  ))

#### Normalization of data ####
###Box cox
#Check data as is
hist(sdat$mean.b20.CT.metric)
shapiro.test(sdat$mean.b20.CT.metric)
normalTest(sdat$mean.b20.CT.metric, method = "da")

#boxcox transformation
b <- boxcox(lm(sdat$mean.b20.CT.metric ~ 1))
lambda <- b$x[which.max(b$y)]
lambda

#check effect of boxcox transformation on data 
hist((sdat$mean.b20.CT.metric^lambda - 1)/lambda)
ggdensity((sdat$mean.b20.CT.metric^lambda - 1)/lambda, y = "density", fill = "grey") 
shapiro.test((sdat$mean.b20.CT.metric^lambda - 1)/lambda)
normalTest((sdat$mean.b20.CT.metric^lambda - 1)/lambda, method = "da")

#Transform data
tdat <- sdat %>%
  mutate(CT.metric = (mean.b20.CT.metric^lambda - 1)/lambda) 

#make sure various variables are in the correct format
tdat$temp <- as.factor(tdat$temp)
tdat$combo <- as.factor(tdat$combo)
tdat$strain <- as.factor(tdat$strain)
tdat$mean.b20.CT.metric <- as.numeric(tdat$mean.b20.CT.metric)

#remove Hep from the general model
hsps <- tdat %>%
  dplyr::filter(!Target %in% c("Hep"))

ggdensity(hsps$CT.metric, fill = "grey")


## Actual plotting ####
###### Primer 1 = ClpB ####
#set target
p <- 1

#make temporary dataframes of the data for set target into subsets
tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

#plot subsets separately to allow axes to vary
p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-7, 5, by=2), 
                     limits=c(-7, 5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 0.889))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-7, 5, by=2), 
                     limits=c(-7, 5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28  

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0, 0.9287))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-7, 5, by=2), 
                     limits=c(-7, 5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32  

#assign a name to the three plots joined together
nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)

Plot_ClpB


###### Primer 2 = DnaJ ####
p <- 2

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 0.93))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 6, by=2), limits=c(-2, 6),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 0.943))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 6, by=2), limits=c(-2, 6),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0, 0.929))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 6, by=2), limits=c(-2, 6),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer 3 = DnaK-fp ####
p <- 3

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-4.5, 4.5, by=1.5), limits=c(-4.5, 4.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-4.5, 4.5, by=1.5), limits=c(-4.5, 4.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0.05, 0.9289))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-4.5, 4.5, by=1.5), limits=c(-4.5, 4.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer 4 = DnaK1 ####
p <- 4

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-9, 9, by=2), 
                     limits=c(-9, 9),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-9, 9, by=2), 
                     limits=c(-9, 9),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-9, 9, by=2), 
                     limits=c(-9, 9),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer 5 = DnaK3 ####
p <- 5

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 4, by=1), 
                     limits=c(-2, 4),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 1))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 4, by=1), 
                     limits=c(-2, 4),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0, 0.93))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 4, by=1), 
                     limits=c(-2, 4),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer 6 = GroEL ####
p <- 6

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0.05, 0.94469))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-6, 6, by=2), 
                     limits=c(-6, 6),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 0.94))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-6, 6, by=2), 
                     limits=c(-6, 6),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0.1, 0.95))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-6, 6, by=2), 
                     limits=c(-6, 6),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer 7 = GroES ####
p <- 7

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 0.99))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 3, by=1), 
                     limits=c(-2, 3),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 0.999))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 3, by=1), 
                     limits=c(-2, 3),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0, 0.925))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-2, 3, by=1), 
                     limits=c(-2, 3),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer 8 = GrpE ####
p <- 8

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0.1, 0.99))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-12, 8, by=4), 
                     limits=c(-12, 8),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 0.95))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-12, 8, by=4), 
                     limits=c(-12, 8),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0.06, 0.99))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-12, 8, by=4), 
                     limits=c(-12, 8),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)

###### Primer  10 = HrcA ####
p <- 10

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 0.94))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-5, 15, by=5), 
                     limits=c(-5, 15.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0.1, 0.94))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-5, 15, by=5), 
                     limits=c(-5, 15.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0, 0.86))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-5, 15, by=5), 
                     limits=c(-5, 15.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer 11 = Hsp20 ####
p <- 11

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0.05, 0.9445))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-12.5, 12.5, by=2.5), 
                     limits=c(-12.5, 12.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0.1, 0.99))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-12.5, 12.5, by=2.5), 
                     limits=c(-12.5, 12.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0.01, 0.98))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-12.5, 12.5, by=2.5), 
                     limits=c(-12.5, 12.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer  12 = HspA ####
p <- 12

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 0.95))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-6.5, 11.5, by=2.5), 
                     limits=c(-6.5, 11.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 0.99))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-6.5, 11.5, by=2.5), 
                     limits=c(-6.5, 11.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$CT.metric, c(0, 0.9))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-6.5, 11.5, by=2.5), 
                     limits=c(-6.5, 11.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Primer  13 = HtpG ####
p <- 13

tmp24 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "24")
tmp28 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "28")
tmp32 <- tdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(temp == "32")

p24 <- tmp24 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp24$CT.metric, c(0, 0.8896))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-5, 7.5, by=2.5), 
                     limits=c(-5, 7.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "") #, title = primer[p])

p24

p28 <- tmp28 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp28$CT.metric, c(0, 0.944))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-5, 7.5, by=2.5), 
                     limits=c(-5, 7.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p28

p32 <- tmp32 %>%
  ggplot(aes(x = combo, y = CT.metric, colour = combo)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  geom_jitter(aes(fill = combo), width = 0.375, height = 0, pch = 21, size = 1, stroke = 2.5, show.legend = FALSE) +
  scale_colour_manual(name="Lake Trophic Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Trophic Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  #coord_cartesian(ylim = quantile(tmp32$mean.b20.CT.metric, c(0.1, 0.8))) +
  theme( axis.text.x = element_blank(), #element_text(angle = 75, vjust = 0.5, hjust=0.35, color="black", size = 10, face = "bold"), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold")) +
  scale_y_continuous(breaks = seq(-5, 7.5, by=2.5), 
                     limits=c(-5, 7.5),
                     labels = label_number(accuracy = 0.1)) +
  labs(y = "", x = "")

p32

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, p24/p28/p32)


###### Put the smaller plots for all the HSPs together ####
Plot_ClpB | Plot_DnaJ | `Plot_DnaK-fp` | Plot_DnaK1 | Plot_DnaK3 | Plot_GroEL | Plot_GroES | Plot_GrpE | Plot_HrcA | Plot_Hsp20 | Plot_HspA | Plot_HtpG

#plot_spacer() attempt
plot_spacer() + Plot_ClpB + plot_spacer() + Plot_DnaJ + plot_spacer() + `Plot_DnaK-fp` + plot_spacer() + Plot_DnaK1 + plot_spacer() + Plot_DnaK3 + plot_spacer() + Plot_GroEL + plot_spacer() + Plot_GroES + plot_spacer() + Plot_GrpE + plot_spacer() + Plot_HrcA + plot_spacer() + Plot_Hsp20 + plot_spacer() + Plot_HspA + plot_spacer() + Plot_HtpG + plot_spacer() + plot_layout(widths = c(-2.05, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -2, 10, -6)) + plot_layout(nrow = 1, ncol = 27)

#add back legends and facet headings in powerpoint for a smoother final effect

tiff("./Final figure files/FigS6_R.tiff", units="px", width=8000, height=4000, res=300)
Plot_ClpB | Plot_DnaJ | `Plot_DnaK-fp` | Plot_DnaK1 | Plot_DnaK3 | Plot_GroEL | Plot_GroES | Plot_GrpE | Plot_HrcA | Plot_Hsp20 | Plot_HspA | Plot_HtpG
dev.off()

#######################################################################################
# Figure S7 ---------------------------------------------------------------
#### Load data #####
dat <- read.csv("cq.data.csv")
met <- read.csv("Microcystis growth and metadata.csv")

#### Process data ####
#process metadata
met <- met %>%
  dplyr::select(!c("lake", "strain", "combo_temp", "quality", "combo", "genotype",
                   "Trophic", "temp"))

sdat <- dat %>%
  #remove blanks
  dplyr::filter(!Sample %in% c("", "Blank", "blank", "BLANK")) %>%
  #the group by factors identifying unique treatments, also include any other factors that need to be kept in the output
  group_by(flask, Target, lake, strain, temp, genotype, Trophic, combo) %>%
  #summarise gene expression within each grouping
  summarise(mean.b20.CT.metric = mean(Final.metric...20.as.control, na.rm = TRUE)) %>%
  #join metadata and gene expression data
  left_join(y = met, by = join_by(flask == flask)) %>%
  #add a column that gives the type of gene
  mutate(Type = case_when(Target == "Hep" ~ "Toxin Gene", .default = "HSP Gene"))  %>%
  #remove unecessary columns
  dplyr::select(!c("X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6", "X.7",))  %>%
  dplyr::filter(temp != "20") %>% 
  #remove the reference condition and gene which were used for calculating the gene expression metric
  dplyr::filter(!Target %in% c("", "rpoA"))

#create labels for the temp facets
temp.labels <- c("24\u00B0C", "28\u00B0C", "32\u00B0C")
names(temp.labels) <- c("24", "28", "32")

#create a vector of gene targets to help with plotting
primer <- unique(sdat$Target)

#make data frames to tally highest genotype across facets
sums <- sdat %>%
  group_by(temp, Target, combo) %>%
  summarise(mean = mean(mean.b20.CT.metric))

#correct genotype label 
sdat <- sdat %>%
  mutate(combo = case_when(
    combo == "HNHG" ~ "HL/HG",
    combo == "HNLG" ~ "HL/LG",
    combo == "LNLG" ~ "LL/LG"
  ))

#### Actual plotting ####
p <- 9

tmpHLHG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/HG")
tmpHLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "HL/LG")
tmpLLLG <- sdat %>% 
  dplyr::filter(Target == primer[p]) %>% 
  dplyr::filter(combo == "LL/LG")

pHLHG <- tmpHLHG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo, group = temp)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  scale_colour_manual(name="Lake Nutrient Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Nutrient Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLHG$mean.b20.CT.metric, c(0, 0.925))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 16), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 16),
         axis.title = element_text(colour = "black", face = "bold", size = 18)) +
  labs(y = "", x = "") #expression(paste("Relative (to RpoA) differential (as compared to 20 ", degree~C, ") gene expression"))

pHLHG

pHLLG <- tmpHLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo, group = temp)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  scale_colour_manual(name="Lake Nutrient Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Nutrient Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpHLLG$mean.b20.CT.metric, c(0, 0.924))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 16), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 16),
         axis.title = element_text(colour = "black", face = "bold", size = 18)) +
  labs(y = "", x = "") #expression(paste("Temperature ", degree~C, ")"))

pHLLG

pLLLG <- tmpLLLG %>%
  ggplot(aes(x = as.factor(temp), y = mean.b20.CT.metric, colour = combo, group = temp)) +
  theme_bw() +
  geom_boxplot(aes(colour = combo), outlier.shape = NA, cex = 1) +
  scale_colour_manual(name="Lake Nutrient Status/Genotype",
                      values=c("HL/HG"="green","LL/LG"="blue","HL/LG"="#56B4E9"))+
  scale_fill_manual(name="Lake Nutrient Status/Genotype",
                    values=c("HL/HG"="#99FF66","LL/LG"="#3399FF","HL/LG"="lightblue")) +
  coord_cartesian(ylim = quantile(tmpLLLG$mean.b20.CT.metric, c(0, 0.8899))) +
  theme( axis.text.x = element_text(color="black", face = "bold", size = 16), 
         legend.position = "",
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         axis.text.y = element_text(color="black", face = "bold", size = 16),
         axis.title = element_text(colour = "black", face = "bold", size = 12)) +
  labs(y = "", x = "") 

pLLLG

nam <- paste("Plot", primer[p], sep = "_")
assign(nam, pHLHG + pHLLG + pLLLG)

Plot_Hep

tiff("./Final figure files/FigS7_R.tiff", units="px", width=2500, height=1500, res=300)
Plot_Hep
dev.off()
#add back legends and facet headings in powerpoint for a smoother final effect


#######################################################################################
# Figure S8 --------------------------------------------------------------
#expression data
tdat <- read.csv("tdat.csv")

#growth data
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# This removes CR19-01, which is the sole LN/HG

#boxcox transformation of data for better normality (precident and testing done previously)
b <- boxcox(lm(data$Growth +1 ~ 1))
lambda <- b$x[which.max(b$y)]
lambda

#Transform data
data <- data %>%
  #while transformation doesn't make a huge difference the histogram and density looks a bit better and it is consistent to use this as used it earlier
  mutate(Growth_transformed = ((Growth + 1)^lambda - 1)/lambda) %>% 
  #also remove unecessary columns of the data
  dplyr::select(flask, strain, Lake, combo_temp, combo, genotype, Trophic, temp, Crash_Numeric, Crash, Growth_transformed, Phase) %>% 
  #and widen data so each Phase has its own growth data column
  dcast(flask + strain + temp + genotype + Trophic + combo + Crash_Numeric ~ Phase, 
        value.var = "Growth_transformed") %>% 
  #rename these to make it clear they are the transformed values
  rename(Exp_transformed = Exp, wk2_transformed = wk2, wk3_transformed = wk3, wk4_transformed = wk4)


#'widen' expression data
tdat_wide <- dcast(data = tdat, flask + combo + temp + genotype + Trophic + lake + strain + Crash_Numeric + ExpPhase + wk2 + wk3 + wk4 ~ Target, 
                   value.var = "CT.metric") %>% 
  left_join(y = data, by = join_by(flask, strain, temp, genotype, Trophic, combo, Crash_Numeric))


full_post_hoc <- lda(combo~cbind(ClpB, DnaJ, `DnaK-fp`, DnaK1, DnaK3, GroEL, GroES, GrpE, Hep, HrcA, Hsp20, HspA, HtpG, 
                                 Exp_transformed, wk2_transformed, wk3_transformed, wk4_transformed), 
                     data = tdat_wide,
                     CV = F)
full_post_hoc


full_plot_lda <- data.frame(tdat_wide[, "combo"], lda = predict(full_post_hoc)$x)
colnames(full_plot_lda) <- c("combo", "lda.LD1", "lda.LD2" )


#biplot with just the expression variables
full_arrows <- data.frame("x0" = 0, 
                          "y0" = 0,
                          "x1" = 4 * coef(full_post_hoc)[,1],
                          "y1" = 4 * coef(full_post_hoc)[,2],
                          "labs" = c("ClpB", "DnaJ", "DnaK-fp", "DnaK1", "DnaK3", "GroEL", 
                                     "GroES", "GrpE", "Hep", "HrcA", "Hsp20", "HspA", "HtpG", "Week 1", 
                                     "Week 2", "Week 3", "Week 4"))

FigS8 <- ggplot(full_plot_lda) + 
  geom_point(aes(x = lda.LD1, y = lda.LD2, fill = combo), size = 5, pch = 21, alpha = 0.75) +
  scale_fill_manual(name="Bacterial Type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),
                     values=c("HNHG"="green","HNLG"="#56B4E9","LNLG"="blue")) +
  labs(fill = "Bacterial Type",
       x = "Linear Discriminate Analysis Axis 1 (74.7%)",
       y = "Linear Discriminate Analysis Axis 2 (25.3%)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.position = "bottom") +
  geom_segment(data = full_arrows[!full_arrows$labs %in% c("Week 1", "Week 2", "Week 3", "Week 4"), ], 
               aes(x = x0, xend = x1, y = y0, yend = y1),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, alpha = 0.75, colour = "black") +
  #geom_text_repel(data = full_arrows[!full_arrows$labs %in% c("Week 1", "Week 2", "Week 3", "Week 4"), ], 
  #          aes(x = x1 + (x1/10), y = y1, label = labs), cex = 5, colour = "black", fontface = "bold",
  #          max.overlaps = 14, segment.size = 1, min.segment.length = Inf, 
  #          segment.curvature = -0.1, point.size = 5) +
  geom_segment(data = full_arrows[full_arrows$labs %in% c("Week 1", "Week 2", "Week 3", "Week 4"), ], 
               aes(x = x0, xend = x1, y = y0, yend = y1),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, alpha = 0.75, colour = "darkgrey") +
  coord_cartesian(ylim = c(-4, 4), xlim = c(-4.2, 3))

FigS8

tiff("./Final figure files/FigS8_R.tiff", units="px", width=2500, height=2000, res=300)
FigS8
dev.off()
  















