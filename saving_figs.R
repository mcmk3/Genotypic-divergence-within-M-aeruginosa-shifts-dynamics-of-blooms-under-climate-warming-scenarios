#######################################################################################
# Setup -------------------------------------------------------------------
setwd("C:/Users/mirte/Documents/UCSD/Projects/Microcystis Warming")

rm(list = ls())

library("dplyr") #For data manipulation
library("reshape2") #for data manipulation
library("MASS") #Needed for LDA

library("ggplot2") #For plotting
library("ggnewscale") #To allow two colour scales in Figure 2
library("ggpubr") #For ggdensity plotting
library("patchwork") #For arranging facets of Figure 4, S3, S4 & S5
library("scales") #For better axis labeling in Figure S4
library("ggrepel") #For labels in Figure 5
library("ggtree") #For maps
library("tiff")

#######################################################################################
# Figure 1 ----------------------------------------------------------------
tree<-read.tree("RAxML_bipartitions.Testing_Aug_16_2024_noCR1901_10kboots.newick")#,node.label = 'support')
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

tiff("./Final figure files/Fig1A.tiff", units="px", width=7500, height=9000, res=300)
Fig_1
dev.off()

#######################################################################################
# Figure 2 ----------------------------------------------------------------
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# Removing CR19-01, which was the sole LN/HG
data_exp<-data[(data$Phase=="Exp"),]#Sorting for exponential phase growth

#Remove outlier
data_exp_remO <-data_exp[(data_exp$Outlier=="n"),]# Removing an LN/LG outlier from the main dataset, but it is added back in below 

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
  ylab("Rate of Exponential Growth")+
  guides(size="none",alpha="none") +
  #FACET
  facet_grid(~combo, labeller=as_labeller(c('HNHG'="HL/HG",'HNLG'="HL/LG",'LNLG'="LL/LG")))

Fig_2

tiff("./Final figure files/Fig2.tiff", units="px", width=2000, height=1250, res=300)
Fig_2
dev.off()

#######################################################################################
# Figure 3 ----------------------------------------------------------------
data<-read.csv("Sara_Revised_New_FixedSM.csv")
data<-data[(data$Keep=="y"),]# This removes CR19-01, which is the sole LN/HG

#Plotting
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
  geom_jitter(aes(fill=NegativePositiveGrowth, colour=combo),
              size=2.5, stroke=2.5, 
              shape=21,width=0.36,height=0.0)+
  theme_bw()+
  theme(strip.text.y=element_text(size=10),strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),axis.text=element_text(size=16),
        axis.title=element_text(size=16), legend.position = "bottom")+
  labs(x=expression(paste("Temperature (", degree~C, ")")),y="Growth Rate") +
  scale_fill_manual(name="Growth Rate",values=c("negative"="black","positive"="white"))+
  scale_color_manual(name="Bacterial Type",labels=c("HNHG"="HL/HG","HNLG"="HL/LG","LNLG"="LL/LG"),values=c("HNHG"="green","HNLG"="#56B4E9","LNHG"="black","LNLG"="blue"))+
  guides(size="none",alpha="none")+
  facet_grid(Phase~combo,labeller = labeller(Phase = week_labeller, combo = combo_labeller))

Fig_3

tiff("./Final figure files/Fig3.tiff", units="px", width=4500, height=3000, res=300)
Fig_3
dev.off()

#######################################################################################
# Figure 4 ----------------------------------------------------------------
##### Load data #####
tdat <- read.csv("tdat.csv")

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
Fig_4 <- Plot_ClpB | Plot_DnaJ | `Plot_DnaK-fp` | Plot_DnaK1 | Plot_DnaK3 | Plot_GroEL | Plot_GroES | Plot_GrpE | Plot_HrcA | Plot_Hsp20 | Plot_HspA | Plot_HtpG

Fig_4

tiff("./Final figure files/Fig4_R.tiff", units="px", width=5500, height=3000, res=300)
Fig_4
dev.off()

#add back legends and facet headings in powerpoint for a smoother final effect

#######################################################################################
# Figure 5 ------------------------------------------------------------------
tdat <- read.csv("tdat.csv")

#first need to 'widen' data
tdat_wide <- dcast(data = tdat, flask + combo + temp + genotype + Trophic + lake + strain + Crash_Numeric + ExpPhase + wk2 + wk3 + wk4 ~ Target, 
                   value.var = "CT.metric")


full_post_hoc <- lda(combo~cbind(ClpB, DnaJ, `DnaK-fp`, DnaK1, DnaK3, GroEL, 
                                 GroES, GrpE, Hep, HrcA, Hsp20, HspA, HtpG, ExpPhase, wk2, wk3, wk4), 
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
       x = "Linear Discriminate Analysis Axis 1 \n(Proportion of trace = 0.69)",
       y = "Linear Discriminate Analysis Axis 2 \n(Proportion of trace = 0.31)") +
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
       x = "Linear Discriminate Analysis Axis 1 (68.7%)",
       y = "Linear Discriminate Analysis Axis 2 (31.2%)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16))

#unlabeled version for powerpoint improved labels
Fig_5 <- ggplot(full_plot_lda) + 
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
       x = "Linear Discriminate Analysis Axis 1 (68.7%)",
       y = "Linear Discriminate Analysis Axis 2 (31.2%)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16),
        plot.title=element_text(hjust=0.5),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.position = "bottom")

Fig_5

tiff("./Final figure files/Fig5_R.tiff", units="px", width=2500, height=2000, res=300)
Fig_5
dev.off()


