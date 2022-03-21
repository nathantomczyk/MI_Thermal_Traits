###################################
### do functional feeding groups have different temperature sensitivities?
###################################


library(ape)
library(ggplot2)
library(nlme)
library(phytools)
library(Rmisc)
library(reshape2)
library(multcomp)
library(MuMIn)
library(cowplot)
library(ggtree)
library(caper)
library(dplyr)
library(tidytree)
library(RColorBrewer)
library(emmeans)
library(lme4)
library(ggsignif)


## EU_data
column.names<-c("Taxon","EU","Grazers","Miners","Xylophagous","Shredders",
                "Gatherers/ collectors", "Active Filter",
                "Passive filter","Predators","Parasites",
                "Other")


ffg_EU<-read.csv("./data/EU_traits_data.csv",col.names = column.names)
ffg_EU[is.na(ffg_EU)] <- 0

ffg_EU.melt<-melt(ffg_EU)
ffg_EU.melt.no.zeros<-ffg_EU.melt[ffg_EU.melt$value!=0,]

ffg.drop<-c("Miners","Xylophagous","Other","Parasites")
ffg_EU.melt.no.zeros<-ffg_EU.melt.no.zeros[ffg_EU.melt.no.zeros$variable %in% ffg.drop==FALSE,]
ffg_EU.melt.no.zeros<-ffg_EU.melt.no.zeros[is.na(ffg_EU.melt.no.zeros$variable)==FALSE,]
ffg_EU.melt.no.zeros$variable <- recode_factor(ffg_EU.melt.no.zeros$variable, Passive.filter = "Filterer", 
                                            Active.Filter = "Filterer")

ffg_EU.melt.no.zeros$variable<-recode_factor(ffg_EU.melt.no.zeros$variable,
                                          Gatherers..collectors="Collector-gatherer",
                                          Predators="Predator",
                                          Shredders="Shredder",
                                          Grazers="Herbivore")

clear.best<-ffg_EU.melt.no.zeros[ffg_EU.melt.no.zeros$value>5,]
generalists<-ffg_EU.melt.no.zeros[ffg_EU.melt.no.zeros$Taxon %in% clear.best$Taxon==FALSE,]
generalists$value<-10
generalists$variable<-"Generalist"
generalists<-unique(generalists)

clean_EU_FFG<-rbind(clear.best,generalists)


T.val.EU<-read.csv("./data/EU_continous_temperature.csv")
T.val.EU[is.na(T.val.EU)]<-0
T.val.EU<-T.val.EU[T.val.EU$temperature>0,]

ffg_EU.melt.no.zeros$GENUS<-sub(" .*", "", ffg_EU.melt.no.zeros$Taxon)


merge.EU.data.val<-merge(T.val.EU,clean_EU_FFG,by="Taxon")



################ 
taxa<-levels(as.factor(as.character(merge.EU.data.val$Taxon))) ### 332 taxa


EU.model<-aov(temperature~variable,data=merge.EU.data.val)
summary(EU.model)
coefficients(EU.model)
TukeyHSD(EU.model)
eu.tukey<-glht(EU.model, linfct=mcp(variable="Tukey"))
letter.label<-cld(eu.tukey)
labels<-c("c","cd","a","ab","bcd","d")
r.squaredGLMM(EU.model)

group.means<-summarySE(merge.EU.data.val,measurevar = "temperature",groupvars="variable")

annotation_layer<-data.frame(x=1:6+0.2,y=12.5,labels)


EU_data_new<-ggplot(merge.EU.data.val,aes(x=variable,y=temperature))+
  geom_boxplot(fill="navyblue",alpha=0.5)+
  #geom_bar(stat="identity",fill="navyblue",alpha=0.5)+
  theme_classic()+xlab("Functional Feeding Group")+
  theme(axis.text.x = element_text(angle = 90))+
  ylab(expression("Temperature preference "^o*"C"))+
  xlab("")+
  theme(text =element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(data=annotation_layer,aes(x=x,y=y,label=labels),size=6)
EU_data_new


############################ Done with EU data. On to US data

ffg<-read.csv("./data/FreshwaterBioTraits_20100927.csv")
ffg$species2<- gsub(" ", "_", ffg$TAXON)
useful_names<-c("Primary functional feeding group abbreviated",
                "Primary feeding mode")

ffg_subset<-ffg[ffg$TRAITS_NAME %in% useful_names,]

levels(as.factor(as.character(ffg_subset$VALUE_TEXT)))


### hand coding things to harmonize the data
ffg_subset<-ffg_subset[!is.na(ffg_subset$VALUE_TEXT),]

ffg_subset[ffg_subset$VALUE_TEXT=="CF","VALUE_TEXT"]<-"Filterer"
ffg_subset[ffg_subset$VALUE_TEXT=="Collector-filterer","VALUE_TEXT"]<-"Filterer"
ffg_subset[ffg_subset$VALUE_TEXT=="CG","VALUE_TEXT"]<-"Collector-gatherer"
ffg_subset[ffg_subset$VALUE_TEXT=="HB","VALUE_TEXT"]<-"Grazer" 
ffg_subset[ffg_subset$VALUE_TEXT=="Piercer herbivore","VALUE_TEXT"]<-"Herbivore"
ffg_subset[ffg_subset$VALUE_TEXT=="Scraper/grazer","VALUE_TEXT"]<-"Herbivore"


groups_to_keep<-c("Shredder","Herbivore","Predator","Collector-gatherer","Filterer")


ffg_subset_clean<-ffg_subset[ffg_subset$VALUE_TEXT %in% groups_to_keep,]

####### identifying which taxa are generalists - here generalists are taxa with more than one identified best FFG

ffg_subset_clean<-ffg_subset_clean[ffg_subset_clean$TAXON!="",]
ffg_subset_clean<-ffg_subset_clean[is.na(ffg_subset_clean$TAXON)==FALSE,]
ffg_subset_clean$numeric.ffg<-NA
ffg_subset_clean[ffg_subset_clean$VALUE_TEXT=="Filterer","numeric.ffg"]<-0
ffg_subset_clean[ffg_subset_clean$VALUE_TEXT=="Collector-gatherer","numeric.ffg"]<-1
ffg_subset_clean[ffg_subset_clean$VALUE_TEXT=="Herbivore","numeric.ffg"]<-2
ffg_subset_clean[ffg_subset_clean$VALUE_TEXT=="Parasite","numeric.ffg"]<-3
ffg_subset_clean[ffg_subset_clean$VALUE_TEXT=="Predator","numeric.ffg"]<-4
ffg_subset_clean[ffg_subset_clean$VALUE_TEXT=="Shredder","numeric.ffg"]<-5

ffg.summary<-summarySE(ffg_subset_clean,measurevar="numeric.ffg",groupvars="TAXON")
ffg.summary[is.na(ffg.summary)] <- 0
reduced.epa.data<-data.frame(TAXON=ffg.summary$TAXON,sd=ffg.summary$sd,numeric.ffg=ffg.summary$numeric.ffg,ffg=NA)

reduced.epa.data[reduced.epa.data$numeric.ffg==0,"ffg"]<-"Filterer"
reduced.epa.data[reduced.epa.data$numeric.ffg==1,"ffg"]<-"Collector-gatherer"
reduced.epa.data[reduced.epa.data$numeric.ffg==2,"ffg"]<-"Herbivore"
reduced.epa.data[reduced.epa.data$numeric.ffg==3,"ffg"]<-"Parasite"
reduced.epa.data[reduced.epa.data$numeric.ffg==4,"ffg"]<-"Predator"
reduced.epa.data[reduced.epa.data$numeric.ffg==5,"ffg"]<-"Shredder"

reduced.epa.data[reduced.epa.data$sd>0,"ffg"]<-"Generalist" ### 



##############################


traits<-read.csv("./data/FreshwaterBioTraits_20100927.csv")

###thermal preferance categories

topt<-traits[traits$TRAITS_NAME=="Thermal optima value",]

maxtemp<-traits[traits$TRAITS_NAME=="Maximum temperature reported",]

ffg.epa<-data.frame(TAXON=reduced.epa.data$TAXON,feeding_group=as.factor(as.character(reduced.epa.data$ffg)))

max.temp.merged<-merge(maxtemp,reduced.epa.data[reduced.epa.data$ffg!="Uncertain",],by="TAXON")
topt.merged<-merge(topt,reduced.epa.data[reduced.epa.data$ffg!="Uncertain",],by="TAXON")

t.opt.merged.cleaned<-data.frame(temperature=topt.merged$VALUE_NUMBER,Taxon=topt.merged$TAXON,
                                 citation=topt.merged$CITATION_COMPLETE,ffg=topt.merged$ffg)

t.opt.merged.cleaned<-unique(t.opt.merged.cleaned)

#### modeling optimum temperatures

optimum_temp_model<-lmer(temperature~ffg+(1|Taxon)+(1|citation),data=t.opt.merged.cleaned)
anova(optimum_temp_model)
summary(optimum_temp_model)
r.squaredGLMM(optimum_temp_model)
optimum_temp_model.contrast<-emmeans(optimum_temp_model, "ffg")
contrast(optimum_temp_model.contrast,"tukey")


max.temp.clean<-data.frame(temperature=max.temp.merged$VALUE_NUMBER,Taxon=max.temp.merged$TAXON,
                           citation=max.temp.merged$CITATION_COMPLETE,ffg=max.temp.merged$ffg)

max.temp.clean<-unique(max.temp.clean)
max.temp.model<-lmer(temperature~ffg+(1|Taxon)+(1|citation),data=max.temp.clean)
anova(max.temp.model)
summary(max.temp.model)
max.temp.model.contrast<-emmeans(max.temp.model, "ffg")
contrast(max.temp.model.contrast,"tukey")
r.squaredGLMM(max.temp.model)

opt.plot<-data.frame(optimum_temp_model.contrast)

t.opt.merged.cleaned$ffg <- factor(t.opt.merged.cleaned$ffg, levels=c("Collector-gatherer","Predator","Shredder","Herbivore",
                                                "Filterer","Generalist"))

contrast.letters<-c("ab","b","a","ab","ab","ab")

annotation.layer<-data.frame(x=1:6+0.2,y=7,labels=contrast.letters)

optimum_temp_bar<-ggplot(t.opt.merged.cleaned,aes(y=temperature,x=ffg))+
  geom_boxplot(fill="navyblue",alpha=0.5)+
  #geom_bar(stat="identity",fill="navyblue",alpha=0.5)+
  theme_classic()+ylab(expression("Optimimum Temperature "^o*"C"))+
  xlab("")+
  #geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0.00001)+
  theme(text =element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(data=annotation.layer,aes(x=x,y=y,label=labels),size=6)
optimum_temp_bar


max.plot<-data.frame(max.temp.model.contrast)

max.temp.clean$ffg <- factor(max.temp.clean$ffg, levels=c("Collector-gatherer","Predator","Shredder","Herbivore",
                                              "Filterer","Generalist"))


contrast.letters.2<-c("ab","ab","a","ab","b","ab")

annotation.layer<-data.frame(x=1:6+0.2,y=5,labels=contrast.letters.2)

max_temp_bar<-ggplot(max.temp.clean,aes(y=temperature,x=ffg))+
  geom_boxplot(fill="navyblue",alpha=0.5)+
  theme_classic()+ylab(expression("Maximum Observed Temperature "^o*"C"))+
  xlab("")+
    theme(text =element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(data=annotation.layer,aes(x=x,y=y,label=labels),size=6)
max_temp_bar


### Reading in tolerance data from Chown 2015, this data was gathered in manipulative lab experiments
tol<-read.csv("./data/bug_tolerance_data.csv")

#Chown Upper thermal tolerance in aquatic insects

#tol$species2<- gsub(" ", "_", tol$Species)
tol$GENUS<-sub(" .*", "", tol$Species)

tol<-tol[tol$Stage=="Immatures",]

tol$Taxon<-tol$Species
tol$TAXON<-tol$Species

#### some taxa are only identified to the level to genus ### separate these one

############## Lab data #################################### Need to clean up further


eu_merge<-merge(tol,clean_EU_FFG,by="Taxon")

eu_merge_tidy<-data.frame(Taxon=eu_merge$Taxon,ffg=eu_merge$variable,ctmax=eu_merge$Ctmax,ult=eu_merge$ULT..C.,citation=eu_merge$Reference)

us_merge<-merge(tol,reduced.epa.data,by="TAXON")

us_merge_tidy<-data.frame(Taxon=us_merge$TAXON,ffg=us_merge$ffg,ctmax=us_merge$Ctmax,ult=us_merge$ULT..C.,citation=us_merge$Reference)

all_tol_data<-rbind(us_merge_tidy,eu_merge_tidy)

all_tol_data<-unique(all_tol_data)



left_out<-tol[tol$Taxon %in% all_tol_data$Taxon==FALSE,]
#left_out$GENUS<-sub(" .*", "", left_out$Species)

reduced.epa.data$GENUS<-sub(" .*", "", reduced.epa.data$TAXON)

genus_data<-reduced.epa.data[reduced.epa.data$GENUS %in% left_out$GENUS,]

genus_summary<-summarySE(genus_data,measurevar = "numeric.ffg",groupvars = c("GENUS","ffg"))

clear_genuses<-c("Anopheles","Celithemis","Centroptilum","Chimarra","Chironomus","Clioperla","Dolophilodes",
                 "Epitheca","Gumaga","Limnephilus","Macromia","Prostoia","Pteronarcella","Pteronarcys","Pycnopsyche",
                 "Taeniopteryx")
genus_summary_good<-genus_summary[genus_summary$GENUS %in% clear_genuses,]

us_genus_good<-merge(left_out,genus_summary_good,by="GENUS")

extra_us_genus<-data.frame(Taxon=us_genus_good$GENUS,ffg=us_genus_good$ffg,ctmax=us_genus_good$Ctmax,
           ult=us_genus_good$ULT..C.,citation=us_genus_good$Reference)

all_tol_data2<-rbind(all_tol_data,extra_us_genus)

########## models

ctmax_model<-lmer(ctmax~ffg+(1|citation),data=all_tol_data2)
anova(ctmax_model)
summary(ctmax_model)
ctmax.model.contrast<-emmeans(ctmax_model, "ffg")
contrast(ctmax.model.contrast,"tukey")
r.squaredGLMM(ctmax_model)
ctmax.model.contrast<-data.frame(ctmax.model.contrast)

ct.max.plot<-all_tol_data2[is.na(all_tol_data2$ctmax)==FALSE,]

ct.max.plot$ffg <- factor(ct.max.plot$ffg, levels=c("Predator","Shredder",
                                              "Filterer","Generalist"))

labeling_letters3<-c("b","ab","ab","a")

annotation.layer<-data.frame(x=1:4+0.2,y=26,labels=labeling_letters3)

ctmax.plot<-ggplot(ct.max.plot,aes(x=ffg,y=ctmax))+
  geom_boxplot(fill="navyblue",alpha=0.5)+
  theme_classic()+ylab(expression("Critical Thermal Maximum "^o*"C"))+
  xlab("Functional Feeding Group")+
  theme(text =element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(data=annotation.layer,aes(x=x,y=y,label=labels),size=6)
ctmax.plot

# ult

ult_model<-lmer(ult~ffg+(1|citation),data=all_tol_data2)
anova(ult_model)
summary(ult_model)
ult.model.contrast<-emmeans(ult_model, "ffg")
contrast(ult.model.contrast,"tukey")
r.squaredGLMM(ult_model)
ult.model.contrast<-data.frame(ult.model.contrast)

ult.plot.data<-all_tol_data2[is.na(all_tol_data2$ult)==FALSE,]

ult.plot.data$ffg <- factor(ult.plot.data$ffg, levels=c("Collector-gatherer","Predator","Shredder","Herbivore",

                                                                  
                                                                                                                "Filterer","Generalist"))
labeling_letters4<-c("ab","b","a","ab","ab","a")
annotation.layer<-data.frame(x=1:6+0.2,y=22,labels=labeling_letters4)

ult.plot<-ggplot(ult.plot.data,aes(x=ffg,y=ult))+
  geom_boxplot(fill="navyblue",alpha=0.5)+
  theme_classic()+ylab(expression("Upper Lethal Temperature "^o*"C"))+
  xlab("Functional Feeding Group")+
  theme(text =element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(data=annotation.layer,aes(x=x,y=y,label=labels),size=6)
ult.plot




tiff(filename="./figures/temperature_ffg_plots.tiff",units="in",res=800,width=10,height=15,compression="lzw")
plot_grid(optimum_temp_bar,EU_data_new,max_temp_bar,ult.plot,ctmax.plot,labels="AUTO",ncol=2,label_x=0.93,label_y=1.0,label_size = 20)
dev.off()


tiff(filename="./figures/topt_poster.tiff",units="in",res=900,width=10,height=6,compression="lzw")
optimum_temp_bar
dev.off()

tiff(filename="./figures/tmax_poster.tiff",units="in",res=900,width=10,height=6,compression="lzw")
max_temp_bar
dev.off()

tiff(filename="./figures/tpref_poster.tiff",units="in",res=900,width=10,height=6,compression="lzw")
EU_data_new
dev.off()

tiff(filename="./figures/ctmax_poster.tiff",units="in",res=900,width=10,height=6,compression="lzw")
ctmax.plot
dev.off()

tiff(filename="./figures/ult_poster.tiff",units="in",res=900,width=10,height=6,compression="lzw")
ult.plot
dev.off()



#################### Phylo analysis #### first is checking how much phylo changes this
library(caper)

tree<-read.tree("./data/species_level_tree.nwk")
#Chesters, Douglas (2016), Data from: Construction of a species-level tree of life for the insects and utility in taxonomic profiling, Dryad, Dataset, https://doi.org/10.5061/dryad.27114
######## Trimming this down to genus level because that is where most of the trait data exists
tips<-tree$tip.label
# genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
# 
# 
# ## now drop all but one of each
# ii<-sapply(genera,function(x,y) grep(x,y)[1],y=tips)
# tree2<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
# tree2$tip.label<-sapply(strsplit(tree2$tip.label,"_"),function(x) x[1])
# str(tree2)

merge.EU.data.val$tip.label<-gsub(" ","_",merge.EU.data.val$Taxon)
merge.EU.data.val$match<-match(merge.EU.data.val$tip.label,tree$tip.label)
clear.best2<-merge.EU.data.val[is.na(merge.EU.data.val$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(clear.best2$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

clear.best.ordered<-merge(trait.ordered,clear.best2,by="tip.label",sort=FALSE)



no_phylo<-lm(temperature~variable,data=clear.best.ordered)
anova(no_phylo)
summary(no_phylo)
r.squaredLR(no_phylo)
AICc(no_phylo)

phylo_model<-gls(temperature~1,data=clear.best.ordered,correlation = corPagel(0.5, pruned.tree))
summary(phylo_model)
r.squaredLR(phylo_model)
AICc(phylo_model)

ffg_phylo_model<-gls(temperature~variable,data=clear.best.ordered,correlation = corPagel(0.5, pruned.tree))
summary(ffg_phylo_model)
anova(ffg_phylo_model)
r.squaredLR(ffg_phylo_model)
AICc(ffg_phylo_model)


T.val.EU$tip.label<-gsub(" ","_",T.val.EU$Taxon)
T.val.EU$match<-match(T.val.EU$tip.label,tree$tip.label)

tree_data<-T.val.EU[is.na(T.val.EU$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data$tip.label, tree$tip.label)])






#topt$matches<-match(topt$GENUS, tree2$tip.label)

## reordering traits for plot

#merge(plot.data.2,plot.data,by="GENUS",sort=FALSE)

####### phylo signature in TPref

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data,by="tip.label",sort=FALSE)

pruned.tree$trait<-trait.ordered2$temperature

trait.values.test<-trait.ordered2$temperature

countSimmap(pruned.tree,trait.ordered2$temperature)


n.tips<-length(trait.values.test)

for (i in 2:(length(trait.values.test))){
  trait.values.test[n.tips+i-1]<-(trait.values.test[i]+trait.values.test[i+1])/2
  
  
}


t.pref.phylo<-ggtree(pruned.tree,layout="rectangular",aes(color=trait.values.test))+ 
  geom_tippoint(,size=2, alpha=.75,show.legend=FALSE)+
  #scale_color_brewer(palette="RdYlBu")+
  scale_color_gradientn(colours=c("blue3","yellow","red3"))+
  labs(fill="Temperature")+
  theme(legend.position = c(0.12,0.9),legend.text=element_text(size=15))+
  #geom_tiplab(size=1)+
  geom_cladelab(node=148, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=181, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=190, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=235, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=139, label="Megaloptera", 
                offset = 0.03, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=137, label="Hymenoptera", 
                offset = 0.01, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=96, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=207, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=201, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')
t.pref.phylo

tiff('./figures/temperature_preferece.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
t.pref.phylo
dev.off()

t.pref.phylo<-ggtree(pruned.tree,layout="rectangular",aes(color=trait.values.test))+ 
  geom_tippoint(,size=2, alpha=.75,show.legend=FALSE)+
  #scale_color_brewer(palette="RdYlBu")+
  scale_color_gradientn(colours=c("blue3","yellow","red3"))+
  labs(fill="Temperature")+
  theme(legend.position = c(0.12,0.9),legend.text=element_text(size=15))+
  geom_tiplab(size=1)+
  geom_cladelab(node=148, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=181, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=190, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=235, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=139, label="Megaloptera", 
                offset = 0.03, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=137, label="Hymenoptera", 
                offset = 0.01, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=96, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=207, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=201, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')
t.pref.phylo

tiff('./figures/names_temperature_preferece.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
t.pref.phylo
dev.off()



tpref_model<-gls(temperature~1,data=trait.ordered2,correlation = corPagel(0.5, pruned.tree))
summary(tpref_model)
anova(tpref_model)
r.squaredLR(tpref_model)

###################### phylogenetic signal in TOpt

t.opt.signal<-topt
shred$tip.label<-gsub(" ","_",shred$Taxon)
shred$match<-match(shred$tip.label,tree$tip.label)

tree_data_shred<-shred[is.na(shred$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data_shred$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data_shred,by="tip.label",sort=FALSE)


trait.values.test<-trait.ordered2$value

n.tips<-length(trait.values.test)

for (i in 2:(length(trait.values.test))){
  trait.values.test[n.tips+i-1]<-(trait.values.test[i]+trait.values.test[i+1])/2
  
  
}

############## Main text figure
clean_EU_FFG$tip.label<-gsub(" ","_",clean_EU_FFG$Taxon)

  
clean_EU_FFG$match<-match(clean_EU_FFG$tip.label,tree$tip.label)

tree_data_main<-clean_EU_FFG[is.na(clean_EU_FFG$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data_main$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data_main,by="tip.label",sort=FALSE)


cat.traits<-c(paste(trait.ordered2$variable),rep(NA,389))

cbPalette <- c( "#E69F00", "#56B4E9","#999999", "#009E73", "#F0E442", "#0072B2", "#000000")

main.ffg.phylo.plot<-ggtree(pruned.tree,layout="rectangular",aes(color=cat.traits))+ 
  geom_tippoint(size=2, alpha=.75,show.legend=FALSE)+
  scale_color_manual(values = cbPalette)+
  theme(legend.position = c(0.12,0.9),legend.title = element_blank(),legend.text=element_text(size=15))+
  #geom_tiplab(size=1)+
  guides(color = guide_legend(override.aes = list(size = 3) ) )+
  geom_cladelab(node=400, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=458, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=490, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=597, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=393, label="Megaloptera", 
                offset = 0.042, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=595, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=555, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=531, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')
main.ffg.phylo.plot



tiff('./figures/main_figure_FFG.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
main.ffg.phylo.plot
dev.off()


############## shredders
ffg_EU.melt<-ffg_EU.melt[ffg_EU.melt$Taxon %in% clean_EU_FFG$Taxon,]

ffg_EU.melt$variable <- recode_factor(ffg_EU.melt$variable, Passive.filter = "Filterer", 
                                            Active.Filter = "Filterer")


shred<-ffg_EU.melt[ffg_EU.melt$variable=="Shredders",]
shred$tip.label<-gsub(" ","_",shred$Taxon)
shred$match<-match(shred$tip.label,tree$tip.label)

tree_data_shred<-shred[is.na(shred$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data_shred$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data_shred,by="tip.label",sort=FALSE)


trait.values.test<-as.numeric(c(paste(trait.ordered2$value),rep(NA,389)))


shredder.plot<-ggtree(pruned.tree,layout="rectangular",aes(color=trait.values.test))+ 
  geom_tippoint(,size=2, alpha=.75,show.legend=FALSE)+
  #scale_color_brewer(palette="RdYlBu")+
  scale_color_gradientn(colours=c("gray70","blue"))+
  theme(legend.position = c(0.12,0.9),legend.title = element_blank())+
  geom_tiplab(size=1)+
  geom_cladelab(node=400, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=458, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=490, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=597, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=393, label="Megaloptera", 
                offset = 0.042, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=595, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=555, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=531, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')

shredder.plot

tpref.signal<-phylosig(pruned.tree,trait.ordered2$value,method="lambda",test=TRUE,)
tpref.signal

shredder_model<-gls(value~1,data=trait.ordered2,correlation = corPagel(0.5, pruned.tree))
summary(shredder_model)
anova(shredder_model)
r.squaredLR(shredder_model)


tiff('./figures/shredders.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
shredder.plot
dev.off()

####################### Herbivore


herb<-ffg_EU.melt[ffg_EU.melt$variable=="Grazers",]
herb$tip.label<-gsub(" ","_",herb$Taxon)
herb$match<-match(herb$tip.label,tree$tip.label)

tree_data_herb<-herb[is.na(herb$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data_herb$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data_herb,by="tip.label",sort=FALSE)


trait.values.test<-as.numeric(c(paste(trait.ordered2$value),rep(NA,389)))





herb.plot<-ggtree(pruned.tree,layout="rectangular",aes(color=trait.values.test))+ 
  geom_tippoint(,size=2, alpha=.75,show.legend=FALSE)+
  #scale_color_brewer(palette="RdYlBu")+
  scale_color_gradientn(colours=c("gray70","blue"))+
  theme(legend.position = c(0.12,0.9),legend.title = element_blank())+
  geom_tiplab(size=1)+
  geom_cladelab(node=400, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=458, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=490, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=597, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=393, label="Megaloptera", 
                offset = 0.042, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=595, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=555, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=531, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')
herb.plot

tpref.signal<-phylosig(pruned.tree,trait.ordered2$value,method="lambda",test=TRUE,)
tpref.signal ### report this value of lambda

herb_model<-gls(value~1,data=trait.ordered2,correlation = corPagel(0.5, pruned.tree))
summary(herb_model) ##
anova(herb_model)
r.squaredLR(herb_model)


tiff('./figures/herbivores.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
herb.plot
dev.off()

####################### Filterer


filter<-ffg_EU.melt[ffg_EU.melt$variable=="Filterer",]
filter$tip.label<-gsub(" ","_",filter$Taxon)
filter<-aggregate(value~tip.label,data=filter,sum)
filter$match<-match(filter$tip.label,tree$tip.label)

tree_data_filter<-filter[is.na(filter$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data_filter$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data_filter,by="tip.label",sort=FALSE)


trait.values.test<-as.numeric(c(paste(trait.ordered2$value),rep(NA,389)))




filter.plot<-ggtree(pruned.tree,layout="rectangular",aes(color=trait.values.test))+ 
  geom_tippoint(size=2, alpha=.75,show.legend=FALSE)+
  #scale_color_brewer(palette="RdYlBu")+
  scale_color_gradientn(colours=c("gray70","blue"))+
  theme(legend.position = c(0.12,0.9),legend.title = element_blank())+
  geom_tiplab(size=1)+
  geom_cladelab(node=400, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=458, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=490, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=597, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=393, label="Megaloptera", 
                offset = 0.042, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=595, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=555, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=531, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')
filter.plot



tpref.signal<-phylosig(pruned.tree,trait.ordered2$value,method="lambda",test=TRUE,)
tpref.signal ### report this value of lambda

filter_model<-gls(value~1,data=trait.ordered2,correlation = corPagel(0.5, pruned.tree))
summary(filter_model) ##
anova(filter_model)
r.squaredLR(filter_model)


tiff('./figures/filterers.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
filter.plot
dev.off()

####################### collector gathers


gather<-ffg_EU.melt[ffg_EU.melt$variable=="Gatherers..collectors",]
gather$tip.label<-gsub(" ","_",gather$Taxon)
gather$match<-match(gather$tip.label,tree$tip.label)

tree_data_gather<-gather[is.na(gather$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data_gather$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data_gather,by="tip.label",sort=FALSE)


trait.values.test<-as.numeric(c(paste(trait.ordered2$value),rep(NA,389)))




gather.plot<-ggtree(pruned.tree,layout="rectangular",aes(color=trait.values.test))+ 
  geom_tippoint(size=2, alpha=.75,show.legend=FALSE)+
  #scale_color_brewer(palette="RdYlBu")+
  scale_color_gradientn(colours=c("gray70","blue"))+
  theme(legend.position = c(0.1,0.9),legend.title = element_blank())+
  geom_tiplab(size=1)+
  geom_cladelab(node=400, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=458, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=490, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=597, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=393, label="Megaloptera", 
                offset = 0.042, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=595, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=555, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=531, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')
gather.plot

ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

ggtree(pruned.tree,layout="rectangular") + geom_text2(aes(label=node), hjust=-.3)

tiff('./figures/node_labels.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
ggtree(pruned.tree,layout="rectangular") + geom_text2(aes(label=node), hjust=-.3)
dev.off()


tpref.signal<-phylosig(pruned.tree,trait.ordered2$value,method="lambda",test=TRUE,)
tpref.signal ### report this value of lambda

gather_model<-gls(value~1,data=trait.ordered2,correlation = corPagel(0.99, pruned.tree),
                  method="ML") ## wouldn't run with REML but did with ML... not too sure why

summary(gather_model) ##
anova(gather_model)
r.squaredLR(gather_model)


tiff('./figures/gatherers.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
gather.plot
dev.off()


####################### Predator


predator<-ffg_EU.melt[ffg_EU.melt$variable=="Predators",]
predator$tip.label<-gsub(" ","_",predator$Taxon)
predator<-aggregate(value~tip.label,data=predator,sum)
predator$match<-match(predator$tip.label,tree$tip.label)

tree_data_predator<-predator[is.na(predator$match)==FALSE,]

pruned.tree<-drop.tip(tree,tree$tip.label[-match(tree_data_predator$tip.label, tree$tip.label)])

trait.ordered<-data.frame(tip.label=pruned.tree$tip.label)

trait.ordered2<-merge(trait.ordered,tree_data_predator,by="tip.label",sort=FALSE)


trait.values.test<-as.numeric(c(paste(trait.ordered2$value),rep(NA,389)))




predator.plot<-ggtree(pruned.tree,layout="rectangular",aes(color=trait.values.test))+ 
  geom_tippoint(size=2, alpha=.75,show.legend=FALSE)+
  #scale_color_brewer(palette="RdYlBu")+
  scale_color_gradientn(colours=c("gray70","blue"))+
  theme(legend.position = c(0.12,0.9),legend.title = element_blank())+
  geom_tiplab(size=1)+
  geom_cladelab(node=400, label="Ephemeroptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  
  geom_cladelab(node=458, label="Odonata", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=490, label="Plecoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=597, label="Diptera", 
                offset = 0.01, textcolor='black', barcolor='black',angle=270,offset.text=0.0015)+
  geom_cladelab(node=393, label="Megaloptera", 
                offset = 0.042, textcolor='black', barcolor='black',vjust=0.2)+
  geom_cladelab(node=595, label="Lepidoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=555, label="Trichoptera", 
                offset = 0.01, textcolor='black', barcolor='black')+
  geom_cladelab(node=531, label="Hemiptera", 
                offset = 0.01, textcolor='black', barcolor='black')
predator.plot



tpref.signal<-phylosig(pruned.tree,trait.ordered2$value,method="lambda",test=TRUE,)
tpref.signal ### report this value of lambda

predator_model<-gls(value~1,data=trait.ordered2,correlation = corPagel(0.5, pruned.tree))
summary(predator_model) ##
anova(predator_model)
r.squaredLR(predator_model)


tiff('./figures/predatorers.tiff', units="in", width=8, height=12, res=1000,compression="lzw")
predator.plot
dev.off()



