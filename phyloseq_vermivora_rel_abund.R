setwd("~/Documents/Toews_Lab/vermivora_microbiome/qiime_output/phyloseq/rel_abund/")
library(phyloseq)
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(vegan)
library(reshape2)
library(decontam)
library(gridExtra)
library(plyr)
library(lubridate)
library(maps)
library(btools) #for faiths pd
library(Maaslin2)
library(geodist)


#read in object from vermivora hz script. Already has contaminants removed, and negatives, and off target spp. is not rarefied.
#note, plumage scores are not corrected here 
physeq3<-readRDS('physeq3.rds')
table(physeq3@sam_data$Species, useNA = 'always') 

#get rid of samples not in ps3k (same ones as other dataset)
ids_for_rel_abund<-readRDS('ids_for_rel_abund.rds')
ps<-prune_samples(rownames(physeq3@sam_data) %in% ids_for_rel_abund, physeq3) 
ps<-prune_taxa(taxa_sums(ps) > 0, ps) #get rid of zero-sum OTUs
table(duplicated(ps@sam_data$Band.)) #no recaps, 123 individuals, good to go!

#add column with admixutre level
ps@sam_data$admx[ps@sam_data$Plumage_score2<0.1]<-'GWWA'
ps@sam_data$admx[ps@sam_data$Plumage_score2>0.9]<-'BWWA'
ps@sam_data$admx[(ps@sam_data$Plumage_score2<0.9 & ps@sam_data$Plumage_score2>0.1)]<-'Intermediate'
table(ps@sam_data$admx, useNA = 'always')
ps@sam_data$Year<-as.factor(ps@sam_data$Year) #make non-continuous

#calculate relative abundance of each ASV
#count/sum total ASVs
ps_ra<-transform_sample_counts(ps, function(x) x / sum(x) )

##------beta diversity-------_###
#make data.frame of metadata
md3k<-data.frame(sample_data(ps), row.names=rownames(sample_data(ps)))

jac3k<-phyloseq::distance(ps, method='jaccard', binary=T)
bray3k<-phyloseq::distance(ps, method='bray')
uni3k<-phyloseq::distance(ps, method='unifrac')
wuni3k<-phyloseq::distance(ps, method='wunifrac')

adonis2(bray3k ~ Species+State+Year, data=md3k, by='margin') #effect of year
adonis2(jac3k ~ Species+State+Year, data=md3k, by='margin') #year and state
adonis2(uni3k ~ Species+State+Year, data=md3k, by='margin') #year, species ns when accounting for state and year
adonis2(wuni3k ~ Species+State+Year, data=md3k, by='margin') #ns

plot(betadisper(bray3k, md3k$State), hull=F, ellipse=T)

#now relative abundance
md3k<-data.frame(sample_data(ps_ra), row.names=rownames(sample_data(ps_ra)))
jac3k_ra<-phyloseq::distance(ps_ra, method='jaccard', binary=T)
bray3k_ra<-phyloseq::distance(ps_ra, method='bray')
uni3k_ra<-phyloseq::distance(ps_ra, method='unifrac')
wuni3k_ra<-phyloseq::distance(ps_ra, method='wunifrac')

adonis2(bray3k_ra ~ Species+State+Year, data=md3k, by='margin') #effect of year
adonis2(jac3k_ra ~ Species+State+Year, data=md3k, by='margin') #year and state
adonis2(uni3k_ra ~ Species+State+Year, data=md3k, by='margin') #year, species ns when accounting for state and year
adonis2(wuni3k_ra ~ Species+State+Year, data=md3k, by='margin') #ns

plot(betadisper(bray3k_ra, md3k$Species), hull=F, ellipse=T)

#this plot is the same whether using ps or ps_ra
ordinate(ps, "PCoA", "unifrac") %>%
  plot_ordination(ps_ra, ., color='Species', shape = "State", title = "unifrac", type='samples') +
  scale_color_manual(values=c('violetred','skyblue','gold1')) + 
  geom_point(size=2) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(shape="Locality", colour="Species")
ordinate(ps_ra, "PCoA", "unifrac") %>%
  plot_ordination(ps_ra, ., color='Species', shape = "State", title = "unifrac-rel abund", type='samples') +
  scale_color_manual(values=c('violetred','skyblue','gold1')) + 
  geom_point(size=2) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(shape="Locality", colour="Species")

#run Maaslin2 on relative abundance object
#####----------are there ASVs that differ between parentals?------
#table S2

#model 1
parental_all<-rownames(ps_ra@sam_data[which(ps_ra@sam_data$admx=='GWWA' | ps_ra@sam_data$admx=='BWWA')]) #79 birds
ps_parental_all<-prune_samples(rownames(ps_ra@sam_data) %in% parental_all, ps_ra) 
ps_parental_all<-prune_taxa(taxa_sums(ps_parental_all) > 0, ps_parental_all) #get rid of zero-sum OTUs

maas_parental_all = Maaslin2(
  input_data = data.frame(ps_parental_all@otu_table), 
  input_metadata = data.frame(ps_parental_all@sam_data), 
  output = "maaslin2_output/parental_all", 
  fixed_effects = c("Species","State","Year"),
  reference = c('Species,BWWA')) 

#model 2
maas_ps_ra = Maaslin2(
  input_data = data.frame(ps_ra@otu_table), 
  input_metadata = data.frame(ps_ra@sam_data), 
  output = "maaslin2_output/ps_ra", 
  fixed_effects = c("Species","State","Year"),
  reference = c('Species,BWWA')) 

#model 3
maas_ps_ra_plum = Maaslin2(
  input_data = data.frame(ps_ra@otu_table), 
  input_metadata = data.frame(ps_ra@sam_data)[which(!(rownames(ps_ra@sam_data) %in% c('UE08L01','UE08L02'))),], 
  output = "maaslin2_output/ps_ra_plumage", 
  fixed_effects = c("Plumage_score2","State","Year"))

#model 4         
#subset to intermediates--table S1d
ps_int<-prune_samples(rownames(ps_ra@sam_data)[ps_ra@sam_data$admx=='Intermediate'], ps_ra) #all intermediates
ps_int<-prune_taxa(taxa_sums(ps_int) > 0, ps_int) #get rid of zero-sum OTUs
table(ps_int@sam_data$Complete_plumage, useNA='always') #all with complete plumage
maas_int_plumage = Maaslin2(
  input_data = data.frame(ps_int@otu_table), 
  input_metadata = data.frame(ps_int@sam_data), 
  output = "maaslin2_output/int_plumage", 
  fixed_effects = c("Plumage_score2", "State","Year"))

#compute carotenoid sum
carot<-ps_ra@sam_data
carot<-carot[,c('Species','admx','Nape','Back','Rump','Breast','Belly')] #pull out body carotenoids
carot$Nape<-as.numeric(carot$Nape)
carot$Back<-as.numeric(carot$Back)
carot$Rump<-as.numeric(carot$Rump)
carot$Breast<-as.numeric(carot$Breast)
carot$Belly<-as.numeric(carot$Belly)
carot$sum<-rowSums(carot[,3:7], na.rm = T)
carot[,3:8]#need to turn sum into NA for birds with missing scores
carot$sum[!complete.cases(carot[,3:7])] <-NA #done (birds w/missing carot plumage have NA)
#add wing bar color, but conver to low score=yellow to match scale for body carotenoids
identical(rownames(carot), rownames(ps_ra@sam_data))
carot$wbar<-5-as.numeric(ps_ra@sam_data$W_bar) #max score is 5, so 5-w_bar reverses the scale
carot$sum2<-rowSums(carot[,c(3:7,9)], na.rm = T) #with wing bar
carot[,c(3:8,9:10)]#need to turn sum into NA for birds with missing scores
carot$sum2[!complete.cases(carot[,c(3:7,9)])] <-NA #done (birds w/missing carot plumage have NA)

ps_ra@sam_data$carot_sum2<-carot$sum2[match(rownames(ps_ra@sam_data), rownames(carot))] #carot sum corrected (none affected anyways)

#model 5
maas_int_carot = Maaslin2(
  input_data = data.frame(ps_int@otu_table), 
  input_metadata = data.frame(ps_int@sam_data), 
  output = "maaslin2_output/int_carot", 
  fixed_effects = c("carot_sum2","State","Year"))

#model 6
maas_ps_ra_carot = Maaslin2(
  input_data = data.frame(ps_ra@otu_table), 
  input_metadata = data.frame(ps_ra@sam_data)[which(!(rownames(ps_ra@sam_data) %in% c('UE08L01','UE08L02'))),], 
  output = "maaslin2_output/ps_ra_carot", 
  fixed_effects = c("carot_sum2","State","Year"))
