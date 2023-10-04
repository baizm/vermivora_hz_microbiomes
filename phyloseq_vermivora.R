setwd("~/Documents/Toews_Lab/vermivora_microbiome/qiime_output/phyloseq/")
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

#read in post-decontam filtered data from qiime2
########-------create physeq object---------------------########
####------read in OTU table, has to be a matrix---------------------
otu_table<-read.csv('../final/feature-table.tsv', sep='\t', header=T, skip=1, check.names = F)
#convert to matrix
otumat<-as.matrix(otu_table[,2:ncol(otu_table)])
#add rownames that are OTU id
rownames(otumat)<-otu_table[,1]
###-------read in taxonomy table-----------------
tax_table<-read.csv('../decontam/silva_taxonomy.tsv', sep='\t', header=F, skip=2)
tax_table2<-separate(data = tax_table, col = V2, 
                     into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = "; ")
colnames(tax_table2)[1]<-'Feature.ID'
colnames(tax_table2)[9]<-'Confidence'
#it wasn't filtered, so i'll have to delete mt,cp,un,euk,decontam...
tax_table2<-tax_table2[which(tax_table2$Feature.ID %in% otu_table$`#OTU ID`),] #9239 OTUs after filtering all mt,cp,un,euk,decontam
#convert to matrix
taxmat<-as.matrix(tax_table2[,2:8])
rownames(taxmat)<-tax_table2$Feature.ID
###-------read in sample data---------------------
sampledata<-read.csv('../../metadata/metadata_16s_vermivora2.tsv', sep='\t', header=T,)
rownames(sampledata)<-sampledata$id
x<-sampledata$id[(which(!(sampledata$id %in% colnames(otumat))))]
sampledata<-sampledata[which(!(sampledata$id %in% x)),] #get rid of 3 samples filtered by qiime (2 neg + 1 low count)
sampledata<-sampledata[,2:ncol(sampledata)]
###-------read in nwk tree with ape package-----------
tree<-read.tree('../final/tree.nwk')
###----------combine into phyloseq object-------------
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sampledata)
physeq = phyloseq(OTU, TAX, SAM, tree)
#sums of OTUs across samples
taxa_sums(physeq)
physeq #9239 taxa and 190 samples

###---housekeeping----###
table(physeq@sam_data$Species, useNA='always')
#remove EABL, PRWA, and negatives
sp_to_remove<-c('EABL','PRWA') #9 Negs/NA, 1 PRWA, 22 EABL=32
#remove these samples
physeq2<-prune_samples(!(grepl(paste(sp_to_remove,collapse='|'),sample_data(physeq)$Species)), physeq)
#pull out all but negatives
exl_negs<-rownames(physeq2@sam_data)[which(!(grepl('neg', sample_data(physeq2)$Type)))]
physeq3<-prune_samples(rownames(physeq2@sam_data) %in% exl_negs, physeq2) 
physeq3<-prune_taxa(taxa_sums(physeq3) > 0, physeq3) #get rid of zero-sum OTUs
table(physeq3@sam_data$Species, useNA='always') 
physeq3 #5923 taxa and 158 samples
length(unique(physeq3@sam_data$Band.)) #144 individuals
#write physeq3 to redo analyses using relative abundance of OTUs
saveRDS(physeq3, "/Users/marcellabaiz/Documents/Toews_Lab/vermivora_microbiome/qiime_output/phyloseq/rel_abund/physeq3.rds")
#id recaps
table(duplicated(physeq3@sam_data$Band.)) #14 duplicated
recap_bands<-physeq3@sam_data$Band.[which(duplicated(physeq3@sam_data$Band.))]
recap_data_all<-physeq3@sam_data[which(physeq3@sam_data$Band. %in% recap_bands),]
recap_data_all<-recap_data_all[order(recap_data_all$Band.),]

####---------find rarefaction depth-------------
#make data frame with total read counts 
sdt = data.table(as(sample_data(physeq3), "data.frame"),
                 TotalReads = sample_sums(physeq3), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
summary(sdt$TotalReads) #range 39-88461
#look at rarecurve
tab <- otu_table(physeq3)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
rarecurve(tab, step=50, lwd=0.8, label=F, xlab='Read count',
          ylab='ASVs',xlim=c(0,7000))

#saveRDS(sdt, file='sdt.RDS')

ps3k<-rarefy_even_depth(physeq3,sample.size=3076,rngseed = 999, replace = F, trimOTUs = T, verbose = T)

#how many females?
table(ps3k@sam_data$Sex, useNA='always') # 7 females #the blank one is probably male (brenda marked cp=3 on datasheet)

#------remove recaps from dataset-----####
recap_bands2<-ps3k@sam_data$Band.[which(ps3k@sam_data$Band. %in% recap_bands)]
recap_data_all2<-ps3k@sam_data[which(ps3k@sam_data$Band. %in% recap_bands2),]
recap_data_all2<-recap_data_all2[order(recap_data_all2$Band.),]
recaps_to_keep<-c('TE27T05','UF01P01','TF40T04','UF03T03','UE28T01','TE27T02','SF11T02', #some are not recaps after filtering
                  'UF10B01','VF28B02','VF22B04','UG02B01','UF10K01','VF02B01') #selected to balance sample numbers across years, etc.
recaps_to_exclude<-rownames(recap_data_all2)[!(rownames(recap_data_all2) %in% recaps_to_keep)]
length(recaps_to_exclude) #10 recap samples
#write vector for rel_abund script
saveRDS(recaps_to_exclude, "/Users/marcellabaiz/Documents/Toews_Lab/vermivora_microbiome/qiime_output/phyloseq/rel_abund/recaps_to_exclude.rds")

ps3k<-prune_samples(!(rownames(ps3k@sam_data) %in% recaps_to_exclude), ps3k) 
ps3k<-prune_taxa(taxa_sums(ps3k) > 0, ps3k) #get rid of zero-sum OTUs
table(duplicated(ps3k@sam_data$Band.))

#write vector for rel_abund script
ids_for_rel_abund<-rownames(ps3k@sam_data)
saveRDS(ids_for_rel_abund, "/Users/marcellabaiz/Documents/Toews_Lab/vermivora_microbiome/qiime_output/phyloseq/rel_abund/ids_for_rel_abund.rds")

ps3k # 4564 taxa and 123 samples

#----check out plumage scores-----#
hist(ps3k@sam_data$Plumage_score2, breaks=30)
table(ps3k@sam_data$Complete_plumage, useNA = 'always') #4 incomplete plumage
hist(ps3k@sam_data$Plumage_score2[which(ps3k@sam_data$Complete_plumage=='Y')], breaks=30) #very similar

#for figure
h = hist(ps3k@sam_data$Plumage_score2, breaks=30,plot=FALSE)
ccat = cut(h$breaks, c(-1, 0.09, 0.89, 1))
pdf('plots2/hist_plumage.pdf', height=3, width=3.5)
plot(h, col=c('gold','gray','cornflowerblue')[ccat], main='', xlab='') 
dev.off()

#add column with admixutre level
ps3k@sam_data$admx[ps3k@sam_data$Plumage_score2<0.1]<-'GWWA'
ps3k@sam_data$admx[ps3k@sam_data$Plumage_score2>0.9]<-'BWWA'
ps3k@sam_data$admx[(ps3k@sam_data$Plumage_score2<0.9 & ps3k@sam_data$Plumage_score2>0.1)]<-'Intermediate'
table(ps3k@sam_data$admx, useNA = 'always') #32 admixed
admixed_ids<-ps3k@sam_data$Band.[which(ps3k@sam_data$admx=='Intermediate')]
saveRDS(admixed_ids, file='admixed_ids.RDS')

###----map loaclities----###
map(database='state',regions=c("Michigan","Ohio","Pennsylvania"),col="black", lwd=.8, resolution=.1)
points(ps3k@sam_data$Longitude, ps3k@sam_data$Latitude, pch=21, bg=24, lwd=.4)

##------beta diversity-------_###
#remove birds without plumage scores (admx=NA) so we can test admixture as a variable
table(is.na(ps3k@sam_data$admx))
keep_nonna1<-rownames(ps3k@sam_data)[which(!(is.na(ps3k@sam_data$admx)))]
ps3k_nar<-prune_samples(rownames(ps3k@sam_data) %in% keep_nonna1, ps3k) 
ps3k_nar<-prune_taxa(taxa_sums(ps3k_nar) > 0, ps3k_nar) #get rid of zero-sum OTUs

#make data.frame of metadata
md3k_nar<-data.frame(sample_data(ps3k_nar), row.names=rownames(sample_data(ps3k_nar)))

jac3k_nar<-phyloseq::distance(ps3k_nar, method='jaccard', binary=T)
bray3k_nar<-phyloseq::distance(ps3k_nar, method='bray')
uni3k_nar<-phyloseq::distance(ps3k_nar, method='unifrac')
wuni3k_nar<-phyloseq::distance(ps3k_nar, method='wunifrac')

#no effect of admx still-strong year effect
adonis2(bray3k_nar ~ admx+State+Year, data=md3k_nar, by='margin') #year
permutest(betadisper(bray3k_nar, md3k_nar$Year)) #F=20.499,p=0.001***
adonis2(jac3k_nar ~ admx+State+Year, data=md3k_nar, by='margin') #state, year
permutest(betadisper(jac3k_nar, md3k_nar$State)) #F=20.653,p=0.001***
permutest(betadisper(jac3k_nar, md3k_nar$Year)) #F=18.373,p=0.001***
adonis2(uni3k_nar ~ admx+State+Year, data=md3k_nar, by='margin') #year
permutest(betadisper(uni3k_nar, md3k_nar$Year)) #F=2.2887,p=0.073
adonis2(wuni3k_nar ~ Species+State+Year, data=md3k_nar, by='margin') #ns

#figs for supplement
pdf(file='plots2/supp/betadisper_bray3k_year.pdf', height=5, width=5)
plot(betadisper(bray3k_nar, md3k_nar$Year), seg.lwd = 0.5, main='')
dev.off()
pdf(file='plots2/supp/betadisper_jac3k_state.pdf', height=5, width=5)
plot(betadisper(jac3k_nar, md3k_nar$State), seg.lwd = 0.5, main='')
dev.off()
pdf(file='plots2/supp/betadisper_jac3k_year.pdf', height=5, width=5)
plot(betadisper(jac3k_nar, md3k_nar$Year), seg.lwd = 0.5, main='')
dev.off()

#calculate alpha values
alpha_r<-estimate_richness(ps3k, measures=c('Shannon','Observed','Chao1'))
alpha_r$PD<-estimate_pd(ps3k)[,1] #grab first column that has faith's pd (second column is species richness--#otus)
alpha_r$Species<-ps3k@sam_data$Species
alpha_r$Plumage_score2<-ps3k@sam_data$Plumage_score2
alpha_r$Complete_plumage<-ps3k@sam_data$Complete_plumage
alpha_r$admx<-ps3k@sam_data$admx
alpha_r$State<-ps3k@sam_data$State
alpha_r$Year<-ps3k@sam_data$Year

##--run linear model on alpha to test multiple variables
alpha_nar<-alpha_r[which(!(is.na(alpha_r$admx))),] #make dataframe of birds w/plumage scores, n=111
#build basic model
lm_chao1 <- lm(Chao1 ~ admx + State + Year, data=alpha_nar)
#run a marginal test Anova (Type II)
lm_chao1 %>% car::Anova() #optional, add test.statistic="F" for F-statistic

lm_shan <- lm(Shannon ~ admx + State + Year, data=alpha_nar)
#run a marginal test Anova (Type II)
lm_shan %>% car::Anova() 

lm_pd <- lm(PD ~ admx + State + Year, data=alpha_nar)
#run a marginal test Anova (Type II)
lm_pd %>% car::Anova() 

#levels for plot
alpha_r$admx<-factor(alpha_r$admx, levels=c('GWWA','Intermediate','BWWA')) #order for plot axis
boxplot(alpha_r$Shannon~alpha_r$admx,xlab='', ylab='Microbiome alpha diversity (Shannon index)',
        main='', boxwex=0.5, outline=F, ylim=c(0,5.5), col=c('gold','gray','cornflowerblue'),
        names=c('GW',"Intermediate",'BW'))
stripchart(alpha_r$Shannon[which(alpha_r$admx=='GWWA')],vertical=T,pch=21,bg='gold', method='jitter',add=T)
stripchart(alpha_r$Shannon[which(alpha_r$admx=='Intermediate')],vertical=T,pch=21,bg='gray', method='jitter',add=T,at=2)
stripchart(alpha_r$Shannon[which(alpha_r$admx=='BWWA')],vertical=T,pch=21,bg='cornflowerblue', method='jitter',add=T,at=3)

#alpha across years
boxplot(alpha_r$Shannon~alpha_r$Year,xlab='', ylab='Microbiome alpha diversity (Shannon index)',
        main='', boxwex=0.5, outline=F, ylim=c(0,5.6))


###-----are there differences between parental populations?-----###
#---pull out birds with extreme plumage scores
parental_all<-rownames(ps3k@sam_data[which(ps3k@sam_data$admx=='GWWA' | ps3k@sam_data$admx=='BWWA')]) #79 birds
ps_parental_all<-prune_samples(rownames(ps3k@sam_data) %in% parental_all, ps3k) 
ps_parental_all<-prune_taxa(taxa_sums(ps_parental_all) > 0, ps_parental_all) #get rid of zero-sum OTUs
table(ps_parental_all@sam_data$Species) 

plot_richness(ps_parental_all, x='Species', color='Species',measures=c('Observed','Chao1','Shannon')) + 
  geom_boxplot() + theme(legend.position="none") 
alpha_parental_all<-alpha_r[which(rownames(alpha_r) %in% rownames(ps_parental_all@sam_data)),]


#beta diversity for parentals
jac_par<-phyloseq::distance(ps_parental_all, method='jaccard', binary=T)
bray_par<-phyloseq::distance(ps_parental_all, method='bray')
uni_par<-phyloseq::distance(ps_parental_all, method='unifrac')
wuni_par<-phyloseq::distance(ps_parental_all, method='wunifrac')

md_par<-data.frame(sample_data(ps_parental_all), row.names=rownames(sample_data(ps_parental_all)))
table(ps_parental_all@sam_data$admx, ps_parental_all@sam_data$admx, useNA='always') #species is the same as admx
adonis2(bray_par ~ Species+State+Year, data=md_par, by='margin') #species ns, year sig
permutest(betadisper(bray_par, md_par$Year)) #F=7.7247,p=0.001***
adonis2(jac_par ~ Species+State+Year, data=md_par, by='margin') #species ns, state and year sig
permutest(betadisper(jac_par, md_par$State)) #F=10.84,p=0.002***
permutest(betadisper(jac_par, md_par$Year)) #F=24.006,p=0.001***
adonis2(uni_par ~ Species+State+Year, data=md_par, by='margin') #species ns, year sig
permutest(betadisper(uni_par, md_par$Year)) #F=4.3919,p=0.006***
adonis2(wuni_par ~ Species+State+Year, data=md_par, by='margin') #species ns

plot(betadisper(bray_par, md_par$Year), seg.lwd = 0.5, main='')
plot(betadisper(jac_par, md_par$State), seg.lwd = 0.5, main='')
plot(betadisper(jac_par, md_par$Year), seg.lwd = 0.5, main='')
plot(betadisper(uni_par, md_par$Year), seg.lwd = 0.5, main='')

#does alpha diversity differ between parental individuals?
#build basic linear model
table(alpha_parental_all$admx, useNA = 'always')
lm_chao1_par <- lm(Chao1 ~ admx + State + Year, data=alpha_parental_all)
#run a marginal test Anova (Type II)
lm_chao1_par %>% car::Anova() #optional, add test.statistic="F" for F-statistic

lm_shan_par <- lm(Shannon ~ admx + State + Year, data=alpha_parental_all)
#run a marginal test Anova (Type II)
lm_shan_par %>% car::Anova() 

lm_pd_par <- lm(PD ~ admx + State + Year, data=alpha_parental_all)
#run a marginal test Anova (Type II)
lm_pd_par %>% car::Anova() 


####-----separate PA and MI-----------
pa_all<-rownames(ps3k@sam_data[which(ps3k@sam_data$State=='PA')]) #73 birds
ps_pa_all<-prune_samples(rownames(ps3k@sam_data) %in% pa_all, ps3k) 
ps_pa_all<-prune_taxa(taxa_sums(ps_pa_all) > 0, ps_pa_all) #get rid of zero-sum OTUs
table(ps_pa_all@sam_data$admx, useNA = 'always') 
ps_pa_all@sam_data[which(is.na(ps_pa_all@sam_data$admx))] #one of the sproul birds without photos
table(ps_pa_all@sam_data$Region)

#--PA parentals just 2021-2022 
pa_2021_2022<-rownames(ps_pa_all@sam_data)[which(ps_pa_all@sam_data$Year=='2021' | ps_pa_all@sam_data$Year=='2022')]
ps_pa_parental_2021_2022<-prune_samples(rownames(ps_pa_all@sam_data) %in% pa_2021_2022, ps_pa_all) #all PA 2021-2022 birds
ps_pa_parental_2021_2022<-prune_samples(rownames(ps_pa_parental_2021_2022@sam_data) %in% parental_all, ps_pa_parental_2021_2022) #all parental PA 2021-2022 birds
ps_pa_parental_2021_2022<-prune_taxa(taxa_sums(ps_pa_parental_2021_2022) > 0, ps_pa_parental_2021_2022) #get rid of zero-sum OTUs
table(ps_pa_parental_2021_2022@sam_data$admx, ps_pa_parental_2021_2022@sam_data$Year)

pdf('plots2/pcoa_pa_parentals_unifrac.pdf',height=3.5, width=4)
ordinate(ps_pa_parental_2021_2022, "PCoA", "unifrac") %>%
  plot_ordination(ps_pa_parental_2021_2022, ., color='Species', title = "Appalachian HZ (UniFrac)", type='samples') +
  scale_color_manual(values=c('cornflowerblue','gold')) + 
  geom_point(size=2) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(colour="Species", shape='Region') +
  stat_ellipse(aes(color=Species, group=Species), level=0.5, lty=2, lwd=0.4)
dev.off()

pdf('plots2/pcoa_pa_parentals_unifrac2.pdf',height=3.5, width=4)
ordinate(ps_pa_parental_2021_2022, "PCoA", "unifrac") %>%
  plot_ordination(ps_pa_parental_2021_2022, ., color='Species', title = "Appalachian HZ (UniFrac)", type='samples') +
  scale_color_manual(values=c('black','black')) +
  scale_fill_manual(values=c('cornflowerblue','gold')) + 
  scale_shape_manual(values=c(21,21)) +
  geom_point(aes(fill = Species, shape = Species), size=3) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + labs(colour="Species") +
  stat_ellipse(aes(fill=Species), geom='polygon', level=0.5, lty=2, lwd=0.4, alpha=0.15) 
#geom_convexhull(aes(fill=Species),alpha=0.3)
dev.off()

#beta diversity for PA parentals
jac_pa_par2122<-phyloseq::distance(ps_pa_parental_2021_2022, method='jaccard', binary=T)
bray_pa_par2122<-phyloseq::distance(ps_pa_parental_2021_2022, method='bray')
uni_pa_par2122<-phyloseq::distance(ps_pa_parental_2021_2022, method='unifrac')
wuni_pa_par2122<-phyloseq::distance(ps_pa_parental_2021_2022, method='wunifrac')

md_pa_par2122<-data.frame(sample_data(ps_pa_parental_2021_2022), row.names=rownames(sample_data(ps_pa_parental_2021_2022)))
adonis2(bray_pa_par2122 ~ Species+Region+Year, data=md_pa_par2122, by='margin') # ns
adonis2(jac_pa_par2122 ~ Species+Region+Year, data=md_pa_par2122, by='margin') # species marginal, region, year
adonis2(uni_pa_par2122 ~ Species+Region+Year, data=md_pa_par2122, by='margin') #ns
adonis2(wuni_pa_par2122 ~ Species+Region+Year, data=md_pa_par2122, by='margin') #ns

permutest(betadisper(jac_pa_par2122, md_pa_par2122$Year)) #F=5.9142,p=0.023**
plot(betadisper(jac_pa_par2122, md_pa_par2122$Year), seg.lwd = 0.5, main='')
permutest(betadisper(jac_pa_par2122, md_pa_par2122$Region)) #F=0.4188,p=0.672
plot(betadisper(jac_pa_par2122, md_pa_par2122$Region), seg.lwd = 0.5, main='')

#does alpha diversity differ between parental individuals--PA only (21-22)?
#build basic linear model
alpha_pa_parental_2122<-alpha_parental_all[which(alpha_parental_all$State=='PA' & alpha_parental_all$Year %in% c('2021','2022')),]
alpha_pa_parental_2122$Region<-ps3k@sam_data$Region[match(rownames(alpha_pa_parental_2122),rownames(ps3k@sam_data))] #add region column
table(alpha_pa_parental_2122$Species, alpha_pa_parental_2122$Year, useNA = 'always')
lm_chao1_par_pa <- lm(Chao1 ~ admx + Region + Year, data=alpha_pa_parental_2122)
#run a marginal test Anova (Type II)
lm_chao1_par_pa %>% car::Anova() #ns

lm_shan_par_pa <- lm(Shannon ~ admx + Region + Year, data=alpha_pa_parental_2122)
#run a marginal test Anova (Type II)
lm_shan_par_pa %>% car::Anova() #ns

lm_pd_par_pa <- lm(PD ~ admx + Region + Year, data=alpha_pa_parental_2122)
#run a marginal test Anova (Type II)
lm_pd_par_pa %>% car::Anova() #ns


###-------michigan--------##
mi_all<-rownames(ps3k@sam_data[which(ps3k@sam_data$State=='MI')]) #50 birds
ps_mi_all<-prune_samples(rownames(ps3k@sam_data) %in% mi_all, ps3k) 
ps_mi_all<-prune_taxa(taxa_sums(ps_mi_all) > 0, ps_mi_all) #get rid of zero-sum OTUs
table(ps_mi_all@sam_data$admx, useNA = 'always') 
ps_mi_all@sam_data[which(is.na(ps_mi_all@sam_data$admx))] #all kzoo birds

table(ps_mi_all@sam_data$Region)

#add lat long to alpha table
alpha_r$Latitude<-ps3k@sam_data$Latitude
alpha_r$Longitude<-ps3k@sam_data$Longitude
alpha_mi<-alpha_r[which(alpha_r$State=='MI'),]

plot(alpha_mi$Latitude~alpha_mi$Shannon, col=factor(alpha_mi$Species))
plot(alpha_mi$Latitude~alpha_mi$PD)
plot(alpha_mi$Latitude~alpha_mi$Chao1)

#michigan parentals--2021, 2022
table(ps_mi_all@sam_data$admx)
ps_mi_parental<-prune_samples(rownames(ps_mi_all@sam_data) %in% parental_all, ps_mi_all) #all mi parental birds
ps_mi_parental<-prune_taxa(taxa_sums(ps_mi_parental) > 0, ps_mi_parental) #get rid of zero-sum OTUs
table(ps_mi_parental@sam_data$admx, ps_mi_parental@sam_data$Region)
table(ps_mi_parental@sam_data$admx, useNA = 'always')

pdf('plots2/pcoa_mi_parentals_unifrac2.pdf',height=3.5, width=4)
ordinate(ps_mi_parental, "PCoA", "unifrac") %>%
  plot_ordination(ps_mi_parental, ., color='Species', title = "Great Lakes HZ (Unifrac)", type='samples') +
  scale_color_manual(values=c('black','black')) +
  scale_fill_manual(values=c('cornflowerblue','gold')) + 
  scale_shape_manual(values=c(21,21)) +
  geom_point(aes(fill = Species, shape = Species), size=3) + theme(axis.text = element_text(size = 14)) +
  theme_bw() + labs(colour="Species") +
  stat_ellipse(aes(fill=Species), geom='polygon', level=0.5, lty=2, lwd=0.4, alpha=0.15) 
#geom_convexhull(aes(fill=Species),alpha=0.3)
dev.off()

#beta diversity for MI parental birds
jac_mi_par<-phyloseq::distance(ps_mi_parental, method='jaccard', binary=T)
bray_mi_par<-phyloseq::distance(ps_mi_parental, method='bray')
uni_mi_par<-phyloseq::distance(ps_mi_parental, method='unifrac')
wuni_mi_par<-phyloseq::distance(ps_mi_parental, method='wunifrac')

md_mi_par<-data.frame(sample_data(ps_mi_parental), row.names=rownames(sample_data(ps_mi_parental)))
adonis2(bray_mi_par ~ Species+Region+Year, data=md_mi_par, by='margin') # ns
adonis2(jac_mi_par ~ Species+Region+Year, data=md_mi_par, by='margin') # ns
adonis2(uni_mi_par ~ Species+Region+Year, data=md_mi_par, by='margin') # species r2=0.05446, 0.046, region r2=0.10604, p=0.033
adonis2(wuni_mi_par ~ Species+Region+Year, data=md_mi_par, by='margin') #ns

permutest(betadisper(uni_mi_par, md_mi_par$Species)) #F=3.7216,p=0.058
plot(betadisper(uni_mi_par, md_mi_par$Species), seg.lwd = 0.5, main='')
permutest(betadisper(uni_mi_par, md_mi_par$Region)) #F=0.1856,p=0.827
plot(betadisper(uni_mi_par, md_mi_par$Region), seg.lwd = 0.5, main='')

#--does alpha diversity differ between parental individuals--MI only (2021-22)?
#build basic linear model
alpha_mi_parental2<-alpha_parental_all[which(alpha_parental_all$State=='MI'),]
alpha_mi_parental2$Region<-ps3k@sam_data$Region[match(rownames(alpha_mi_parental2),rownames(ps3k@sam_data))] #add region column
table(alpha_mi_parental2$admx, useNA = 'always')
lm_chao1_par_mi <- lm(Chao1 ~ admx + Region + Year, data=alpha_mi_parental2)
#run a marginal test Anova (Type II)
lm_chao1_par_mi %>% car::Anova() #region p=0.03982, F=3.7453

plot_richness(ps_mi_parental, x="Region", color="Species",measures=c("Chao1","Shannon")) #more ASVs in the north

boxplot(alpha_mi_parental2$Chao1~ alpha_mi_parental2$Region + alpha_mi_parental2$Species, xlab='', ylab="Alpha diversity (Chao1 index)",
        at = c(2,1,3,5,4,6), boxwex=0.5, outline=F, col=c('steelblue2','steelblue4','steelblue3','gold','gold3','gold4'))
stripchart(alpha_mi_parental2$Chao1[which(alpha_mi_parental2$Species=='BWWA' & alpha_mi_parental2$Region=='LPS')],vertical=T,pch=21,bg='steelblue4',method='jitter',add=T)
stripchart(alpha_mi_parental2$Chao1[which(alpha_mi_parental2$Species=='BWWA' & alpha_mi_parental2$Region=='LPN')],vertical=T,pch=21,bg='steelblue2',method='jitter',add=T,at=2)
stripchart(alpha_mi_parental2$Chao1[which(alpha_mi_parental2$Species=='BWWA' & alpha_mi_parental2$Region=='UP')],vertical=T,pch=21,bg='steelblue3',method='jitter',add=T,at=3)
stripchart(alpha_mi_parental2$Chao1[which(alpha_mi_parental2$Species=='GWWA' & alpha_mi_parental2$Region=='LPN')],vertical=T,pch=21,bg='gold',method='jitter',add=T,at=5)
stripchart(alpha_mi_parental2$Chao1[which(alpha_mi_parental2$Species=='GWWA' & alpha_mi_parental2$Region=='UP')],vertical=T,pch=21,bg='gold4',method='jitter',add=T,at=6)

alpha_mi_parental2$Latitude<-ps3k@sam_data$Latitude[match(rownames(alpha_mi_parental2),rownames(ps3k@sam_data))] #add lat column
plot(alpha_mi_parental2$Chao1~alpha_mi_parental2$Latitude,pch=21,bg=c('cornflowerblue','gold')[factor(alpha_mi_parental2$Species)],
     xlab='Latitude',ylab='Chao1 index', main='Great Lakes HZ',cex=1.5)
lm(Chao1 ~ Latitude + admx + Year, data=alpha_mi_parental2) %>% car::Anova()  #F=5.4880, P=0.02817

pdf('plots2/cor_lat_chao1_greatlakes2.pdf', height=4, width=5)
ggplot(alpha_mi_parental2, aes(x = Latitude, y = Chao1)) + 
  geom_point(aes(fill = Species, shape = Species), size=3) + 
  scale_color_manual(values=c('black','black')) + scale_fill_manual(values=c('cornflowerblue','gold')) +
  scale_shape_manual(values=c(21,21)) +
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2) +
  coord_cartesian(ylim = c(-0.2,max(alpha_mi_parental2$Chao1)))
dev.off()

lm_shan_par_mi <- lm(Shannon ~ admx + Region + Year, data=alpha_mi_parental2)
#run a marginal test Anova (Type II)
lm_shan_par_mi %>% car::Anova() 

lm_pd_par_mi <- lm(PD ~ admx + Region + Year, data=alpha_mi_parental2)
#run a marginal test Anova (Type II)
lm_pd_par_mi %>% car::Anova() 


##-------just parental BWWA--------
ps_bwwa_all<-prune_samples(rownames(ps_parental_all@sam_data)[ps_parental_all@sam_data$Species=='BWWA'], ps_parental_all) 
ps_bwwa_all<-prune_taxa(taxa_sums(ps_bwwa_all) > 0, ps_bwwa_all) #get rid of zero-sum OTUs
table(ps_bwwa_all@sam_data$State, useNA='always')

#bwwa parentals in 2021 and 2022 only 
keep_bwwa_2122<-rownames(ps_bwwa_all@sam_data)[which(ps_bwwa_all@sam_data$Year=='2021' | ps_bwwa_all@sam_data$Year=='2022')]
ps_bwwa_2122<-prune_samples(rownames(ps_bwwa_all@sam_data) %in% keep_bwwa_2122, ps_bwwa_all) 
ps_bwwa_2122<-prune_taxa(taxa_sums(ps_bwwa_2122) > 0, ps_bwwa_2122) #get rid of zero-sum OTUs
table(ps_bwwa_2122@sam_data$State, useNA='always')

pdf('plots2/pcoa_bwwa_par_state2122_jacc.pdf', height=3, width=3.5)
ordinate(ps_bwwa_2122, "PCoA", "jaccard") %>%
  plot_ordination(ps_bwwa_2122, ., shape='State', title = "BWWA parentals 2021-22, Jaccard", type='samples') +
  scale_fill_manual(values=c('cornflowerblue')) +scale_shape_manual(values=c(21, 24)) +
  geom_point(size=3,  aes(fill='cornflowerblue'), show.legend = F) + 
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(colour="State") +
  stat_ellipse(geom = 'path' , aes(color = State, group=State,linetype=factor(ps_bwwa_2122@sam_data$State)), level=0.5, lwd=0.6) +
  scale_linetype_manual(values=c(2,2), guide='none')
dev.off()

#beta diversity for bw 2122 birds
jac_bwwa_2122<-phyloseq::distance(ps_bwwa_2122, method='jaccard', binary=T)
bray_bwwa_2122<-phyloseq::distance(ps_bwwa_2122, method='bray')
uni_bwwa_2122<-phyloseq::distance(ps_bwwa_2122, method='unifrac')
wuni_bwwa_2122<-phyloseq::distance(ps_bwwa_2122, method='wunifrac')

md_bw_2122<-data.frame(sample_data(ps_bwwa_2122), row.names=rownames(sample_data(ps_bwwa_2122)))
adonis2(bray_bwwa_2122 ~ State+Year, data=md_bw_2122, by='margin') # year r2=0.03465, p=0.038
adonis2(jac_bwwa_2122 ~ State+Year, data=md_bw_2122, by='margin') #year r2=0.02862, p=0.006, state=0.02634 p=0.041
adonis2(uni_bwwa_2122 ~ State+Year, data=md_bw_2122, by='margin') #ns
adonis2(wuni_bwwa_2122 ~ State+Year, data=md_bw_2122, by='margin') #ns

permutest(betadisper(bray_bwwa_2122, md_bw_2122$Year)) #F=0.451,p=0.492
permutest(betadisper(jac_bwwa_2122, md_bw_2122$Year)) #F=1.2464,p=0.262
permutest(betadisper(jac_bwwa_2122, md_bw_2122$State)) #F=4.2242,p=0.05

#alpha bw parental 21-22
#build basic linear model
alpha_bw2122<-alpha_r[which(rownames(alpha_r) %in% rownames(ps_bwwa_2122@sam_data)),]
lm_chao1_bw2122 <- lm(Chao1 ~ State + Year, data=alpha_bw2122)
#run a marginal test Anova (Type II)
lm_chao1_bw2122 %>% car::Anova() #ns
plot_richness(ps_bwwa_2122, x="State", color="Year",measures=c("Chao1","Shannon")) #may be effect of UP-3 of 4 samples collected in 2022

lm_shan_bw2122 <- lm(Shannon ~ State + Year, data=alpha_bw2122)
#run a marginal test Anova (Type II)
lm_shan_bw2122 %>% car::Anova() #ns

lm_pd_bw2122 <- lm(PD ~ State + Year, data=alpha_bw2122)
#run a marginal test Anova (Type II)
lm_pd_bw2122 %>% car::Anova() #ns

##-------just parental GWWA--------
ps_gwwa_all<-prune_samples(rownames(ps_parental_all@sam_data)[ps_parental_all@sam_data$Species=='GWWA'], ps_parental_all) 
ps_gwwa_all<-prune_taxa(taxa_sums(ps_gwwa_all) > 0, ps_gwwa_all) #get rid of zero-sum OTUs
table(ps_gwwa_all@sam_data$State, ps_gwwa_all@sam_data$Year, useNA='always')

#gwwa parentals 2021-2022 only
ps_gwwa_2122<-prune_samples(rownames(ps_gwwa_all@sam_data)[which(ps_gwwa_all@sam_data$Year %in% c('2021','2022'))] , ps_gwwa_all) 
ps_gwwa_2122<-prune_taxa(taxa_sums(ps_gwwa_2122) > 0, ps_gwwa_2122) #get rid of zero-sum OTUs
table(ps_gwwa_2122@sam_data$State, ps_gwwa_2122@sam_data$Year, useNA='always')

#beta diversity for gwwa 21-22 parental birds
jac_gwwa_2122<-phyloseq::distance(ps_gwwa_2122, method='jaccard', binary=T)
bray_gwwa_2122<-phyloseq::distance(ps_gwwa_2122, method='bray')
uni_gwwa_2122<-phyloseq::distance(ps_gwwa_2122, method='unifrac')
wuni_gwwa_2122<-phyloseq::distance(ps_gwwa_2122, method='wunifrac')

md_gw_2122<-data.frame(sample_data(ps_gwwa_2122), row.names=rownames(sample_data(ps_gwwa_2122)))
adonis2(bray_gwwa_2122 ~ State+Year, data=md_gw_2122, by='margin') # ns
adonis2(jac_gwwa_2122 ~ State+Year, data=md_gw_2122, by='margin') # state r2=0.06207, p=0.002, year r2=0.05939, p=0.001
adonis2(uni_gwwa_2122 ~ State+Year, data=md_gw_2122, by='margin') #ns
adonis2(wuni_gwwa_2122 ~ State+Year, data=md_gw_2122, by='margin') #ns

permutest(betadisper(bray_gwwa_2122, md_gw_2122$State)) #F=0.7661,p=0.359
plot(betadisper(bray_gwwa_2122, md_gw_2122$State), seg.lwd = 0.5, main='')
permutest(betadisper(bray_gwwa_2122, md_gw_2122$Year)) #F=0.0571,p=0.831
plot(betadisper(bray_gwwa_2122, md_gw_2122$Year), seg.lwd = 0.5, main='')

pdf('plots2/pcoa_gwwa_par_state_jacc2122.pdf', height=3, width=3.5)
ordinate(ps_gwwa_2122, "PCoA", "jaccard") %>%
  plot_ordination(ps_gwwa_2122, ., shape='State', title = "GWWA parentals 2021-22, Jaccard", type='samples') +
  scale_fill_manual(values=c('gold','gold')) +  scale_shape_manual(values=c(21, 24)) +
  geom_point(size=3,  aes(fill='State'), show.legend = F) + 
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(colour="State") +
  stat_ellipse(geom = 'path' , aes(color = State, group=State,linetype=factor(ps_gwwa_2122@sam_data$State)), level=0.5, lwd=0.6) +
  scale_linetype_manual(values=c(2,2), guide='none')
dev.off()

#alpha gw parental 21-22
#build basic linear model
alpha_gw2122<-alpha_r[which(rownames(alpha_r) %in% rownames(ps_gwwa_2122@sam_data)),]
lm_chao1_gw2122 <- lm(Chao1 ~ State + Year, data=alpha_gw2122)
#run a marginal test Anova (Type II)
lm_chao1_gw2122 %>% car::Anova() #year p=0.009
plot_richness(ps_gwwa_2122, x="State", color="Year",measures=c("Chao1","Shannon")) #may be effect of UP-3 of 4 samples collected in 2022

lm_shan_gw2122 <- lm(Shannon ~ State + Year, data=alpha_gw2122)
#run a marginal test Anova (Type II)
lm_shan_gw2122 %>% car::Anova() #year=0.010

lm_pd_gw2122 <- lm(PD ~ State + Year, data=alpha_gw2122)
#run a marginal test Anova (Type II)
lm_pd_gw2122 %>% car::Anova() #year, P=0.008


#----intermediates 2021-2022-----
#subset to intermediates
ps_int<-prune_samples(rownames(ps3k@sam_data)[ps3k@sam_data$admx=='Intermediate'], ps3k) #all intermediates
ps_int<-prune_taxa(taxa_sums(ps_int) > 0, ps_int) #get rid of zero-sum OTUs
table(ps_int@sam_data$Complete_plumage, useNA='always') #all with complete plumage

ps_int_2122<-prune_samples(rownames(ps_int@sam_data)[which(ps_int@sam_data$Year %in% c('2021','2022'))] , ps_int) 
ps_int_2122<-prune_taxa(taxa_sums(ps_int_2122) > 0, ps_int_2122) #get rid of zero-sum OTUs
table(ps_int_2122@sam_data$State, ps_int_2122@sam_data$Year, useNA='always')

#beta diversity for intermediates 21-22 
jac_int_2122<-phyloseq::distance(ps_int_2122, method='jaccard', binary=T)
bray_int_2122<-phyloseq::distance(ps_int_2122, method='bray')
uni_int_2122<-phyloseq::distance(ps_int_2122, method='unifrac')
wuni_int_2122<-phyloseq::distance(ps_int_2122, method='wunifrac')

md_int_2122<-data.frame(sample_data(ps_int_2122), row.names=rownames(sample_data(ps_int_2122)))
adonis2(bray_int_2122 ~ State+Year, data=md_int_2122, by='margin') # ns
adonis2(jac_int_2122 ~ State+Year, data=md_int_2122, by='margin') # state r2=0.05539, p=0.040
adonis2(uni_int_2122 ~ State+Year, data=md_int_2122, by='margin') #ns
adonis2(wuni_int_2122 ~ State+Year, data=md_int_2122, by='margin') # year r2=0.11857, p=0.024

permutest(betadisper(jac_int_2122, md_int_2122$State)) #F=0.3283,p=0.614
permutest(betadisper(wuni_int_2122, md_int_2122$Year)) #F=18.554,p=0.003**
plot(betadisper(wuni_int_2122, md_int_2122$Year), seg.lwd = 0.5, main='')

pdf('plots2/pcoa_int_state_jacc2122.pdf', height=3, width=3.5)
ordinate(ps_int_2122, "PCoA", "jaccard") %>%
  plot_ordination(ps_int_2122, ., shape='State', title = "Intermediates 2021-22, Jaccard", type='samples') +
  scale_fill_manual(values=c('gray','gray')) +  scale_shape_manual(values=c(21, 24)) +
  geom_point(size=3,  aes(fill='State'), show.legend = F) + 
  theme_bw() + theme(axis.text = element_text(size = 11)) + labs(colour="State") +
  stat_ellipse(geom = 'path' , aes(color = State, group=State,linetype=factor(ps_int_2122@sam_data$State)), level=0.5, lwd=0.6) +
  scale_linetype_manual(values=c(2,2), guide='none')
dev.off()

#alpha for intermediates 2021-2022
#build basic linear model
alpha_int_2122<-alpha_r[which(alpha_r$admx=='Intermediate' & alpha_r$Year %in% c('2021','2022')),]
table(alpha_int_2122$admx, alpha_int_2122$State, useNA = 'always')
lm_chao1_int2122 <- lm(Chao1 ~ State + Year, data=alpha_int_2122)
#run a marginal test Anova (Type II)
lm_chao1_int2122 %>% car::Anova() #ns
plot_richness(ps_int_2122, x="State", color="Species",measures=c("Chao1","Shannon")) #more ASVs in the north

lm_shan_int2122 <- lm(Shannon ~ State + Year, data=alpha_int_2122)
#run a marginal test Anova (Type II)
lm_shan_int2122 %>% car::Anova() #ns

lm_pd_int2122 <- lm(PD ~ State + Year, data=alpha_int_2122)
#run a marginal test Anova (Type II)
lm_pd_int2122 %>% car::Anova() #year sum sq=31.662, P=0.031


#####----------are there ASVs that differ between parentals?------
#model 1
maas_parental_all = Maaslin2(
  input_data = data.frame(ps_parental_all@otu_table), 
  input_metadata = data.frame(ps_parental_all@sam_data), 
  output = "maaslin2_output2/parental_all", 
  fixed_effects = c("Species","State","Year")) 
#model 2
maas_ps3k = Maaslin2(
  input_data = data.frame(ps3k@otu_table), 
  input_metadata = data.frame(ps3k@sam_data), 
  output = "maaslin2_output2/ps3k", 
  fixed_effects = c("Species","State","Year"),
  reference = c('Species,BWWA')) 
#model 3
maas_ps3k_plumage2 = Maaslin2(
  input_data = data.frame(ps3k@otu_table), 
  input_metadata = data.frame(ps3k@sam_data), 
  output = "maaslin2_output2/ps3k_plumage2", 
  fixed_effects = c("Plumage_score2",'State','Year'))
#model 4
maas_int_plumage2 = Maaslin2(
  input_data = data.frame(ps_int@otu_table), 
  input_metadata = data.frame(ps_int@sam_data), 
  output = "maaslin2_output2/int_plumage2", 
  fixed_effects = c("Plumage_score2",'State','Year'))

#compute carotenoid score sum
carot<-ps3k@sam_data
carot<-carot[,c('Species','admx','Nape','Back','Rump','Breast','Belly')] #pull out body carotenoids
carot$Nape<-as.numeric(carot$Nape)
carot$Back<-as.numeric(carot$Back)
carot$Rump<-as.numeric(carot$Rump)
carot$Breast<-as.numeric(carot$Breast)
carot$Belly<-as.numeric(carot$Belly)
carot$sum<-rowSums(carot[,3:7], na.rm = T) 
carot[,3:8]#need to turn sum into NA for birds with missing scores
carot$sum[!complete.cases(carot[,3:7])] <-NA #done (birds w/missing carot plumage have NA)
#add wing bar color, but convert to low score=yellow to match scale for body carotenoids
identical(rownames(carot), rownames(ps3k@sam_data))
carot$wbar<-5-as.numeric(ps3k@sam_data$W_bar) #max score is 5, so 5-w_bar reverses the scale
carot$sum2<-rowSums(carot[,c(3:7,9)], na.rm = T) #with wing bar
carot[,c(3:8,9:10)]#need to turn sum into NA for birds with missing scores
carot$sum2[!complete.cases(carot[,c(3:7,9)])] <-NA #done (birds w/missing carot plumage have NA)

#model 5
ps_int@sam_data$carot_sum2<-carot$sum2[match(rownames(ps_int@sam_data), rownames(carot))] 
maas_int_carot2 = Maaslin2(
  input_data = data.frame(ps_int@otu_table), 
  input_metadata = data.frame(ps_int@sam_data), 
  output = "maaslin2_output2/int_carot2", 
  fixed_effects = c("carot_sum2",'State','Year'))

#model 6
ps3k@sam_data$carot_sum2<-carot$sum2[match(rownames(ps3k@sam_data), rownames(carot))] 
maas_ps3k_carot3 = Maaslin2(
  input_data = data.frame(ps3k@otu_table), 
  input_metadata = data.frame(ps3k@sam_data),
  output = "maaslin2_output2/ps3k_carot3", 
  fixed_effects = c("carot_sum2",'State','Year'))


#---counts of ASVs that differ between species
# golden-wing taxon ca0365837eefeb9bf94c58dd59c0e8c3
rick<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='ca0365837eefeb9bf94c58dd59c0e8c3')]
rick_counts<-data.frame(ps3k@sam_data)
rick_counts$asv_count=as.vector(rick)
summary(rick_counts$asv_count[which(rick_counts$admx=='GWWA')]) #range 0-123, mean 8.4
summary(rick_counts$asv_count[which(rick_counts$admx=='BWWA')]) #range 0-8.00, mean 0.16
length(which(rick_counts$asv_count > 0)) #16 out of 123

rick_counts$admx<-factor(rick_counts$admx, levels=c('GWWA','Intermediate','BWWA')) #order for plot axis
pdf(file='plots2/counts_rick.pdf', height=3.8, width=3)
boxplot(rick_counts$asv_count~rick_counts$admx, xlab='', ylab='Counts of Rickettsia ASV',
        main='', boxwex=0.3, outline=F,ylim=c(0,140), col=c('gold','gray','cornflowerblue'))
stripchart(rick_counts$asv_count[which(rick_counts$admx=='GWWA')],vertical=T,pch=21,bg='gold', method='jitter',add=T)
stripchart(rick_counts$asv_count[which(rick_counts$admx=='Intermediate')],vertical=T,pch=21,bg='gray', method='jitter',add=T,at=2)
stripchart(rick_counts$asv_count[which(rick_counts$admx=='BWWA')],vertical=T,pch=21,bg='cornflowerblue', method='jitter',add=T,at=3)
dev.off()

ggplot(rick_counts, aes(x = Plumage_score2, y = asv_count)) + 
  geom_point(aes(fill = admx, shape = admx), size=3) + 
  scale_color_manual(values=c('black','black','black')) + scale_fill_manual(values=c('gold','gray','cornflowerblue')) +
  scale_shape_manual(values=c(21,21,21)) +
  #  geom_point(size=3, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2)


# bwwa-associated Kineococcus: b1dde0ed7ab1f036845b5643cfaa699e
kine<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='b1dde0ed7ab1f036845b5643cfaa699e')]
kine_counts<-data.frame(ps3k@sam_data)
kine_counts$asv_count=as.vector(kine)
summary(kine_counts$asv_count[which(kine_counts$admx=='GWWA')]) #range 0-27, mean 1
summary(kine_counts$asv_count[which(kine_counts$admx=='BWWA')]) #range 0-33, mean 2
length(which(kine_counts$asv_count > 0)) #20 out of 123
kine_counts[which(kine_counts$asv_count > 0),]

kine_counts$admx<-factor(kine_counts$admx, levels=c('GWWA','Intermediate','BWWA')) #order for plot axis, note not truly admix, but species
boxplot(kine_counts$asv_count~kine_counts$admx, xlab='', ylab='Counts of Kineococcus ASV',
        main='', boxwex=0.3, outline=F, ylim=c(0,35), col=c('gold','gray','cornflowerblue'))
stripchart(kine_counts$asv_count[which(kine_counts$admx=='GWWA')],vertical=T,pch=21,bg='gold', method='jitter',add=T)
stripchart(kine_counts$asv_count[which(kine_counts$admx=='Intermediate')],vertical=T,pch=21,bg='gray', method='jitter',add=T,at=2)
stripchart(kine_counts$asv_count[which(kine_counts$admx=='BWWA')],vertical=T,pch=21, bg='cornflowerblue', method='jitter',add=T,at=3)

ggplot(kine_counts, aes(x = Plumage_score2, y = asv_count)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2)

# brewster's associated ASV: dd070aa4a7da33679f1038d652140369, Jatrophihabitans
jatr<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='dd070aa4a7da33679f1038d652140369')]
jatr_counts<-data.frame(ps3k@sam_data)
jatr_counts$asv_count=as.vector(jatr)
summary(jatr_counts$asv_count[which(jatr_counts$admx=='GWWA')]) #range=0-25, mean 1
summary(jatr_counts$asv_count[which(jatr_counts$admx=='BWWA')]) #range 0-48, mean 1.56
summary(jatr_counts$asv_count[which(jatr_counts$admx=='Intermediate')]) #range 0-87, mean 7.69

length(which(jatr_counts$asv_count > 0)) #13 out of 123
jatr_counts[which(jatr_counts$asv_count > 0),]
table(jatr_counts$admx, useNA = 'always') #percentages for paper

jatr_counts$admx<-factor(jatr_counts$admx, levels=c('GWWA','Intermediate','BWWA')) #order for plot axis
pdf('plots2/counts_jatr.pdf', height=3.8,width=3)
boxplot(jatr_counts$asv_count~jatr_counts$admx, xlab='', ylab='Counts of Jatrophihabitans ASV',
        main='', boxwex=0.5, outline=F,ylim=c(0,90), col=c('gold','gray','cornflowerblue'),
        names=c('GW','Intermediate','BW'))
stripchart(jatr_counts$asv_count[which(jatr_counts$admx=='GWWA')],vertical=T,pch=21,bg='gold', method='jitter',add=T)
stripchart(jatr_counts$asv_count[which(jatr_counts$admx=='Intermediate')],vertical=T,pch=21,bg='gray', method='jitter',add=T,at=2)
stripchart(jatr_counts$asv_count[which(jatr_counts$admx=='BWWA')],vertical=T,pch=21,bg='cornflowerblue', method='jitter',add=T,at=3)
dev.off()

jatr_counts[which(jatr_counts$asv_count > 0),] #includes PA, MI, 2020-2022 -- cool result


####----check if species-specific bateria occur in same birds--are they *incompatible*??---
candidate<-data.frame(ps3k@sam_data)
candidate$kine<-kine_counts$asv_count
candidate$rick<-rick_counts$asv_count
candidate$jatr<-jatr_counts$asv_count
candidate[which(candidate$Species=='BRWA'),]
candidate[which(candidate$admx=='Intermediate'),c(48,52:54)]


#-----write otu table to biom format for picrust2
library(biomformat)
biom_ps3k<-as(otu_table(ps3k),"matrix")
biom_ps3k<-make_biom(data=biom_ps3k)
write_biom(biom_ps3k,"biom_ps3k.biom")
#write metadata to file
write.table(ps3k@sam_data, 'ps3k_metadata_forpicrust.tsv', quote=F, sep='\t')


###---scatter plots
#spingomonas
sphing<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='b5c8ed4496d79b916859402cdf11f0a7')]
sphing_counts<-data.frame(ps3k@sam_data)
sphing_counts$asv_count=as.vector(sphing)
sphing_counts_int<-sphing_counts[which(sphing_counts$admx=='Intermediate'),]

pdf(file='plots2/counts_sphing_scatter.pdf', height=4, width=3.5)
ggplot(sphing_counts_int, aes(x = Plumage_score2, y = asv_count)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2) +
  coord_cartesian(ylim = c(-0.2,max(sphing_counts_int$asv_count)))
dev.off()

#pseudomonas
psued<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='5aec9bd35889489e4a05c78a82358060')]
psued_counts<-data.frame(ps3k@sam_data)
psued_counts$asv_count=as.vector(psued)
psued_counts_int<-psued_counts[which(psued_counts$admx=='Intermediate'),]

pdf(file='plots2/counts_pseud_scatter.pdf', height=4, width=3.5)
ggplot(psued_counts_int, aes(x = Plumage_score2, y = asv_count)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2) +
  coord_cartesian(ylim = c(-0.2,max(psued_counts_int$asv_count)))
dev.off()

#Allorhizobium
allo<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='fc6622c636a5210293fb2873fc4761d9')]
allo_counts<-data.frame(ps3k@sam_data)
allo_counts$asv_count=as.vector(allo)
allo_counts_int<-allo_counts[which(allo_counts$admx=='Intermediate'),]

pdf(file='plots2/counts_allo_scatter.pdf', height=4, width=3.5)
ggplot(allo_counts_int, aes(x = Plumage_score2, y = asv_count)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2) +
  coord_cartesian(ylim = c(-0.2,max(allo_counts_int$asv_count)))
dev.off()

#Methylobacterium
meth<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='04ecfad5772d2e09a84a0f5ef460536c')]
meth_counts<-data.frame(ps3k@sam_data)
meth_counts$asv_count=as.vector(meth)
meth_counts_int<-meth_counts[which(meth_counts$admx=='Intermediate'),]

pdf(file='plots2/counts_meth_scatter.pdf', height=4, width=3.5)
ggplot(meth_counts_int, aes(x = Plumage_score2, y = asv_count)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2) +
  coord_cartesian(ylim = c(-0.2,max(meth_counts_int$asv_count)))
dev.off()

#microbacteriaceae
micro<-ps3k@otu_table[which(rownames(ps3k@otu_table)=='3ab43fcc04514ef20f0fe8689a8c5ea7')]
micro_counts<-data.frame(ps3k@sam_data)
micro_counts$asv_count=as.vector(micro)
micro_counts[which(micro_counts$asv_count>0 & micro_counts$carot_sum2>15),]

pdf(file='plots2/counts_micro_scatter.pdf', height=4, width=3.5)
ggplot(micro_counts, aes(x = carot_sum2, y = asv_count)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_text(size = 14)) +
  stat_smooth(method = "lm", col = "darkgray", lty=2) +
  coord_cartesian(ylim = c(-0.2,max(micro_counts$asv_count)))
dev.off()

##-----top taxa-----
topp3k<-as.data.frame(ps3k@tax_table)
table(topp3k$Phylum, useNA = 'always')
sort(table(topp3k$Phylum), decreasing=T)
sort(table(topp3k$Phylum)/sum(table(topp3k$Phylum)), decreasing=T) #top phyla by proportion

topp<-names(sort(table(topp3k$Phylum), decreasing=T)[1:8]) #top 8 phyla

glom3k<-tax_glom(ps3k, taxrank = "Phylum", NArm=F) #glom by phylum
mglom3k<-psmelt(glom3k) #31 OTUs x 123 individuals=3813 rows
mglom3k$Phylum[which(!(mglom3k$Phylum %in% topp))]<-'other' #change non-common phyla to other 
table(mglom3k$Phylum)

ggplot(mglom3k, aes(x = sample_Species, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('cornsilk3','darkblue',
                                        'cornflowerblue','coral3','darkgreen',
                                        'darkolivegreen3','black','darkgoldenrod','gold2')) +
                                          ggtitle('Rarefied dataset') +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  facet_wrap(vars(State))

###---how many OTUs are unique to each parental species?-----
#parentals--all
library(VennDiagram)
venn.diagram(
  x = list(rownames(ps_bwwa_all@otu_table),rownames(ps_gwwa_all@otu_table)),
  category.names = c("BWWA" , "GWWA"), #parentals only
  filename = 'venn_par.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c('cornflowerblue','gold')
) #343 OTUs overlap in both species

#all individuals--by admx
venn.diagram(
  x = list(rownames(ps_bwwa_all@otu_table),rownames(ps_gwwa_all@otu_table),rownames(ps_int@otu_table)),
  category.names = c("BWWA" , "GWWA","Intermediate"), #including intermediates
  filename = 'venn_par_int.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c('cornflowerblue','gold','gray')
) #343 OTUs overlap in both species

##all individuals mi--by admx
#michigan only nas removed
ps3k_nar_mi<-prune_samples(rownames(ps3k_nar@sam_data) %in% rownames(ps_mi_all@sam_data), ps3k_nar) 
ps3k_nar_mi<-prune_taxa(taxa_sums(ps3k_nar_mi) > 0, ps3k_nar_mi) #get rid of zero-sum OTUs
table(ps3k_nar_mi@sam_data$admx, ps3k_nar_mi@sam_data$Species)
#make data.frame of metadata
md3k_nar_mi<-data.frame(sample_data(ps3k_nar_mi), row.names=rownames(sample_data(ps3k_nar_mi)))

table(ps3k_nar_mi@sam_data$admx)
#subset into groups for venn
ps_bwwa_mi<-prune_samples(rownames(ps3k_nar_mi@sam_data)[which(ps3k_nar_mi@sam_data$admx=='BWWA')], ps3k_nar_mi) 
ps_bwwa_mi<-prune_taxa(taxa_sums(ps_bwwa_mi) > 0, ps_bwwa_mi) #get rid of zero-sum OTUs

ps_gwwa_mi<-prune_samples(rownames(ps3k_nar_mi@sam_data)[which(ps3k_nar_mi@sam_data$admx=='GWWA')], ps3k_nar_mi) 
ps_gwwa_mi<-prune_taxa(taxa_sums(ps_gwwa_mi) > 0, ps_gwwa_mi) #get rid of zero-sum OTUs

ps_int_mi<-prune_samples(rownames(ps3k_nar_mi@sam_data)[which(ps3k_nar_mi@sam_data$admx=='Intermediate')], ps3k_nar_mi) 
ps_int_mi<-prune_taxa(taxa_sums(ps_int_mi) > 0, ps_int_mi) #get rid of zero-sum OTUs

venn.diagram(
  x = list(rownames(ps_bwwa_mi@otu_table),rownames(ps_gwwa_mi@otu_table),rownames(ps_int_mi@otu_table)),
  main='Great Lakes HZ',
  category.names = c("BWWA" , "GWWA","Intermediate"), #including intermediates
  filename = 'venn_par_mi.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c('cornflowerblue','gold','gray')
) #49 OTUs overlap in 3 groups

##all individuals pa--by admx
#pennsylvania only nas removed
ps3k_nar_pa<-prune_samples(rownames(ps3k_nar@sam_data) %in% rownames(ps_pa_all@sam_data), ps3k_nar) 
ps3k_nar_pa<-prune_taxa(taxa_sums(ps3k_nar_pa) > 0, ps3k_nar_pa) #get rid of zero-sum OTUs
table(ps3k_nar_pa@sam_data$admx, ps3k_nar_pa@sam_data$Species)
#make data.frame of metadata
md3k_nar_pa<-data.frame(sample_data(ps3k_nar_pa), row.names=rownames(sample_data(ps3k_nar_pa)))

table(ps3k_nar_pa@sam_data$admx)
#subset into groups for venn
ps_bwwa_pa<-prune_samples(rownames(ps3k_nar_pa@sam_data)[which(ps3k_nar_pa@sam_data$admx=='BWWA')], ps3k_nar_pa) 
ps_bwwa_pa<-prune_taxa(taxa_sums(ps_bwwa_pa) > 0, ps_bwwa_pa) #get rid of zero-sum OTUs

ps_gwwa_pa<-prune_samples(rownames(ps3k_nar_pa@sam_data)[which(ps3k_nar_pa@sam_data$admx=='GWWA')], ps3k_nar_pa) 
ps_gwwa_pa<-prune_taxa(taxa_sums(ps_gwwa_pa) > 0, ps_gwwa_pa) #get rid of zero-sum OTUs

ps_int_pa<-prune_samples(rownames(ps3k_nar_pa@sam_data)[which(ps3k_nar_pa@sam_data$admx=='Intermediate')], ps3k_nar_pa) 
ps_int_pa<-prune_taxa(taxa_sums(ps_int_pa) > 0, ps_int_pa) #get rid of zero-sum OTUs

venn.diagram(
  x = list(rownames(ps_bwwa_pa@otu_table),rownames(ps_gwwa_pa@otu_table),rownames(ps_int_pa@otu_table)),
  main='Appalachian HZ',
  category.names = c("BWWA" , "GWWA","Intermediate"), #including intermediates
  filename = 'venn_par_pa.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c('cornflowerblue','gold','gray')
) #174 OTUs overlap in 3 groups

#manually plot to scale (used numbers from above)
library(eulerr)
pdf('plots2/venn_appalachian.pdf',height=3.5, width=4)
plot(euler(c(A=1193, B=559, C=1513,'A&B'=67, 'B&C'=68, 'A&C'=189, 'A&B&C'=174)),
     fill = c('cornflowerblue','gold','gray'), main='N ASVs, Appalachian HZ',
     alpha=0.7,quantities = list(cex =1.2), legend = list(labels = c("BWWA", "GWWA",'Admixed')))
dev.off()
pdf('plots2/venn_greatlakes.pdf',height=3.5, width=4)
plot(euler(c(A=126, B=352, C=223,'A&B'=23, 'B&C'=33, 'A&C'=13, 'A&B&C'=49)),
     fill = c('cornflowerblue','gold','gray'), main='N ASVs, Great Lakes HZ',
     alpha=0.7,quantities = list(cex =1.2), legend = list(labels = c("BWWA", "GWWA",'Admixed')), pt.cex=1.2)
dev.off()

#manual check!
otus_bwwa_all<-rownames(ps_bwwa_all@otu_table)
otus_gwwa_all<-rownames(ps_gwwa_all@otu_table)
#otus unique to bwwa
otus_bwwa_only<-setdiff(otus_bwwa_all, otus_gwwa_all)
length(otus_bwwa_only) #1399
otus_gwwa_only<-setdiff(otus_gwwa_all, otus_bwwa_all)
length(otus_gwwa_only) #899
#otus shared by both species
otus_shared<-intersect(otus_bwwa_all, otus_gwwa_all)
length(otus_shared) #343

#unique ASVs
#make a ps object with only the unique bwwa ASVs, include all parental birds
#BWWA
ps_bw_otus<-subset_taxa(ps_bwwa_all, rownames(ps_bwwa_all@otu_table) %in% otus_bwwa_only)
ps_bw_otus<-prune_taxa(taxa_sums(ps_bw_otus) > 0, ps_bw_otus) #get rid of zero-sum OTUs
bwotu <- as.data.frame(ps_bw_otus@otu_table)
bwotu$asv <- rownames(ps_bw_otus@otu_table) #add column of asv ids
#calculate prevalence of private bwwa otus
bwotu[bwotu==0]<-NA #turn zeros to NA
bwotu$sumNA<-rowSums(is.na(bwotu)) #add column of sums of NAs per row
summary(bwotu$sumNA) #min is 44, meaning present in (50-44=6) 6/50, or 12% of individuals
head(bwotu[order(bwotu$sumNA),])
hist(bwotu$sumNA) #range 1-6 individuals
bwotu$prev<-c((50-bwotu$sumNA)/50) #add column of prevalence (proportion of individuals with this ASV)
head(bwotu$prev[order(bwotu$sumNA)], 100)

##make a ps object with only the unique gwwa ASVs, include all parental birds
ps_gw_otus<-subset_taxa(ps_gwwa_all, rownames(ps_gwwa_all@otu_table) %in% otus_gwwa_only)
ps_gw_otus<-prune_taxa(taxa_sums(ps_gw_otus) > 0, ps_gw_otus) #get rid of zero-sum OTUs

gwotu <- as.data.frame(ps_gw_otus@otu_table)
gwotu$asv <- rownames(ps_gw_otus@otu_table) #add column of asv ids
#calculate prevalence of private gwwa otus
gwotu[gwotu==0]<-NA #turn zeros to NA
gwotu$sumNA<-rowSums(is.na(gwotu)) #add column of sums of NAs per row
summary(gwotu$sumNA) #min is 27, meaning present in (29-27=2) 2/29, or 7% of individuals
head(gwotu[order(gwotu$sumNA),])
hist(gwotu$sumNA) #range 1-2 individuals
gwotu$prev<-c((29-gwotu$sumNA)/29) #add column of prevalence (proportion of individuals with this ASV)
head(gwotu$prev[order(gwotu$sumNA)], 100)

##make a ps object with only the shared bwwa-gwwa ASVs, include all parental birds
ps_shared_otus<-subset_taxa(ps_parental_all, rownames(ps_parental_all@otu_table) %in% otus_shared)
ps_shared_otus<-prune_taxa(taxa_sums(ps_shared_otus) > 0, ps_shared_otus) #get rid of zero-sum OTUs

sharedotu <- as.data.frame(ps_shared_otus@otu_table)
sharedotu$asv <- rownames(ps_shared_otus@otu_table) #add column of asv ids
#calculate prevalence of share otus
sharedotu[sharedotu==0]<-NA #turn zeros to NA
sharedotu$sumNA<-rowSums(is.na(sharedotu)) #add column of sums of NAs per row
summary(sharedotu$sumNA) #min is 28, meaning present in (79-28=51) 51/79, or 65% of individuals
head(sharedotu[order(sharedotu$sumNA),])
hist(sharedotu$sumNA) 
sharedotu$prev<-c((79-sharedotu$sumNA)/79) #add column of prevalence (proportion of individuals with this ASV)
head(sharedotu$prev[order(sharedotu$sumNA)], 100)
head(sharedotu[order(sharedotu$sumNA),])
tail(sharedotu[order(sharedotu$sumNA),])


#histograms of ASV prevalences 
hist(79-sharedotu$sumNA, breaks=50, xlab='N individuals', ylab='N ASVs', main='Shared ASVs')
hist(29-gwotu$sumNA, breaks=50, xlab='N individuals', ylab='N ASVs', main='GWWA private ASVs')
hist(50-bwotu$sumNA, breaks=50, xlab='N individuals', ylab='N ASVs', main='BWWA private ASVs')

#what proportion of reads are shared, versus private?
sum(taxa_sums(ps_shared_otus)) #181375
sum(taxa_sums(ps_shared_otus))/sum(taxa_sums(ps_parental_all)) #75% of all reads are shared OTUs
summary(sharedotu$prev) #6.6%
sum(taxa_sums(ps_bw_otus)) #41091
sum(taxa_sums(ps_bw_otus))/sum(taxa_sums(ps_parental_all)) #17% of all reads are private BWWA OTUs
summary(bwotu$prev) #2.1%
sum(taxa_sums(ps_gw_otus)) #20538
sum(taxa_sums(ps_gw_otus))/sum(taxa_sums(ps_parental_all)) #8% of all reads are private GWWA OTUs
summary(gwotu$prev) #3.6%
