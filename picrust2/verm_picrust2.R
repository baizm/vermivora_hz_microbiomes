setwd("~/Documents/Toews_Lab/vermivora_microbiome/picrust2")

library(Maaslin2)
library(data.table)

load('verm_picrust2.RData')
#read in metadata
m<-read.csv("../qiime_output/phyloseq/ps3k_metadata_forpicrust2.tsv", sep='\t', header=T, check.names = F) 
m$Year<-as.character(m$Year)

p2EC = as.data.frame(fread('./picrust2_out_pipeline/EC/pred_metagenome_unstrat_descrip.tsv'))
rownames(p2EC) = p2EC$"function"
ec = as.matrix(p2EC[,c(-1,-2)])
ec = round(ec)
#look at descriptions for EC
ed<-read.csv('picrust2_out_pipeline/EC/EC_descrip.tsv', sep='\t', header=T, check.names = F)

p2KO = as.data.frame(fread('./picrust2_out_pipeline/KO/pred_metagenome_unstrat_descrip.tsv'))
rownames(p2KO) = p2KO$"function"
ko = as.matrix(p2KO[,c(-1,-2)])
ko = round(ko)
#look at descriptions for KO
kd<-read.csv('picrust2_out_pipeline/KO/KO_descrip.tsv', sep='\t', header=T, check.names = F)


#-----run maaslin2 on KO counts------
#subset to intermediates to look at plumage assocations
m_int<-m[which(m$admx=='Intermediate'),]
ko_int<-ko[,rownames(m_int)]
ko_int<-data.frame(ko_int)
dim(ko_int) 

maas_pi_plumage_int2 = Maaslin2(
  input_data = ko_int,
  input_metadata = m_int,
  output = "redo/maaslin2_output_all/pi_plumage_int", 
  fixed_effects = c("Plumage_score2",'State','Year'))

maas_pi_carot_int2 = Maaslin2(
  input_data = ko_int,
  input_metadata = m_int,
  output = "redo/maaslin2_output_all/pi_carot_int", 
  fixed_effects = c("carot_sum2",'State','Year'))

##how many KOs in intermediates?
dim(ko) #7579 KOs in total
dim(ko_int[rowSums(ko_int[])>0,]) #6986 KOs in intermediates (rows that sum to greater than zero)
dim(ko_int[rowSums(ko_int[])>0,])/dim(ko) #92% of the KOs present in intermediates


#-----run maaslin2 on all EC counts------
#subset to intermediates to look at plumage assocations
ec_int<-ec[,rownames(m_int)]
ec_int<-data.frame(ec_int)

maas_pi_plumage_int2_ec = Maaslin2(
  input_data = ec_int,
  input_metadata = m_int,
  output = "redo/maaslin2_output_ec/all/pi_plumage_int", 
  fixed_effects = c("Plumage_score2","State","Year"))

maas_pi_carot_int2_ec = Maaslin2(
  input_data = ec_int,
  input_metadata = m_int,
  output = "redo/maaslin2_output_ec/all/pi_carot_int", 
  fixed_effects = c("carot_sum2","State","Year"))

##how many ECs in intermediates?
dim(ec) #2314 total
dim(ec_int[rowSums(ec_int[])>0,]) #2158 ECs in intermediates (rows that sum to greater than zero)
dim(ec_int[rowSums(ec_int[])>0,])/dim(ec) #93% of the ECs present in intermediates

