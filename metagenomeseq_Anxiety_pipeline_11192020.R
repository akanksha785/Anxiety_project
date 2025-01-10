library(readxl)
library(dada2)
library(DECIPHER)
library(phangorn)
library(phyloseq)
library(gridExtra)
library(dplyr)
library(ape)
library(plyr)                                     #Check dada2 sequences obtained from otu table
library(ggplot2)
library(tidyverse)
library(csv)
library(decontam)
library(metagenomeSeq)
library(splinectomeR)
library(biomformat)
library(superheat)############################heatmaps
library(reshape2)
library(tibble)
library(tidyr)
library(vegan)
library(gplots) ##########################################for heatmap.2 function#######################
library(gdata) ##########merging columns into one
library(igraph)
library(splinectomeR)
library(ALDEx2)
library(Maaslin2)
library(HMP)
library(dendextend)
library(ppcor)     #########spearman partial correlation coefficient
library(rbin)
library(stats)
library(ggpubr)
library(nloptr)
library(coda)
library(gss)
library(remotes)
library(devtools)
library(SpiecEasi)
library(randomForest)
library(caret)
library(RColorBrewer) ############color palatte
library(plotrix)
library(robCompositions)
library(Maaslin2)
library(psych)    ##################for correlation test###################
library(rmcorr)
library(clr)
library(taxonomizr)
library(compositions)   ###########################for clr transformation#########################
#install_github("hallucigenia-sparsa/seqtime") 
#library(seqtime)
devtools::install_github("zdk123/SpiecEasi")

library(devtools)
devtools::install_github("biobakery/maaslin2")


path<-'C:/Akanksha/LAB WORK/fastq_metagenomeseq'
list.files(path)

# Identify all the files within the above directory
#fns<- list.files(path, recursive = TRUE)


fnFs<-sort(list.files(path, pattern = "_L001_R1_001.fastq"))

fnRs<-sort(list.files(path, pattern = "_L001_R2_001.fastq"))


##Extracting sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names<-paste("Sample",sample.names,sep="")


plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtfwd_df <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtrvs_df<- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtfwd_df) <- sample.names
names(filtrvs_df) <- sample.names

#See quality of samples

out<-filterAndTrim(fnFs,filtfwd_df,fnRs,filtrvs_df,verbose=TRUE,trimLeft=20)####### make multihread = TRue

errF <- learnErrors(filtfwd_df, nbases=1e8, multithread=TRUE,randomize=TRUE,verbose=1)

# Learn reverse error rates
errR <- learnErrors(filtrvs_df, nbases=1e8, multithread=TRUE,randomize=TRUE,verbose=1)

dadaFwd <- dada(filtfwd_df, err=errF, multithread = TRUE)
dadaRvs <- dada(filtrvs_df, err=errR, multithread = TRUE)
dadaFwd[[1]]

mergers <- mergePairs(dadaFwd, filtfwd_df, dadaRvs, filtrvs_df, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


seqtab <- makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
dim(seqtab)         
table(nchar(getSequences(seqtab)))

###Output of table shows 5 sequences are 230 basepairs in length, 20 sequences are 231 basepairs in length######
####This DADA2 sequence table is analogous to a traditional OTU table#######################


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)   
sum(seqtab.nochim)/sum(seqtab)

###91 % of reads are non-chimeric############

for (i in 1:nrow(seqtab.nochim)) {
  for (j in 1:ncol(seqtab.nochim)) {
    if(seqtab.nochim[i,j]<=10){
      seqtab.nochim[i,j]<-0
    }
    
  }
  
}

dim(seqtab.nochim)

seqtabnochim1<-seqtab.nochim[rowSums(seqtab.nochim)>0,]
dim(seqtabnochim1)
seqtabnochim2<-seqtabnochim1[,colSums(seqtabnochim1)>0]

dim(seqtabnochim2)  
#seqtabnochim2<-t(seqtabnochim2)
table(nchar(getSequences(seqtabnochim2)))

#seqs_length<-rownames(seqtabnochim2)
#seqtabnochim3<-apply(seqs_length,2,function(x){(nchar(getSequences(x)))})

taxa <- assignTaxonomy(seqtabnochim2, "C:/Akanksha/LAB WORK/fastq_metagenomeseq/Silva/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

taxa_1<- addSpecies(taxa,"C:/Akanksha/LAB WORK/fastq_metagenomeseq/Silva/silva_species_assignment_v138.fa.gz")

taxonomy<-array(taxa_1)
taxonomy0<-data.frame(taxonomy)

#Dividing taxonomy into kingdom,phylum,family,order,class,genus,species

#Kingdom<-taxonomy0[1:3298, ]
#Phylum<-taxonomy0[3299:6596, ]
#Class<-taxonomy0[6597:9894, ]
#Order<-taxonomy0[9895:13192, ]
#Family<-taxonomy0[13193:16490, ]
#Genus<-taxonomy0[16491:19788, ]
#Species<-taxonomy0[19789:23086, ]

mapfile0<-read.csv('C:/Akanksha/LAB WORK/fastq_metagenomeseq/mapFile_UNC_anx_04082021.csv', header = TRUE, stringsAsFactors = FALSE)
#sample_EPDS<-mapfile0$SC0
#sample_PSS<-mapfile0$SC4
#sample_GAD7<-mapfile0$SC5

#samdf<-mapfile0[1:72,]
samdf1 <- data.frame(mapfile0[,-1], row.names = mapfile0[,1])
samdf2<-sample_data(samdf1)

seqs<-getSequences(seqtabnochim2)
names(seqs)<-seqs
alignment<-AlignSeqs(DNAStringSet(seqs), anchor = NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR_1 <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
#tree <- read.tree(file = "fitGTR_1$tree")

ps<-phyloseq(tax_table(taxa_1),sample_data(samdf2),otu_table(seqtabnochim2,taxa_are_rows = FALSE), phy_tree(fitGTR_1$tree))
rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL)

ps0<- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))  
table(phyloseq::tax_table(ps0)[, "Phylum"])
plot_bar(ps0, fill = "Phylum") + theme(legend.position="bottom")

tree<-plot_tree(ps0, nodelabf = nodeplotblank, color = NULL,
                shape = NULL, label.tips = "Phylum", sizebase = 5,
                text.size = 4, base.spacing = 0.01,
                ladderize = FALSE, plot.margin = 0.1)
plot(tree)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps0),
                      tax_table(ps0))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

filterPhyla = c("Campilobacterota"," Cyanobacteria"," Euryarchaeota","Fusobacteriota","Spirochaetota","Synergistota","Patescibacteria")              ##Removing these two phyla because their mean is too low and also its sum

###########ps2--------->Removing phyla because their mean and sum  is too low.############

tax_table_ps0<-data.frame(tax_table(ps0))
tax_table_ps0 <-tax_table_ps0[(!(tax_table_ps0$Phylum=="Cyanobacteria") & !(tax_table_ps0$Phylum=="Euryarchaeota")& 
               !(tax_table_ps0$Phylum=="Campilobacterota")& 
                 !(tax_table_ps0$Phylum=="Fusobacteriota")& 
                 !(tax_table_ps0$Phylum=="Spirochaetota")& 
                 !(tax_table_ps0$Phylum=="Synergistota")& 
                 !(tax_table_ps0$Phylum=="Patescibacteria")),]


#rownames(tax_table_ps0)<-paste("Seq",1:nrow(tax_table_ps0), sep="_")
tax_table_ps0<-as.matrix(tax_table_ps0)

otu_table_ps0<-data.frame(otu_table(ps0))
#colnames(otu_table_ps0)<-paste("Seq",1:ncol(otu_table_ps0), sep="_")
otu_table_ps0<-t(otu_table_ps0)
otu_table_ps0<-data.frame(otu_table_ps0)
seqs_ps0<-rownames(tax_table_ps0)
otu_table_ps0<-otu_table_ps0[c(seqs_ps0),]
otu_table_ps0<-t(otu_table_ps0)

samdf3<-data.frame(sample_data(ps0))

ps1<-phyloseq(tax_table(tax_table_ps0),sample_data(samdf3),otu_table(otu_table_ps0, taxa_are_rows = FALSE))


#ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)

# Subset to the remaining phyla
prevdf_1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))

#rownames(prevdf_1)<-paste("Seq",1:nrow(prevdf_1), sep="_")

# Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(ps1)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf_1)[(prevdf_1$Prevalence >= prevalenceThreshold)]

##########ps3-------->keeping taxa which is more than 10% of total samples####################################
ps2 = prune_taxa(keepTaxa, ps1)
#write.csv(otu_table(ps2), file = "otu_table.csv")


prevdf_2= apply(X = otu_table(ps2),
                 MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

prevdf_2 = data.frame(Prevalence = prevdf_2,
                      TotalAbundance = taxa_sums(ps2),
                      tax_table(ps2))

prevdf_2<-subset(prevdf_2, Phylum %in% get_taxa_unique(ps2, "Phylum"))

TotalAbundance<-prevdf_2$TotalAbundance

abund_phyla_plot<-ggplot(prevdf_2, aes(TotalAbundance, Prevalence ,color=Phylum)) +geom_point(size = 4)+
  scale_x_log10()+ xlab("Total Abundance") + ylab("ASV Prevalence") +
  facet_wrap(~Phylum)+theme(legend.position="none")


abund_phyla_plot<-abund_phyla_plot+theme(aspect.ratio = 1)
#abund_phyla_plot<-abund_phyla_plot+theme_bw()
plot(abund_phyla_plot,cex.lab = 0.5)

write.csv(otu_table_ps0,"Abundance.csv")

set.seed(100)
featureData =data.frame(tax_table(ps2))
rownames(featureData)<-paste("Seq",1:nrow(featureData), sep="_")
matrixData<-matrix(otu_table(ps2),ncol=ncol(otu_table(ps2)))
rownames(matrixData)<-rownames(otu_table(ps2))
colnames(matrixData)<-colnames(otu_table(ps2))
matrixData<-t(matrixData)
rownames(matrixData)<-paste("Seq",1:nrow(matrixData), sep="_")

metadata<-sample_data(ps2)[match(colnames(matrixData),rownames(sample_data(ps2))), ]
metadata<-metadata[match(colnames(matrixData),rownames(metadata)),]

otus_metagenomeSeq<-newMRexperiment(matrixData, 
                                    phenoData = AnnotatedDataFrame(metadata), 
                                    featureData = AnnotatedDataFrame(featureData))

otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq,present=1, depth = 1)
p = cumNormStatFast(otus_metagenomeSeq_filter)    #percentile by which to normalize counts
anxietyData = cumNorm(otus_metagenomeSeq, p = p)
otu_norm <- metagenomeSeq::MRcounts(anxietyData, norm = TRUE, log = FALSE)
otu_norm<-apply(otu_norm,c(1,2),function(x){round(x,0)})
otu_norm<-t(otu_norm)

#samdf3<-data.frame(sample_data(ps2))
metadata$log_SC5<-log((metadata$SC5)+1)
#samdf3<-metadata[match(colnames(otu_norm),rownames(metadata)),]
featureData<-as.matrix(featureData)

ps_norm<-phyloseq(tax_table(featureData),sample_data(metadata),otu_table(otu_norm,taxa_are_rows = FALSE))

#################################################################################################################

tab<- microbiome::alpha(ps_norm, index = c("Shannon"))
shannon<-tab$diversity_shannon
metadata<-cbind.data.frame(metadata,shannon)
#SC5<-metadata$SC5
#Age<-metadata$Age
#Years_of_Education<-metadata$Years.of.Education
#BMI<-metadata$BMI
#normalized_gestational_weeks<-metadata$Normalized_Gestationalweeks
write.csv(metadata, "metadata.csv")

#featureData<-as.matrix(featureData)
metadata_v1<-metadata[(metadata$Visit=="T1"),]
#colnames_v1<-rownames(metadata_v1)
#otu_norm<-t(otu_norm)
otu_norm_v1<-otu_norm[match(rownames(metadata_v1),rownames(otu_norm)),]
#matrixData_v1<-t(matrixData_v1)
otu_norm_v1<-otu_norm_v1[,colSums(otu_norm_v1)>0]
otu_norm_v1<-otu_norm_v1[rowSums(otu_norm_v1)>0,]
otu_norm_v1<-as.matrix(otu_norm_v1)

taxa_norm_v1<-featureData[match(colnames(otu_norm_v1),rownames(featureData)),]

metadata_v1_alpha<-metadata_v1[!(metadata_v1$Subject==1),]
otu_norm_v1_alpha<-otu_norm[match(rownames(metadata_v1_alpha),rownames(otu_norm)),]
otu_norm_v1_alpha<-otu_norm_v1_alpha[,colSums(otu_norm_v1_alpha)>0]
otu_norm_v1_alpha<-otu_norm_v1_alpha[rowSums(otu_norm_v1_alpha)>0,]
otu_norm_v1_alpha<-as.matrix(otu_norm_v1_alpha)
taxa_norm_v1_alpha<-as.matrix(featureData[match(colnames(otu_norm_v1_alpha),rownames(featureData)),])



metadata_v2<-metadata[(metadata$Visit=="T2"),]
#colnames_v2<-rownames(metadata_v2)
#otu_norm<-t(otu_norm)
otu_norm_v2<-otu_norm[match(rownames(metadata_v2),rownames(otu_norm)),]
#matrixData_v1<-t(matrixData_v1)
otu_norm_v2<-otu_norm_v2[,colSums(otu_norm_v2)>0]
otu_norm_v2<-otu_norm_v2[rowSums(otu_norm_v2)>0,]
otu_norm_v2<-as.matrix(otu_norm_v2)

taxa_norm_v2<-featureData[match(colnames(otu_norm_v2),rownames(featureData)),]

metadata_v4<-metadata[(metadata$Visit=="T4"),]

metadata_v4_alpha<-metadata_v4[!(metadata_v4$Subject==25),]
#colnames_v4<-rownames(metadata_v4)
#otu_norm<-t(otu_norm)


otu_norm_v4<-otu_norm[match(rownames(metadata_v4),rownames(otu_norm)),]
otu_norm_v4<-otu_norm_v4[,colSums(otu_norm_v4)>0]
otu_norm_v4<-otu_norm_v4[rowSums(otu_norm_v4)>0,]
otu_norm_v4<-as.matrix(otu_norm_v4)



otu_norm_v4_alpha<-otu_norm[match(rownames(metadata_v4_alpha),rownames(otu_norm)),]
#matrixData_v1<-t(matrixData_v1)
#otu_norm_v4<-as.matrix(otu_norm_v4)
#Since one of the rowsums is 0, needed to do this additional step in order to make the fit model work
otu_norm_v4_alpha<-otu_norm_v4_alpha[,colSums(otu_norm_v4_alpha)>0]
otu_norm_v4_alpha<-otu_norm_v4_alpha[rowSums(otu_norm_v4_alpha)>0,]
otu_norm_v4_alpha<-as.matrix(otu_norm_v4_alpha)

taxa_norm_v4_alpha<-as.matrix(featureData[match(colnames(otu_norm_v4_alpha),rownames(featureData)),])

taxa_norm_v4<-featureData[match(colnames(otu_norm_v4),rownames(featureData)),]
#race<-k1$Race
#SC5_1<-SC5+1
#log_SC5<-log(SC5_1)
#BMI<-k1$BMI
#k1<-cbind.data.frame(k1,log_SC5)

#Using spearman correlation to find basic correlation between socio-demographics and anxiety scores.
#Getting warning:Cannot compute exact p-value with ties: because  spearman correlation coefficient considers the rank of values, 
#the correlation test ignores the same ranks to find the p-values as a result we get the warning

model_1<-cor.test(SC5,Age,method = "spearman", data=metadata,exact = FALSE )
model_1

model_2<-cor.test(SC5,Years_of_Education,method = "spearman",data=metadata,exact = FALSE )
model_2

model_3<-cor.test(SC5,BMI,method = "spearman",data=metadata,exact = FALSE)
model_3

model_4<-cor.test(SC5,normalized_gestational_weeks,method = "spearman",data=metadata,exact = FALSE)
model_4

#No statistical significant difference in the GAD-7(log) values between whites and non-whites.

metadata$WhitesVsNOn_whites<-rep(0, nrow(metadata))
metadata$WhitesVsNOn_whites[!(metadata$Race==0)]<-1
model_5<-lm(SC5~WhitesVsNOn_whites, data = metadata)
summary(model_5)

#Statistical significant difference in the  GAD-7 values between T1 and T2 visit

data_1<-metadata[!(metadata$Visit=="T4"),]
model_6<-aov(SC5~Visit,data = data_1)
summary(model_6)

# No Statistical significant difference in the  GAD-7 values between T2 and T4 visit

data_2<-metadata[!(metadata$Visit=="T1"),]
model_7<-aov(SC5~Visit, data = data_2) 
summary(model_7)

#No Statistical significant difference in the  GAD-7 values between T1 and T4 visit

data_3<-metadata[!(metadata$Visit=="T2"),]
model_8<-aov(SC5~Visit,data = data_3)
summary(model_8)
##############Statistically significant figures##########################

fig_1 <- ggplot(metadata, aes(x=Visit, y=SC5,fill=Visit)) + 
  geom_violin()+
  geom_jitter(shape=16,   position=position_jitterdodge(jitter.width = 0.1, 
                                                        dodge.width = 0.9, seed = NA)) +labs(x="Visit", y = "GAD-7 score")
fig_1<-fig_1 + scale_fill_brewer(palette="Dark2")+guides(fill=guide_legend(title='Visit'))+theme_bw()+theme(aspect.ratio = 1)
fig_1<-fig_1+theme(text = element_text(size = 20))
plot(fig_1)

fig_2 <- ggplot(metadata, aes(x=BMI, y=SC5)) + geom_point(size=4,color="lightslateblue")+labs(x="BMI", y = "GAD-7 score")
fig_2<-fig_2+geom_smooth(method = "lm",col="black",size=2)+theme_bw()+theme(aspect.ratio=1)
fig_2<-fig_2+theme(text = element_text(size = 20))
plot(fig_2)
####################################################################################################################
#Visit 1
#Participant Characteristics
participant_ch_v1<-metadata[(metadata$Visit=="T1"),]
anx_low_v1<-count(participant_ch_v1$Type=="Low")
anx_mild_v1<-count(participant_ch_v1$Type=="Mild")
anx_moderate_v1<-count(participant_ch_v1$Type=="Moderate")
anx_severe_v1<-count(participant_ch_v1$Type=="Severe")

gad_participant_ch_v1<-rbind(anx_low_v1[2,],anx_mild_v1[2,],anx_moderate_v1[2,],anx_severe_v1[2,])

gad_participant_ch_v1<-gad_participant_ch_v1[,2]

categories<-c("Low","Mild","Moderate","Severe")

gad_participant_ch_v1<-cbind.data.frame(categories,gad_participant_ch_v1)

names(gad_participant_ch_v1)[2] <- "Count"
gad_participant_ch_v1$Total<-rep(23, nrow(gad_participant_ch_v1))

gad_participant_ch_v1$percent_v1<-(gad_participant_ch_v1$Count/gad_participant_ch_v1$Total)*100

fig_3 <- ggplot(gad_participant_ch_v1, aes(x=categories, y=percent_v1,fill=categories)) + 
  geom_col()+labs(x="Anxiety categories", y = "Percentage")
fig_3<-fig_3+scale_fill_brewer(palette="Dark2")
plot(fig_3)

#Visit 2

participant_ch_v2<-metadata[(metadata$Visit=="T2"),]
anx_low_v2<-count(participant_ch_v2$Type=="Low")
anx_mild_v2<-count(participant_ch_v2$Type=="Mild")
anx_moderate_v2<-count(participant_ch_v2$Type=="Moderate")
anx_severe_v2<-count(participant_ch_v2$Type=="Severe")

gad_participant_ch_v2<-rbind(anx_low_v2[2,],anx_mild_v2[2,],anx_moderate_v2[2,],anx_severe_v2[2,])

gad_participant_ch_v2<-gad_participant_ch_v2[,2]

categories<-c("Low","Mild","Moderate","Severe")

gad_participant_ch_v2<-cbind.data.frame(categories,gad_participant_ch_v2)

names(gad_participant_ch_v2)[2] <- "Count"
gad_participant_ch_v2$Total<-rep(24, nrow(gad_participant_ch_v2))

gad_participant_ch_v2$percent_v2<-(gad_participant_ch_v2$Count/gad_participant_ch_v2$Total)*100

fig_4 <- ggplot(gad_participant_ch_v2, aes(x=categories, y=percent_v2,fill=categories)) + 
  geom_col()+labs(x="Anxiety categories", y = "Percentage")
fig_4<-fig_4+scale_fill_brewer(palette="Dark2")
plot(fig_4)

#Visit 4
participant_ch_v4<-metadata[(metadata$Visit=="T4"),]
anx_low_v4<-count(participant_ch_v4$Type=="Low")
anx_mild_v4<-count(participant_ch_v4$Type=="Mild")
anx_moderate_v4<-count(participant_ch_v4$Type=="Moderate")
anx_severe_v4<-count(participant_ch_v4$Type=="Severe")
gad_participant_ch_v4<-rbind(anx_low_v4[2,],anx_mild_v4[2,],anx_moderate_v4[2,])
gad_participant_ch_v4<-gad_participant_ch_v4[,2]
categories<-c("Low","Mild","Moderate")
gad_participant_ch_v4<-cbind.data.frame(levels,gad_participant_ch_v4)
names(gad_participant_ch_v4)[2] <- "Count"
gad_participant_ch_v4$Total<-rep(24, nrow(gad_participant_ch_v4))
gad_participant_ch_v4$percent_v4<-(gad_participant_ch_v4$Count/gad_participant_ch_v4$Total)*100

fig_5<- ggplot(gad_participant_ch_v4, aes(x=categories, y=percent_v4,fill=categories)) + 
  geom_col()+labs(x="Anxiety categories", y = "Percentage")
fig_5<-fig_5+scale_fill_brewer(palette="Dark2")
plot(fig_5)

new_row<-c("Severe",0,24,0.00)
gad_participant_ch_v4_new<-rbind.data.frame(gad_participant_ch_v4,new_row)

gad_participant_ch<-rep(c("Low","Mild","Moderate","Severe"),times=3)
gad_participant_ch<-cbind(gad_participant_ch, rep(c("T1","T2","T4"), each = 4))
colnames(gad_participant_ch)[1]<-"categories"
colnames(gad_participant_ch)[2]<-"Visit"

percent_1<-gad_participant_ch_v1$percent_v1
percent_2<-gad_participant_ch_v2$percent_v2
percent_4<-as.numeric(gad_participant_ch_v4_new$percent_v4)

percent<-c(percent_1,percent_2,percent_4)
percent<-round(percent)

gad_participant_ch<-cbind(gad_participant_ch,percent)
gad_participant_ch<-data.frame(gad_participant_ch)

data_4<- data.frame(Low = c(16, 9, 12),
                         Mild = c(5, 6, 7),
                         Moderate = c(1, 7, 5),
                         Severe = c(1, 2, 0)
)
rownames(data_4)<-c("T1","T2","T4")

data_4_v1<-data_4[c(1,2),]

model_9<-chisq.test(data_4_v1)
print(model_9)

data_4_v2<-data_4[c(2,3),]

model_10<-chisq.test(data_4_v2)
print(model_10)

data_4_v4<-data_4[c(1,3)]

model_11<-chisq.test(data_4_v4)
print(model_11)

# One box per treatment
fig_6<- ggplot(gad_participant_ch, aes(x=categories, y=as.numeric(percent), fill=categories)) +
  geom_bar(stat = "identity", position = "dodge")+facet_wrap(~Visit)+ labs(x="Anxiety categories", y = "Percentage")
fig_6<-plot(fig_6)+scale_fill_brewer(palette="Dark2")+theme_bw()
fig_6<-fig_6+theme(text = element_text(size = 20))
plot(fig_6)

############################################################################################################################

#Alpha diversity by GAD-7 scores  [###????ADD P-VALUES TO THE PLOT]

ps_v1<-phyloseq(tax_table(taxa_norm_v1),sample_data(metadata_v1),otu_table(otu_norm_v1, taxa_are_rows = FALSE))
ps_v1_alpha<-phyloseq(tax_table(taxa_norm_v1_alpha),sample_data(metadata_v1_alpha),otu_table(otu_norm_v1_alpha, taxa_are_rows = FALSE))


ps_v2<-phyloseq(tax_table(taxa_norm_v2),sample_data(metadata_v2),otu_table(otu_norm_v2, taxa_are_rows = FALSE))

ps_v4<-phyloseq(tax_table(taxa_norm_v4),sample_data(metadata_v4),otu_table(otu_norm_v4, taxa_are_rows = FALSE))
ps_v4_alpha<-phyloseq(tax_table(taxa_norm_v4_alpha),sample_data(metadata_v4_alpha),otu_table(otu_norm_v4_alpha, taxa_are_rows = FALSE))

set.seed(100)
model_12<-adonis2(shannon~SC5+BMI,data=metadata, permutations = 10000)
print(model_12)

alpha_plot_1<-plot_richness(ps_norm, x="SC5",measures =c("Shannon"))+ geom_point(color="Maroon1",size=4)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_1<-alpha_plot_1+labs(x="GAD-7 score", y = "Alpha Diversity")+ggtitle("All Visits")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_1<-alpha_plot_1+theme(text = element_text(size = 20))
plot(alpha_plot_1)


shannon_v1_alpha<-metadata_v1_alpha$shannon
SC5_v1_alpha<-metadata_v1_alpha$SC5  
BMI_v1_alpha<-metadata_v1_alpha$BMI    #min value is 21 and max is 41.1

set.seed(100)
model_13<-adonis2(shannon_v1_alpha~SC5_v1_alpha+BMI_v1_alpha,data=metadata_v1_alpha, permutations = 10000)
print(model_13)

alpha_plot_2<-plot_richness(ps_v1_alpha, x="SC5",measures =c("Shannon"))+ geom_point(color="Darkgreen",size=5)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_2<-alpha_plot_2+labs(x="GAD-7 score", y = "Alpha Diversity")+ggtitle("Visit 1")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_2<-alpha_plot_2+theme(text = element_text(size = 20))
plot(alpha_plot_2)


shannon_v2<-metadata_v2$shannon
SC5_v2<-metadata_v2$SC5
BMI_v2<-metadata_v2$BMI  #min value is 23 and max value is 44.3
set.seed(100)
model_14<-adonis2(shannon_v2~SC5_v2+BMI_v2,data=metadata_v2, permutations = 10000)
print(model_14)

alpha_plot_3<-plot_richness(ps_v2, x="SC5",measures =c("Shannon"))+ geom_point(color="chocolate",size=5)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_3<-alpha_plot_3+labs(x="GAD-7 score", y = "Alpha Diversity")+ggtitle("Visit 2")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_3<-alpha_plot_3+theme(text = element_text(size = 20))
plot(alpha_plot_3)


normalized_gestational_weeks_v4_alpha<-metadata_v4_alpha$Normalized_Gestationalweeks
shannon_v4_alpha<-metadata_v4_alpha$shannon
BMI_v4_alpha<-metadata_v4_alpha$BMI
SC5_v4_alpha<-metadata_v4_alpha$SC5

set.seed(100)


model_15<-adonis2(shannon_v4_alpha~SC5_v4_alpha+BMI_v4_alpha,data=metadata_v4_alpha, permutations = 10000)
print(model_15)


alpha_plot_4<-plot_richness(ps_v4_alpha, x="SC5_v4_alpha",measures =c("Shannon"))+ geom_point(color="lightslateblue",size=5)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_4<-alpha_plot_4+labs(x="GAD-7 score", y = "Alpha Diversity")+ggtitle("Visit 4")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_4<-alpha_plot_4+theme(text = element_text(size = 20))
plot(alpha_plot_4)



########With respect to normalized gestational weeks###################

set.seed(100)
model_16<-adonis2(shannon~normalized_gestational_weeks+BMI,data=metadata, permutations = 10000)
print(model_16)
alpha_plot_5<-plot_richness(ps_norm,x = "Normalized_Gestationalweeks",measures = c("Shannon"))+geom_point(color ="Maroon1",size = 4)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_5<-alpha_plot_5+labs(x="Normalized Gestational weeks", y = "Alpha Diversity")+ggtitle("All Visits")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_5<-alpha_plot_5+theme(text = element_text(size = 20))
plot(alpha_plot_5)

normalized_gestational_weeks_v1_alpha<-metadata_v1_alpha$Normalized_Gestationalweeks

set.seed(100)
model_17<-adonis2(shannon_v1_alpha~normalized_gestational_weeks_v1_alpha+BMI_v1_alpha,data=metadata_v1_alpha, permutations = 10000)
print(model_17)
alpha_plot_6<-plot_richness(ps_v1_alpha, x="normalized_gestational_weeks_v1_alpha",measures =c("Shannon"))+ geom_point(color="Darkgreen",size=5)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_6<-alpha_plot_6+labs(x="Normalized Gestationalweeks", y = "Alpha Diversity")+ggtitle("Visit 1")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_6<-alpha_plot_6+theme(text = element_text(size = 20))
plot(alpha_plot_6)

set.seed(100)
normalized_gestational_weeks_v2<-metadata_v2$Normalized_Gestationalweeks
model_18<-adonis2(shannon_v2~normalized_gestational_weeks_v2+BMI_v2,data=metadata_v2, permutations = 10000)
print(model_18)
alpha_plot_7<-plot_richness(ps_v2, x="normalized_gestational_weeks_v2",measures =c("Shannon"))+ geom_point(color="chocolate",size=5)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_7<-alpha_plot_7+labs(x="Normalized Gestationalweeks", y = "Alpha Diversity")+ggtitle("Visit 2")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_7<-alpha_plot_7+theme(text = element_text(size = 20))
plot(alpha_plot_7)

set.seed(100)


model_19<-adonis2(shannon_v4_alpha~normalized_gestational_weeks_v4_alpha+BMI_v4_alpha,data=metadata_v4_alpha, permutations = 10000)
print(model_19)

alpha_plot_8<-plot_richness(ps_v4_alpha, x="normalized_gestational_weeks_v4_alpha",measures =c("Shannon"))+ geom_point(color="lightslateblue",size=5)+geom_smooth(method = "lm",col="black",size=2)
alpha_plot_8<-alpha_plot_8+labs(x="Gestational weeks", y = "Alpha Diversity")+ggtitle("Visit 4")+theme_bw()+theme(aspect.ratio=1)
alpha_plot_8<-alpha_plot_8+theme(text = element_text(size = 20))
plot(alpha_plot_8)

######Perform PCoA using weighted and uweighted Uni-frac distances   ########Beta Diversity##############
set.seed(100)

bray_curtis<-vegdist((otu_table(ps_norm)),method="bray")  #####Calculating distance 
bray_curtis<-as.matrix(bray_curtis)
#rownames(bray_curtis)<-paste("Seq",1:nrow(bray_curtis), sep="_")
#colnames(bray_curtis)<-paste("Seq",1:ncol(bray_curtis), sep ="_")

bray_curtis_melt<-melt(bray_curtis)        ############sample 1- 2 --dist
colnames(bray_curtis_melt)<-c("Sample1","Sample2", "BrayCurtisDistance")
bray_curtis_melt_Sample1<-merge(bray_curtis_melt,sample_data(ps_norm)[,c("SC5")], by.x="Sample1",
      by.y=0, all.x=TRUE)

bray_curtis_melt_both<-bray_curtis_melt_Sample1[-which(bray_curtis_melt_Sample1$Sample1==bray_curtis_melt_Sample1$Sample2),]

average_bc_anx<-data.frame(cbind(unique(bray_curtis_melt_both$SC5), 
                                 sapply(unique(bray_curtis_melt_both$SC5), function(x){
                                   average=mean(bray_curtis_melt_both[which(bray_curtis_melt_both$SC5==x),"BrayCurtisDistance"])})))
colnames(average_bc_anx)<-c("SC5", "AverageBrayCurtis")

AverageBrayCurtis_GAD<-average_bc_anx$AverageBrayCurtis
gadscores<-average_bc_anx$SC5

set.seed(100)
model_20<-adonis2(AverageBrayCurtis_GAD~gadscores,data=average_bc_anx,permutations = 10000)
print(model_20)

beta_plot_1<-ggplot(average_bc_anx, aes(x=SC5, y=AverageBrayCurtis))+geom_point(color="Maroon1",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("All Visits")
beta_plot_1<-beta_plot_1+labs(x="GAD-7 score", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)+theme(text = element_text(size = 20))
plot(beta_plot_1)


#####Beta diversity by Visit1######################

bray_curtis_v1<-vegdist((otu_table(ps_v1)),method="bray")  #####Calculating distance 
bray_curtis_v1<-as.matrix(bray_curtis_v1)
#rownames(bray_curtis)<-paste("Seq",1:nrow(bray_curtis), sep="_")
#colnames(bray_curtis)<-paste("Seq",1:ncol(bray_curtis), sep ="_")

bray_curtis_melt_v1<-melt(bray_curtis_v1)        ############sample 1- 2 --dist
colnames(bray_curtis_melt_v1)<-c("Sample1","Sample2", "BrayCurtisDistance")
bray_curtis_melt_v1_Sample1<-merge(bray_curtis_melt_v1,sample_data(ps_v1)[,c("SC5")], by.x="Sample1",
                                by.y=0, all.x=TRUE)
#bray_curtis_melt_both<-merge(bray_curtis_melt_Sample1,sample_data(ps2)[,c("SC5")], by.x="Sample1",
#    by.y=0, all.x=TRUE)

bray_curtis_melt_both_v1<-bray_curtis_melt_v1_Sample1[-which(bray_curtis_melt_v1_Sample1$Sample1==bray_curtis_melt_v1_Sample1$Sample2),]

average_bc_anx_v1<-data.frame(cbind(unique(bray_curtis_melt_both_v1$SC5), 
                                 sapply(unique(bray_curtis_melt_both_v1$SC5), function(x){
                                   average=mean(bray_curtis_melt_both_v1[which(bray_curtis_melt_both_v1$SC5==x),"BrayCurtisDistance"])})))
colnames(average_bc_anx_v1)<-c("SC5", "AverageBrayCurtis")


AverageBrayCurtis_GAD_v1<-average_bc_anx_v1$AverageBrayCurtis
gadscores_v1<-average_bc_anx_v1$SC5
set.seed(100)
model_21<-adonis2(AverageBrayCurtis_GAD_v1~gadscores_v1,data=average_bc_anx_v1,permutations = 10000)
print(model_21)

beta_plot_2<-ggplot(average_bc_anx_v1, aes(x=SC5, y=AverageBrayCurtis_v1))+geom_point(color="Darkgreen",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("Visit 1")
beta_plot_2<-beta_plot_2+labs(x="GAD-7 score", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)
beta_plot_2<-beta_plot_2+theme(text = element_text(size = 20))
plot(beta_plot_2)

#######Beta diversity by Visit 2####################

bray_curtis_v2<-vegdist((otu_table(ps_v2)),method="bray")  #####Calculating distance 
bray_curtis_v2<-as.matrix(bray_curtis_v2)

bray_curtis_melt_v2<-melt(bray_curtis_v2)        ############sample 1- 2 --dist
colnames(bray_curtis_melt_v2)<-c("Sample1","Sample2", "BrayCurtisDistance")
bray_curtis_melt_v2_Sample1<-merge(bray_curtis_melt_v2,sample_data(ps_v2)[,c("SC5")], by.x="Sample1",
                                   by.y=0, all.x=TRUE)

bray_curtis_melt_both_v2<-bray_curtis_melt_v2_Sample1[-which(bray_curtis_melt_v2_Sample1$Sample1==bray_curtis_melt_v2_Sample1$Sample2),]

average_bc_anx_v2<-data.frame(cbind(unique(bray_curtis_melt_both_v2$SC5), 
                                    sapply(unique(bray_curtis_melt_both_v2$SC5), function(x){
                                      average=mean(bray_curtis_melt_both_v2[which(bray_curtis_melt_both_v2$SC5==x),"BrayCurtisDistance"])})))
colnames(average_bc_anx_v2)<-c("SC5", "AverageBrayCurtis")

AverageBrayCurtis_GAD_v2<-average_bc_anx_v2$AverageBrayCurtis
gadscores_v2<-average_bc_anx_v2$SC5
set.seed(100)
model_22<-adonis2(AverageBrayCurtis_GAD_v2~gadscores_v2,data=average_bc_anx_v2,permutations = 10000)
print(model_22)

beta_plot_3<-ggplot(average_bc_anx_v2, aes(x=SC5, y=AverageBrayCurtis_v2))+geom_point(color="chocolate",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("Visit 2")
beta_plot_3<-beta_plot_3+labs(x="GAD-7 score", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)+theme(text = element_text(size = 20))
plot(beta_plot_3)

########Beta diversity by Visit 4################

bray_curtis_v4<-vegdist((otu_table(ps_v4)),method="bray")  #####Calculating distance 
bray_curtis_v4<-as.matrix(bray_curtis_v4)
#rownames(bray_curtis)<-paste("Seq",1:nrow(bray_curtis), sep="_")
#colnames(bray_curtis)<-paste("Seq",1:ncol(bray_curtis), sep ="_")

bray_curtis_melt_v4<-melt(bray_curtis_v4)        ############sample 1- 2 --dist
colnames(bray_curtis_melt_v4)<-c("Sample1","Sample2", "BrayCurtisDistance")
bray_curtis_melt_v4_Sample1<-merge(bray_curtis_melt_v4,sample_data(ps_v4)[,c("SC5")], by.x="Sample1",
                                   by.y=0, all.x=TRUE)
#bray_curtis_melt_both<-merge(bray_curtis_melt_Sample1,sample_data(ps2)[,c("SC5")], by.x="Sample1",
#    by.y=0, all.x=TRUE)

bray_curtis_melt_both_v4<-bray_curtis_melt_v4_Sample1[-which(bray_curtis_melt_v4_Sample1$Sample1==bray_curtis_melt_v4_Sample1$Sample2),]

average_bc_anx_v4<-data.frame(cbind(unique(bray_curtis_melt_both_v4$SC5), 
                                    sapply(unique(bray_curtis_melt_both_v4$SC5), function(x){
                                      average=mean(bray_curtis_melt_both_v4[which(bray_curtis_melt_both_v4$SC5==x),"BrayCurtisDistance"])})))
colnames(average_bc_anx_v4)<-c("SC5", "AverageBrayCurtis")

AverageBrayCurtis_GAD_v4<-average_bc_anx_v4$AverageBrayCurtis
gadscores_v4<-average_bc_anx_v4$SC5
set.seed(100)
model_23<-adonis2(AverageBrayCurtis_GAD_v4~gadscores_v4,data=average_bc_anx_v4,permutations = 10000)
print(model_23)

beta_plot_4<-ggplot(average_bc_anx_v4, aes(x=SC5, y=AverageBrayCurtis_v4))+geom_point(color="lightslateblue",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("Visit 4")
beta_plot_4<-beta_plot_4+labs(x="GAD-7 score", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)+theme(text = element_text(size = 20))
plot(beta_plot_4)

#####################################################################################################

#By gestational weeks

bray_curtis_melt_weeks<-merge(bray_curtis_melt,sample_data(ps_norm)[,c("Normalized_Gestationalweeks")], by.x="Sample1",
                              by.y=0, all.x=TRUE)
bray_curtis_melt_weeks<-bray_curtis_melt_weeks[-which(bray_curtis_melt_weeks$Sample1==bray_curtis_melt_weeks$Sample2),]

average_bc_weeks<-data.frame(cbind(unique(bray_curtis_melt_weeks$Normalized_Gestationalweeks), 
                                   sapply(unique(bray_curtis_melt_weeks$Normalized_Gestationalweeks), function(x){
                                     average=mean(bray_curtis_melt_weeks[which(bray_curtis_melt_weeks$Normalized_Gestationalweeks==x),"BrayCurtisDistance"])})))
colnames(average_bc_weeks)<-c("Normalized_Gestationalweeks", "AverageBrayCurtis")

AverageBrayCurtis_weeks<-average_bc_weeks$AverageBrayCurtis
gestational_weeks_all<-average_bc_weeks$Normalized_Gestationalweeks
set.seed(100)
model_24<-adonis2(AverageBrayCurtis_weeks~gestational_weeks_all,data=average_bc_weeks,permutations = 10000)
print(model_24)

beta_plot_5<-ggplot(average_bc_weeks, aes(x=Normalized_Gestationalweeks, y=AverageBrayCurtis))+geom_point(color="Maroon1",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("All Visits")
beta_plot_5<-beta_plot_5+labs(x="Normalized Gestationalweeks", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)+theme(text = element_text(size = 20))
plot(beta_plot_5)

bray_curtis_melt_v1_weeks<-merge(bray_curtis_melt_v1,sample_data(ps_v1)[,c("Normalized_Gestationalweeks")], by.x="Sample1",
                                 by.y=0, all.x=TRUE)
bray_curtis_melt_v1_weeks<-bray_curtis_melt_v1_weeks[-which(bray_curtis_melt_v1_weeks$Sample1==bray_curtis_melt_v1_weeks$Sample2),]
average_bc_weeks_v1<-data.frame(cbind(unique(bray_curtis_melt_v1_weeks$Normalized_Gestationalweeks), 
                                      sapply(unique(bray_curtis_melt_v1_weeks$Normalized_Gestationalweeks), function(x){
                                        average=mean(bray_curtis_melt_v1_weeks[which(bray_curtis_melt_v1_weeks$Normalized_Gestationalweeks==x),"BrayCurtisDistance"])})))
colnames(average_bc_weeks_v1)<-c("Normalized_Gestationalweeks", "AverageBrayCurtis")

AverageBrayCurtis_weeks_v1<-average_bc_weeks_v1$AverageBrayCurtis
gestational_weeks_v1<-average_bc_weeks_v1$Normalized_Gestationalweeks


set.seed(100)
model_25<-adonis2(AverageBrayCurtis_weeks_v1~gestational_weeks_v1,data=average_bc_weeks_v1,permutations = 10000)
print(model_25)

beta_plot_6<-ggplot(average_bc_weeks_v1, aes(x=Normalized_Gestationalweeks, y=AverageBrayCurtis))+geom_point(color="Darkgreen",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("Visit 1")
beta_plot_6<-beta_plot_6+labs(x="Normalized Gestationalweeks", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)+theme(text = element_text(size = 20))
plot(beta_plot_6)

bray_curtis_melt_v2_weeks<-merge(bray_curtis_melt_v2,sample_data(ps_v2)[,c("Normalized_Gestationalweeks")], by.x="Sample1",
                                 by.y=0, all.x=TRUE)

bray_curtis_melt_v2_weeks<-bray_curtis_melt_v2_weeks[-which(bray_curtis_melt_v2_weeks$Sample1==bray_curtis_melt_v2_weeks$Sample2),]

average_bc_weeks_v2<-data.frame(cbind(unique(bray_curtis_melt_v2_weeks$Normalized_Gestationalweeks), 
                                      sapply(unique(bray_curtis_melt_v2_weeks$Normalized_Gestationalweeks), function(x){
                                        average=mean(bray_curtis_melt_v2_weeks[which(bray_curtis_melt_v2_weeks$Normalized_Gestationalweeks==x),"BrayCurtisDistance"])})))
colnames(average_bc_weeks_v2)<-c("Normalized_Gestationalweeks", "AverageBrayCurtis")

AverageBrayCurtis_weeks_v2<-average_bc_weeks_v2$AverageBrayCurtis
gestational_weeks_v2<-average_bc_weeks_v2$Normalized_Gestationalweeks


set.seed(100)
model_26<-adonis2(AverageBrayCurtis_weeks_v2~gestational_weeks_v2,data=average_bc_weeks_v2,permutations = 10000)
print(model_26)
beta_plot_7<-ggplot(average_bc_weeks_v2, aes(x=Normalized_Gestationalweeks, y=AverageBrayCurtis))+geom_point(color="chocolate",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("Visit 2")
beta_plot_7<-beta_plot_7+labs(x="Normalized Gestationalweeks", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)+theme(text = element_text(size = 20))
plot(beta_plot_7)




bray_curtis_melt_v4_weeks<-merge(bray_curtis_melt_v4,sample_data(ps_v4)[,c("Normalized_Gestationalweeks")], by.x="Sample1",
                                 by.y=0, all.x=TRUE)

bray_curtis_melt_v4_weeks<-bray_curtis_melt_v4_weeks[-which(bray_curtis_melt_v4_weeks$Sample1==bray_curtis_melt_v4_weeks$Sample2),]


average_bc_weeks_v4<-data.frame(cbind(unique(bray_curtis_melt_v4_weeks$Normalized_Gestationalweeks), 
                                      sapply(unique(bray_curtis_melt_v4_weeks$Normalized_Gestationalweeks), function(x){
                                        average=mean(bray_curtis_melt_v4_weeks[which(bray_curtis_melt_v4_weeks$Normalized_Gestationalweeks==x),"BrayCurtisDistance"])})))
colnames(average_bc_weeks_v4)<-c("Normalized_Gestationalweeks", "AverageBrayCurtis")


AverageBrayCurtis_weeks_v4<-average_bc_weeks_v4$AverageBrayCurtis
gestational_weeks_v4<-average_bc_weeks_v4$Normalized_Gestationalweeks

set.seed(100)
model_27<-adonis2(AverageBrayCurtis_weeks_v4~gestational_weeks_v4,data=average_bc_weeks_v4,permutations = 10000)
print(model_27)

beta_plot_8<-ggplot(average_bc_weeks_v4, aes(x=Normalized_Gestationalweeks, y=AverageBrayCurtis))+geom_point(color="lightslateblue",size=5)+geom_smooth(method="lm",color="black",size=2)+ggtitle("Visit 4")
beta_plot_8<-beta_plot_8+labs(x="Gestational weeks", y = "AverageBrayCurtis")+theme_bw()+theme(aspect.ratio=1)+theme(text = element_text(size = 20))
plot(beta_plot_8)

############################################################################################################################################
#metagenomeSeq is designed to determine features (be it Operational Taxanomic Unit (OTU), species, etc.) that are differentially abundant between two or more groups of multiple samples. metagenomeSeq is designed to address the effects of both normalization and under-sampling of microbial communities on disease association detection and the testing of feature correlations.

#All visits together

set.seed(100)
featureData =data.frame(tax_table(ps2))
rownames(featureData)<-paste("Seq",1:nrow(featureData), sep="_")
matrixData<-matrix(otu_table(ps2),ncol=ncol(otu_table(ps2)))
rownames(matrixData)<-rownames(otu_table(ps2))
colnames(matrixData)<-colnames(otu_table(ps2))
matrixData<-t(matrixData)
rownames(matrixData)<-paste("Seq",1:nrow(matrixData), sep="_")

metadata<-sample_data(ps2)[match(colnames(matrixData),rownames(sample_data(ps2))), ]
metadata<-metadata[match(colnames(matrixData),rownames(metadata)),]

otus_metagenomeSeq<-newMRexperiment(matrixData, 
                                    phenoData = AnnotatedDataFrame(metadata), 
                                    featureData = AnnotatedDataFrame(featureData))

otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq,present=1, depth = 1)
p = cumNormStatFast(otus_metagenomeSeq_filter)    #percentile by which to normalize counts
anxietyData = cumNorm(otus_metagenomeSeq, p = p)
anxiety_sample_data<-pData(anxietyData)
#anxietyType<-anxiety_sample_data$Type
#visitStatus<-anxiety_sample_data$Visit
anxietyscores<-anxiety_sample_data$SC5
#anxietyResult<-anxiety_sample_data$Result
gestational_weeks_meta<-anxiety_sample_data$Normalized_Gestationalweeks
#age_meta<-anxiety_sample_data$Age
#education_meta<-anxiety_sample_data$Education
bmi<-anxiety_sample_data$BMI

mod<-model.matrix(~anxietyscores+gestational_weeks_meta+bmi)
settings<-zigControl(maxit = 10, verbose = TRUE)
fit<-fitZig(obj = anxietyData, mod = mod, control = settings, useCSSoffset = TRUE)

z1<-median(calculateEffectiveSamples(fit))

mean_eff_size = round(z1)    #####################mean effective size is 15##############
otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq, present = mean_eff_size, 
                                       depth = 1) 
p = cumNormStatFast(otus_metagenomeSeq_filter)    #percentile by which to normalize counts
anxietyData = cumNorm(otus_metagenomeSeq_filter, p = p)

#Calculating count matrix
anxietyData_norm<-MRcounts(anxietyData, norm=TRUE, log=FALSE)
anxiety_sample_data<-pData(anxietyData)
#visitStatus<-anxiety_sample_data$Visit
#anxietyType<-anxiety_sample_data$Type
#anxietyResult<-anxiety_sample_data$Result
anxietyscores<-anxiety_sample_data$SC5
#anxiety_study<-anxiety_sample_data$Study
gestational_weeks_meta<-anxiety_sample_data$Normalized_Gestationalweeks
#age_meta<-anxiety_sample_data$Age
#education_meta<-anxiety_sample_data$Education
#gestational_weeks<-as.factor(anxiety_sample_data$GestationalWeeks_sinceConception)
bmi<-anxiety_sample_data$BMI

#Define the normalisation factor
norm.factor <- normFactors(anxietyData)
norm.factor <- log2(norm.factor/median(norm.factor)+1)       #min value is 0.028 and max is 1.53

mod<-model.matrix(~anxietyscores+gestational_weeks_meta+bmi+norm.factor)

settings<-zigControl(maxit = 10, verbose = TRUE)
fit<-fitZig(obj = anxietyData, mod = mod, control = settings, useCSSoffset = FALSE)

#coef<-MRcoefs(fit1,by=2)
#head(MRcoefs(fit1))

zigFit<-slot(fit,"fit")
finalMod<-slot(fit,"fit")$design

fit<-eBayes(zigFit, trend = TRUE, robust = TRUE)
topTable(fit)

pvalues<-apply(fit$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit$coefficients,pvalues)
otus_mG<-data.frame(otus_mG)

otus_mG_filter_anxiety<-otus_mG[(otus_mG$anxietyscores.1<0.1),]
otus_mG_filter_anxiety<-otus_mG_filter_anxiety[(otus_mG_filter_anxiety$gestational_weeks_meta.1>0.1),]
otus_mG_filter_anxiety<-otus_mG_filter_anxiety[(otus_mG_filter_anxiety$bmi.1>0.1),]
otus_mG_filter_anxiety<-otus_mG_filter_anxiety[(otus_mG_filter_anxiety$norm.factor.1>0.1),]
otus_mG_filter_anxiety<-cbind(otus_mG_filter_anxiety, Type = rep("anxiety only", nrow(otus_mG_filter_anxiety)))

otus_mG_filter_anxiety_weeks<-otus_mG[(otus_mG$anxietyscores.1<0.1),]
otus_mG_filter_anxiety_weeks<-otus_mG_filter_anxiety_weeks[(otus_mG_filter_anxiety_weeks$gestational_weeks_meta.1<0.1),]
otus_mG_filter_anxiety_weeks<-cbind(otus_mG_filter_anxiety_weeks, Type = rep("anxiety and weeks", nrow(otus_mG_filter_anxiety_weeks)))

otus_mG_filter_anxiety_bmi<-otus_mG[(otus_mG$anxietyscores.1<0.1),]
otus_mG_filter_anxiety_bmi<-otus_mG_filter_anxiety_bmi[(otus_mG_filter_anxiety_bmi$bmi.1<0.1),]
otus_mG_filter_anxiety_bmi<-cbind(otus_mG_filter_anxiety_bmi, Type = rep("anxiety and bmi", nrow(otus_mG_filter_anxiety_bmi)))

#otus_all <- merge(otus_mG_filter_anxiety,otus_mG_filter_anxiety_weeks,otus_mG_filter_anxiety_bmi, by = "anxietyscores", all = TRUE)

otus_all<-rbind.data.frame(otus_mG_filter_anxiety,otus_mG_filter_anxiety_weeks,otus_mG_filter_anxiety_bmi)

#otus_all<-rownames(otus_all)
h<-rownames(otus_all)

#m<-data.frame(matrixData)
#n1<-data.frame(featureData)

meta_taxa<-data.frame(featureData[c(h),])
#meta_otu<-m[c(h),]
#data<-p[c(h),]

#otu_normalized_all<-anxietyData_norm[c(otu_all),]
#rownames(otus_all) <- paste(rownames(otus_all), " ", "(", meta_taxa$Genus, ")", sep = "")


t<-cbind.data.frame(meta_taxa$Genus,otus_all$anxietyscores,otus_all$Type)
#rownames(t)<-rownames(meta_otu)
#t<-na.omit(t)

colnames(t)<-c("Genus","Coefficient","Type")
sign<-sapply(as.numeric(t[,"Coefficient"]),sign)
t<-cbind.data.frame(t,sign)
t$Visit<-rep("All", nrow(t))

#rownames(otu_normalized_all)<-t$Genus

ggplot(t,aes(y=Coefficient, x=Genus, col= as.factor(sign)))+xlab(label="Genus")+ylab(label=paste("Coefficient"))+
  coord_flip()+
  geom_point(aes(shape=as.factor(Type)),size = 5)+scale_colour_manual(values = c("-1" = "Maroon1", "1" = "lightslateblue"))+ggtitle("All Visits")+
  geom_abline(slope =0)+guides(col=guide_legend("Sign"),
                               shape=guide_legend("Type"))+
  theme_bw()+theme(aspect.ratio = 1)+
  theme(text = element_text(size = 20))


write.csv(t, file="Significant_allvisits.csv")
#For Visit 1############################################

set.seed(100)
matrixData<-t(matrixData)
matrixData_v1<-matrixData[c(colnames_v1),]
matrixData_v1<-matrixData_v1[,colSums(matrixData_v1)>0]
matrixData_v1<-matrixData_v1[rowSums(matrixData_v1)>0,]
matrixData_v1<-t(matrixData_v1)
featureData_v1<-featureData[match(rownames(matrixData_v1),rownames(featureData)),]
featureData_v1<-data.frame(featureData_v1)

otus_metagenomeSeq_v1<-newMRexperiment(matrixData_v1, 
                                    phenoData = AnnotatedDataFrame(metadata_v1), 
                                    featureData = AnnotatedDataFrame(featureData_v1))

otus_metagenomeSeq_filter_v1 = filterData(otus_metagenomeSeq_v1,present=1, depth = 1)
p_v1= cumNormStatFast(otus_metagenomeSeq_filter_v1) 

anxietyData_v1 = cumNorm(otus_metagenomeSeq_v1, p = p_v1)
anxiety_sample_data_v1<-pData(anxietyData_v1)

#anxietyResult_v1<-anxiety_sample_data_v1$Result
bmi_v1<-anxiety_sample_data_v1$BMI
#age_v1<-anxiety_sample_data_v1$Age
#gestational_weeks_v1<-anxiety_sample_data_v1$Normalized_Gestationalweeks
anxietyscores_v1<-anxiety_sample_data_v1$SC5

mod_v1<-model.matrix(~anxietyscores_v1+bmi_v1)
settings<-zigControl(maxit = 10, verbose = TRUE)
fit_v1<-fitZig(obj = anxietyData_v1, mod = mod_v1,control = settings, useCSSoffset = FALSE)

mean_eff_size_v1<-median(calculateEffectiveSamples(fit_v1))
mean_eff_size_v1 = round(mean_eff_size_v1)    #####################mean effective size is 7##############

otus_metagenomeSeq_filter_v1 = filterData(otus_metagenomeSeq_v1, present = mean_eff_size_v1, 
                                       depth = 1) 
p_v1 = cumNormStatFast(otus_metagenomeSeq_filter_v1)    #percentile by which to normalize counts

anxietyData_v1 = cumNorm(otus_metagenomeSeq_filter_v1, p = p_v1)

anxietyData_v1_norm<-MRcounts(anxietyData_v1, norm=TRUE, log=FALSE)
anxiety_sample_data_v1<-pData(anxietyData_v1)

#anxietyResult_v1<-anxiety_sample_data_v1$Result
anxietyscores_v1<-anxiety_sample_data_v1$SC5
bmi_v1<-anxiety_sample_data_v1$BMI
#age_v1<-anxiety_sample_data_v1$Age
#gestational_weeks_v1<-anxiety_sample_data_v1$Normalized_Gestationalweeks
# Define the normalisation factor
norm.factor_v1 <- normFactors(anxietyData_v1)                           ### min value is 0.03    and max value is 1.33
norm.factor_v1 <- log2(norm.factor_v1/median(norm.factor_v1) + 1)

mod_v1<-model.matrix(~anxietyscores_v1+bmi_v1+norm.factor_v1)

settings<-zigControl(maxit = 10, verbose = TRUE)
fit_v1<-fitZig(obj = anxietyData_v1, mod = mod_v1, control = settings, useCSSoffset = FALSE)

zigFit_v1<-slot(fit_v1,"fit")
finalMod_v1<-slot(fit_v1,"fit")$design

fit_v1<-eBayes(zigFit_v1, trend = TRUE, robust = TRUE)
topTable(fit_v1)

#coef_v1<-MRcoefs(fit_v1,by=2)

pvalues_v1<-apply(fit_v1$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG_v1<-cbind(fit_v1$coefficients,pvalues_v1)
otus_mG_v1<-data.frame(otus_mG_v1)

otus_mG_filter_anxiety_v1<-otus_mG_v1[(otus_mG_v1$anxietyscores_v1.1<0.1),]
otus_mG_filter_anxiety_v1<-otus_mG_filter_anxiety_v1[(otus_mG_filter_anxiety_v1$bmi_v1.1>0.1),]
otus_mG_filter_anxiety_v1<-otus_mG_filter_anxiety_v1[(otus_mG_filter_anxiety_v1$norm.factor_v1.1>0.1),]
otus_mG_filter_anxiety_v1<-cbind(otus_mG_filter_anxiety_v1, Type = rep("anxiety only", nrow(otus_mG_filter_anxiety_v1)))

otus_mG_filter_anxiety_bmi_v1<-otus_mG_v1[(otus_mG_v1$anxietyscores_v1.1<0.1),]
otus_mG_filter_anxiety_bmi_v1<-otus_mG_filter_anxiety_bmi_v1[(otus_mG_filter_anxiety_bmi_v1$bmi_v1.1<0.1),]
otus_mG_filter_anxiety_bmi_v1<-cbind(otus_mG_filter_anxiety_bmi_v1, Type = rep("anxiety and bmi", nrow(otus_mG_filter_anxiety_bmi_v1)))

otus_v1<-rbind(otus_mG_filter_anxiety_v1,otus_mG_filter_anxiety_bmi_v1)

h_v1<-rownames(otus_v1)

meta_taxa_v1<-data.frame(featureData[c(h_v1),])

t_v1<-cbind.data.frame(meta_taxa_v1$Genus,otus_v1$anxietyscores_v1,otus_v1$Type)

colnames(t_v1)<-c("Genus","Coefficient","Type")
sign_v1<-sapply(as.numeric(t_v1[,"Coefficient"]),sign)
t_v1<-cbind.data.frame(t_v1,sign_v1)
t_v1$Visit<-rep("Visit 1", nrow(t_v1))



ggplot(t_v1,aes(y=Coefficient, x=Genus, col= as.factor(sign_v1)))+xlab(label="Genus")+ylab(label=paste("Coefficient"))+
  coord_flip()+
  geom_point(aes(shape=as.factor(Type)),size = 5)+scale_colour_manual(values = c("-1" = "chocolate", "1" = "Darkgreen"))+ggtitle("Visit 1")+
  geom_abline(slope =0)+guides(col=guide_legend("Sign"),
                               shape=guide_legend("Type"))+
  theme_bw()+theme(aspect.ratio = 1)+
  theme(text = element_text(size = 20))


write.csv(t_v1, file="Significant_visit1.csv")


#For visit 2#############################################

set.seed(100)
#matrixData<-t(matrixData)
matrixData_v2<-matrixData[c(colnames_v2),]
matrixData_v2<-matrixData_v2[,colSums(matrixData_v2)>0]
matrixData_v2<-matrixData_v2[rowSums(matrixData_v2)>0,]
matrixData_v2<-t(matrixData_v2)
featureData_v2<-featureData[match(rownames(matrixData_v2),rownames(featureData)),]
featureData_v2<-data.frame(featureData_v2)

otus_metagenomeSeq_v2<-newMRexperiment(matrixData_v2, 
                                       phenoData = AnnotatedDataFrame(metadata_v2), 
                                       featureData = AnnotatedDataFrame(featureData_v2))

otus_metagenomeSeq_filter_v2 = filterData(otus_metagenomeSeq_v2,present=1, depth = 1)
p_v2= cumNormStatFast(otus_metagenomeSeq_filter_v2) 

anxietyData_v2 = cumNorm(otus_metagenomeSeq_v2, p = p_v2)

anxietyData_v2_norm<-MRcounts(anxietyData_v2, norm=TRUE, log=FALSE)

anxiety_sample_data_v2<-pData(anxietyData_v2)

anxietyResult_v2<-anxiety_sample_data_v2$Result
bmi_v2<-anxiety_sample_data_v2$BMI
#gestational_weeks_v2<-anxiety_sample_data_v2$Normalized_Gestationalweeks
anxietyscores_v2<-anxiety_sample_data_v2$SC5

mod_v2<-model.matrix(~anxietyscores_v2+bmi_v2)
settings<-zigControl(maxit = 10, verbose = TRUE)
fit_v2<-fitZig(obj = anxietyData_v2, mod = mod_v2,control = settings, useCSSoffset = FALSE)

mean_eff_size_v2<-median(calculateEffectiveSamples(fit_v2))
mean_eff_size_v2 = round(mean_eff_size_v2)    #####################mean effective size is 7##############

otus_metagenomeSeq_filter_v2 = filterData(otus_metagenomeSeq_v2, present = mean_eff_size_v2, 
                                          depth = 1) 
p_v2 = cumNormStatFast(otus_metagenomeSeq_filter_v2)    #percentile by which to normalize counts

anxietyData_v2 = cumNorm(otus_metagenomeSeq_filter_v2, p = p_v2)
anxiety_sample_data_v2<-pData(anxietyData_v2)

#anxietyResult_v2<-anxiety_sample_data_v2$Result
bmi_v2<-anxiety_sample_data_v2$BMI
#gestational_weeks_v2<-anxiety_sample_data_v2$Normalized_Gestationalweeks
anxietyscores_v2<-anxiety_sample_data_v2$SC5

norm.factor_v2<- normFactors(anxietyData_v2)                           ### min value is 0.03    and max value is 1.47
norm.factor_v2<- log2(norm.factor_v2/median(norm.factor_v2) + 1)

mod_v2<-model.matrix(~anxietyscores_v2+bmi_v2+norm.factor_v2)

settings<-zigControl(maxit = 10, verbose = TRUE)
fit_v2<-fitZig(obj = anxietyData_v2, mod = mod_v2, control = settings, useCSSoffset = FALSE)

zigFit_v2<-slot(fit_v2,"fit")
finalMod_v2<-slot(fit_v2,"fit")$design

fit_v2<-eBayes(zigFit_v2, trend = TRUE, robust = TRUE)
topTable(fit_v2)

pvalues_v2<-apply(fit_v2$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG_v2<-cbind(fit_v2$coefficients,pvalues_v2)
otus_mG_v2<-data.frame(otus_mG_v2)

#coef_v2<-MRcoefs(fit_v2,by=2)

otus_mG_filter_anxiety_v2<-otus_mG_v2[(otus_mG_v2$anxietyscores_v2.1<0.1),]
otus_mG_filter_anxiety_v2<-otus_mG_filter_anxiety_v2[(otus_mG_filter_anxiety_v2$bmi_v2.1>0.1),]
otus_mG_filter_anxiety_v2<-otus_mG_filter_anxiety_v2[(otus_mG_filter_anxiety_v2$norm.factor_v2.1>0.1),]
otus_mG_filter_anxiety_v2<-cbind(otus_mG_filter_anxiety_v2, Type = rep("anxiety only", nrow(otus_mG_filter_anxiety_v2)))

otus_mG_filter_anxiety_bmi_v2<-otus_mG_v2[(otus_mG_v2$anxietyscores_v2.1<0.1),]
otus_mG_filter_anxiety_bmi_v2<-otus_mG_filter_anxiety_bmi_v2[(otus_mG_filter_anxiety_bmi_v2$bmi_v2.1<0.1),]
otus_mG_filter_anxiety_bmi_v2<-cbind(otus_mG_filter_anxiety_bmi_v2, Type = rep("anxiety and bmi", nrow(otus_mG_filter_anxiety_bmi_v2)))

otus_v2<-rbind(otus_mG_filter_anxiety_v2,otus_mG_filter_anxiety_bmi_v2)

h_v2<-rownames(otus_v2)

meta_taxa_v2<-data.frame(featureData[c(h_v2),])

t_v2<-cbind.data.frame(meta_taxa_v2$Genus,otus_v2$anxietyscores_v2,otus_v2$Type)
t_v2<-na.omit(t_v2)

colnames(t_v2)<-c("Genus","Coefficient","Type")
sign_v2<-sapply(as.numeric(t_v2[,"Coefficient"]),sign)
t_v2<-cbind.data.frame(t_v2,sign_v2)
t_v2$Visit<-rep("Visit 2", nrow(t_v2))

ggplot(t_v2,aes(y=Coefficient, x=Genus, col= as.factor(sign_v2)))+xlab(label="Genus")+ylab(label=paste("Coefficient"))+
  coord_flip()+
  geom_point(aes(shape=as.factor(Type)),size = 5)+scale_colour_manual(values = c("-1" = "lightslateblue", "1" = "Maroon1"))+ggtitle("Visit 2")+
  geom_abline(slope =0)+guides(col=guide_legend("Sign"),
                               shape=guide_legend("Type"))+
  theme_bw()+theme(aspect.ratio = 1)+
  theme(text = element_text(size = 20))


write.csv(t_v2, file="Significant_visit2.csv")

##################For visit 4##########################

set.seed(100)
#matrixData<-t(matrixData)
matrixData_v4<-matrixData[c(colnames_v4),]
matrixData_v4<-matrixData_v4[,colSums(matrixData_v4)>0]
matrixData_v4<-matrixData_v4[rowSums(matrixData_v4)>0,]
matrixData_v4<-t(matrixData_v4)
featureData_v4<-featureData[match(rownames(matrixData_v4),rownames(featureData)),]
featureData_v4<-data.frame(featureData_v4)

otus_metagenomeSeq_v4<-newMRexperiment(matrixData_v4, 
                                       phenoData = AnnotatedDataFrame(metadata_v4), 
                                       featureData = AnnotatedDataFrame(featureData_v4))

otus_metagenomeSeq_filter_v4 = filterData(otus_metagenomeSeq_v4,present=1, depth = 1)
p_v4= cumNormStatFast(otus_metagenomeSeq_filter_v4) 

anxietyData_v4 = cumNorm(otus_metagenomeSeq_v4, p = p_v4)

anxiety_sample_data_v4<-pData(anxietyData_v4)

#anxietyResult_v4<-anxiety_sample_data_v4$Result
bmi_v4<-anxiety_sample_data_v4$BMI
anxietyscores_v4<-anxiety_sample_data_v4$SC5
gestational_weeks_v4<-anxiety_sample_data_v4$Normalized_Gestationalweeks

mod_v4<-model.matrix(~anxietyscores_v4+bmi_4)
settings<-zigControl(maxit = 10, verbose = TRUE)
fit_v4<-fitZig(obj = anxietyData_v4, mod = mod_v4,control = settings,useCSSoffset = FALSE) 

mean_eff_size_v4<-median(calculateEffectiveSamples(fit_v4))
mean_eff_size_v4 = round(mean_eff_size_v4)    #####################mean effective size is 14##############

otus_metagenomeSeq_filter_v4 = filterData(otus_metagenomeSeq_v4, present = mean_eff_size_v4, 
                                          depth = 1) 
p_v4 = cumNormStatFast(otus_metagenomeSeq_filter_v4)    #percentile by which to normalize counts

anxietyData_v4 = cumNorm(otus_metagenomeSeq_filter_v4, p = p_v4)
anxiety_sample_data_v4<-pData(anxietyData_v4)

#anxietyResult_v4<-anxiety_sample_data_v4$Result
bmi_v4<-anxiety_sample_data_v4$BMI
#gestational_weeks_v4<-anxiety_sample_data_v4$Normalized_Gestationalweeks
anxietyscores_v4<-anxiety_sample_data_v4$SC5
gestational_weeks_v4<-anxiety_sample_data_v4$Normalized_Gestationalweeks

norm.factor_v4<- normFactors(anxietyData_v4)                           ### min value is 0.61   and max value is 1.37
norm.factor_v4<-log2(norm.factor_v4/median(norm.factor_v4) + 1)


mod_v4<-model.matrix(~anxietyscores_v4+bmi_v4+norm.factor_v4)

settings<-zigControl(maxit = 10, verbose = TRUE)
fit_v4<-fitZig(obj = anxietyData_v4, mod = mod_v4, control = settings, useCSSoffset = FALSE)

zigFit_v4<-slot(fit_v4,"fit")
finalMod_v4<-slot(fit_v4,"fit")$design

fit_v4<-eBayes(zigFit_v4, trend = TRUE, robust = TRUE)
topTable(fit_v4)

pvalues_v4<-apply(fit_v4$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG_v4<-cbind(fit_v4$coefficients,pvalues_v4)
otus_mG_v4<-data.frame(otus_mG_v4)

#coef_v4<-MRcoefs(fit_v4)

otus_mG_filter_anxiety_v4<-otus_mG_v4[(otus_mG_v4$anxietyscores_v4.1<0.1),]
otus_mG_filter_anxiety_v4<-otus_mG_filter_anxiety_v4[(otus_mG_filter_anxiety_v4$bmi_v4.1>0.1),]
otus_mG_filter_anxiety_v4<-cbind(otus_mG_filter_anxiety_v4, Type = rep("anxiety only", nrow(otus_mG_filter_anxiety_v4)))

otus_mG_filter_anxiety_bmi_v4<-otus_mG_v4[(otus_mG_v4$anxietyscores_v4.1<0.1),]
otus_mG_filter_anxiety_bmi_v4<-otus_mG_filter_anxiety_bmi_v4[(otus_mG_filter_anxiety_bmi_v4$bmi_v4.1<0.1),]
otus_mG_filter_anxiety_bmi_v4<-cbind(otus_mG_filter_anxiety_bmi_v4, Type = rep("anxiety and bmi", nrow(otus_mG_filter_anxiety_bmi_v4)))

otus_v4<-rbind(otus_mG_filter_anxiety_v4,otus_mG_filter_anxiety_bmi_v4)

h_v4<-rownames(otus_v4)

meta_taxa_v4<-data.frame(featureData[c(h_v4),])

t_v4<-cbind.data.frame(meta_taxa_v4$Genus,otus_v4$anxietyscores_v4,otus_v4$Type)
#t_v2<-na.omit(t_v2)

colnames(t_v4)<-c("Genus","Coefficient","Type")
sign_v4<-sapply(as.numeric(t_v4[,"Coefficient"]),sign)
t_v4<-cbind.data.frame(t_v4,sign_v4)
t_v4$Visit<-rep("Visit 4", nrow(t_v4))

ggplot(t_v4,aes(y=Coefficient, x=Genus, col= as.factor(sign_v4)))+xlab(label="Genus")+ylab(label=paste("Coefficient"))+
  coord_flip()+
  geom_point(aes(shape=as.factor(Type)),size = 5)+scale_colour_manual(values = c("-1" = "Darkgreen", "1" = "chocolate"))+ggtitle("Visit 4")+
  geom_abline(slope =0)+guides(col=guide_legend("Sign"),
                               shape=guide_legend("Type"))+
  theme_bw()+theme(aspect.ratio = 1)+
  theme(text = element_text(size = 20))



write.csv(t_v4, file="Significant_visit4.csv")

##########################################plots over time##########################

######################################################################################################################################


################################################################################################################################################################################################################
#Microbial Networks#########################################################

# Co-abundance networks
#--------------------------------------
# summarize at the genus level
#--------------------------------------------
############################################################################################

ps3 = tax_glom(ps2, "Genus", NArm = TRUE)


log_otu<-data.frame(log(otu_table(ps3)+1,2))
colnames(log_otu)<-paste("Seq",1:ncol(log_otu), sep="_")
log_otu<-t(log_otu)

taxa_transformed<-data.frame(tax_table(ps3))
rownames(taxa_transformed)<-paste("Seq",1:nrow(taxa_transformed), sep="_")
#log_otu<-t(log_otu)
#otu_agg<-data.frame(otu_table(ps3))

#samdf3<-metadata[between(metadata$SC5, 0,9),]
#nonanxious_names<-rownames(samdf3)

#otu_transformed_nonanx<-data.frame(log_otu[match(rownames(samdf3),rownames(log_otu)),])
#colnames(otu_transformed_nonanx)<-paste("Seq",1:ncol(otu_transformed_nonanx), sep="_")
#otu_transformed_nonanx<-t(otu_transformed_nonanx)
#taxa_transformed_nonanx<-data.frame(taxa_table[match(rownames(otu_transformed_nonanx),rownames(taxa_table)),])

# Using SparCC

nboot=100
sparCc_list<-list()
sparCc_rand_list<-list()
matrix_raw<-matrix(0, ncol=nboot, nrow=(nrow(log_otu)*(nrow(log_otu)-1))/2)
matrix_random<-matrix_raw
for (iboot in 1:nboot){
  print(iboot)
  iboot_samples<-sample(1:ncol(log_otu),0.9*ncol(log_otu))
  iboot_out<-log_otu[,iboot_samples]
  iboot_out_rand<-apply(iboot_out,2,sample)
  sparcc_iboot<-sparcc(t(iboot_out), iter = 200, inner_iter = 50, th = 0.05)
  sparcc_rand_iboot<-sparcc(t(iboot_out_rand), iter = 200, inner_iter = 50, th = 0.05)
  corr_iboot<-sparcc_iboot$Cor
  corr_rand_iboot<-sparcc_rand_iboot$Cor
  colnames(corr_iboot)<-rownames(corr_iboot)<-colnames(corr_rand_iboot)<-rownames(corr_rand_iboot)<-rownames(log_otu)
  #sparCc_list[[iboot]]<-corr_iboot
  #sparCc_rand_list[[iboot]]<-corr_rand_iboot
  iboot_raw<-melt(corr_iboot)[melt(upper.tri(corr_iboot))[,3]==TRUE,]
  iboot_random<-melt(corr_rand_iboot)[melt(upper.tri(corr_rand_iboot))[,3]==TRUE,]
  matrix_raw[,iboot]<-iboot_raw[,3]
  matrix_random[,iboot]<-iboot_random[,3]
}

rownames(matrix_raw)<-paste(iboot_raw[,1],iboot_raw[,2],sep=".")   #
rownames(matrix_random)<-paste(iboot_random[,1],iboot_random[,2],sep="_")  #

significant<-matrix(0, nrow=nrow(matrix_raw), ncol=3)
for (i in 1:nrow(matrix_raw)){
  print(i)
  ttest_i<-t.test(matrix_raw[i,], matrix_random[i,])
  significant[i,]<-c(ttest_i$estimate, ttest_i$p.value)
}
colnames(significant)<-c("rho", "rho_rand", "pvalue")
rownames(significant)<-rownames(matrix_raw)
padj<-p.adjust(significant[,"pvalue"], method="fdr")
significant<-cbind(significant,padj)

significant_filter<-significant[abs(significant[,"rho"])>=0.3 & significant[,"padj"]<0.001,]
seqNames<-t(sapply(rownames(significant_filter), function(x){r<-strsplit(x,".",fixed=TRUE)[[1]]; return(r)}))
significant_seqonly<-cbind(seqNames, significant_filter)
seqNames<-data.frame(seqNames)

significant_seqonly<-cbind(taxa_transformed[seqNames[,1],"Genus"],
                                  taxa_transformed[seqNames[,2],"Genus"],
                                  significant_seqonly)

colnames(significant_seqonly)<-c("From","To","FromSeq","ToSeq","rho","rho_random","p","padj")
sign_edge<-sapply(as.numeric(significant_seqonly[,"rho"]),sign)
nodes_genus<-significant_seqonly[,3]

phylum_all<-taxa_transformed[c(nodes_genus),]
phylum_all<-phylum_all$Phylum

significant_seqonly<-cbind.data.frame(significant_seqonly,sign_edge,phylum_all)


write.table(significant_seqonly,"C:/Akanksha/LAB WORK/fastq_metagenomeseq/Edges_SparCc_Phenotypes.txt",
            sep="\t", col.names =TRUE,row.names = FALSE)

write.csv(significant_seqonly, file="microbiome_networks.csv")

######################################################################################################################

##########################################################################################################################################################################################


##########################picrust2 Tutorial########################################
seqtab.nochim_filtered<-data.frame(otu_table(ps2))
#seqtab.nochim_filtered<-data.frame(seqtab.nochim_filtered)

# Table with the sequences that passed all the filtering (low abundance, contamination)
seqtab.nochim_filtered<-t(seqtab.nochim_filtered)
seqtab.nochim_filtered<-seqtab.nochim_filtered[rowSums(seqtab.nochim_filtered)>0,]
seqtab.nochim_filtered<-seqtab.nochim_filtered[,colSums(seqtab.nochim_filtered)>0]
#colnames(seqtab.nochim_filtered)<-substring(colnames(seqtab.nochim_filtered),3,17)
#colnames(seqtab.nochim_filtered)<-paste("Sample", colnames(seqtab.nochim_filtered), sep="")
#colnames(seqtab.nochim_filtered)<-gsub(".","",colnames(seqtab.nochim_filtered), fixed=TRUE)

# This is the table that you need to input to Picrust2
subOTU_table_filtered<-cbind(OTUID=paste("seq", seq(1,nrow(seqtab.nochim_filtered),by=1), sep="_"),
                                                  seqtab.nochim_filtered)
rownames(subOTU_table_filtered)<-NULL
subOTU_table_filtered<-data.frame(subOTU_table_filtered)
names(subOTU_table_filtered)[1]<-"OTU_ID"
#rownames(subOTU_table_filtered)<-NULL

#subOTU_table_filtered_PIPHILLIN<-t(subOTU_table_filtered_PIPHILLIN)
#subOTU_filtered_biom<-make_biom(subOTU_table_filtered, sample_metadata = NULL, observation_metadata = NULL,
                               # id = NULL, matrix_element_type = "int")
#write_biom(subOTU_filtered_biom, "C:/Akanksha/LAB WORK/fastq_metagenomeseq/subOTU_filtered_biom")

write.table(subOTU_table_filtered,
            file="subOTU_filtered.txt",sep = "\t")

# This will generate the Fasta file that you need

#seqtab.nochim.fasta<-data.frame(otu_table(ps2))
seqtab.nochim.normalized_fasta<-c()
linenumber=0
for (i in 1:nrow(seqtab.nochim_filtered)){
  linenumber=linenumber+1
  seqtab.nochim.normalized_fasta[linenumber]<-paste(">seq",i,sep = "_")
  linenumber=linenumber+1
  seqtab.nochim.normalized_fasta[linenumber]<-rownames(seqtab.nochim_filtered)[i]
  seqtab.nochim.normalized_fasta<-noquote(seqtab.nochim.normalized_fasta)
  
}
write(seqtab.nochim.normalized_fasta,
      file="C:/Akanksha/LAB WORK/fastq_metagenomeseq/seqtabnochimnormalized.fasta")

#Save metadata into tsv file
write.table(metadata,"C:/Akanksha/LAB WORK/fastq_metagenomeseq/metadata.tsv")


############This is not present in my output#########------rxn_abund_table_unnorm----------------------------------

abund_table_unnormal<-read.csv("C:/Akanksha/LAB WORK/fastq_metagenomeseq/pred_metagenome_unstrat_GeneID.csv")
abund_table_unnormal<- data.frame(abund_table_unnormal[,-1], row.names = abund_table_unnormal[,1])


# I selected the IDs of the samples I was interested in comparing using the phyloseq container that has the CSS normalized
# samples
anx_IDs<-rownames(sample_data(ps2))[which(sample_data(ps2)[,"Result"]== "ANC")]
Noanx_IDs<-rownames(sample_data(ps2))[which(sample_data(ps2)[,"Result"]== "AHC")]

#Normalized data

#abund_table_unnorm_filtered_matrix<-matrix(as.numeric(unlist(abund_table_unnorm_filtered)),
 #                                             ncol=ncol(abund_table_unnorm_filtered),
  #                                            nrow=nrow(abund_table_unnorm_filtered))
#rownames(abund_table_unnorm_filtered_matrix)<-rownames(abund_table_unnorm_filtered)
#colnames(abund_table_unnorm_filtered_matrix)<-colnames(abund_table_unnorm_filtered)

#abund_table_norm_filtered_den<-matrix(rep(colSums(abund_table_unnorm_filtered), nrow=nrow(abund_table_unnorm_filtered)),
 #                                        ncol=ncol(abund_table_unnorm_filtered), nrow=nrow(abund_table_unnorm_filtered), 
  #                                       byrow = TRUE)
#abund_table_norm_filtered_matrix<-abund_table_unnorm_filtered_matrix/abund_table_norm_filtered_den
#abund_table_norm_filtered_matrix<-abund_table_norm_filtered_matrix[rowSums(abund_table_norm_filtered_matrix)>0,]
#abund_table_norm_filtered_matrix<-abund_table_norm_filtered_matrix[,colSums(abund_table_norm_filtered_matrix)>0]

abund_table_clr<-data.frame(clr(abund_table_unnormal))


#abund_table_clr<-abund_table_clr[,colSums(abund_table_clr)>0]
#abund_table_unnorm_filtered<-abund_table_unnorm_filtered[rowSums(abund_table_unnorm_filtered)>0,]



# I used a wilcox-test to identify the genes that encode for enzymes that are difference 
# between the cases (MDD) and controls (noMDD)
#ko_wilcox.test<-NULL
#for (iko in 1:nrow(abund_table_norm_filtered_matrix)){
 # print(iko)
  #iko_wilcox<-wilcox.test(as.numeric(abund_table_norm_filtered_matrix[iko,match(anx_IDs,colnames(abund_table_norm_filtered_matrix))]), 
    #                      as.numeric(abund_table_norm_filtered_matrix[iko,match(Noanx_IDs,colnames(abund_table_norm_filtered_matrix))]),
   #                       conf.int = TRUE,conf.level = 0.9)
  #ko_wilcox.test<-rbind(ko_wilcox.test,
   #                     cbind(rownames(abund_table_norm_filtered_matrix)[iko], iko_wilcox[[8]][1],
    #                          iko_wilcox[[8]][2],iko_wilcox[3]))
#}
#rownames(ko_wilcox.test)<-NULL
#colnames(ko_wilcox.test)<-c("KO_ID", "LCI", "UCI","pvalue")
#pvalue_fdrcorrected<-p.adjust(ko_wilcox.test[,"pvalue"], method="fdr")
#ko_wilcox.test<-data.frame(ko_wilcox.test,pvalue_fdrcorrected)
#ko_wilcox.test<-ko_wilcox.test[order(ko_wilcox.test$pvalue_fdrcorrected),]


gene_ID<-c("K01667","K04103","K07130","K00128")
abund_table_sp<-abund_table_clr[c(gene_ID),]
abund_table_sp<-abund_table_sp[,colSums(abund_table_sp)>0]  #####Since most values are 0, it is necessary to do this step


anx_IDs_abund<-anx_IDs[anx_IDs%in%colnames(abund_table_sp)]
Noanx_IDs_abund<-Noanx_IDs[Noanx_IDs%in%colnames(abund_table_sp)]

ko_wilcox.test_1<-NULL
for (iko_1 in 1:nrow(abund_table_sp)){
  print(iko_1)
  iko_wilcox_1<-wilcox.test(as.numeric(abund_table_sp[iko_1,match(anx_IDs_abund,colnames(abund_table_sp))]), 
                          as.numeric(abund_table_sp[iko_1,match(Noanx_IDs_abund,colnames(abund_table_sp))]),
                          conf.int = TRUE,conf.level = 0.9)
  ko_wilcox.test_1<-rbind(ko_wilcox.test_1,
                        cbind(rownames(abund_table_sp)[iko_1], iko_wilcox_1[[8]][1],
                              iko_wilcox_1[[8]][2],iko_wilcox_1[3]))
}

rownames(ko_wilcox.test_1)<-NULL
colnames(ko_wilcox.test_1)<-c("KO_ID", "LCI", "UCI","pvalue")
pvalue_fdrcorrected_1<-p.adjust(ko_wilcox.test_1[,"pvalue"], method="fdr")
ko_wilcox.test_1<-data.frame(ko_wilcox.test_1,pvalue_fdrcorrected_1)
ko_wilcox.test_1<-ko_wilcox.test_1[order(ko_wilcox.test_1$pvalue_fdrcorrected_1),]

abund_table_mg<-data.frame(t(abund_table_sp))

SampleID<-rownames(abund_table_mg)
rownames(abund_table_mg)<-NULL
abund_table_mg<-cbind.data.frame(SampleID,abund_table_mg)


abund_table_mg<-merge(abund_table_mg,sample_data(ps2)[,c("SC5")], by.x="SampleID",
                      by.y=0, all.x=TRUE)
abund_table_mg$Result[abund_table_mg$SC5 > 9] <- "ANC"
abund_table_mg$Result[abund_table_mg$SC5 <= 9] <- "AHC"

fig_7<-ggplot(abund_table_mg,aes(x = Result, y = K01667))+geom_boxplot(width=0.5,lwd=1.5,color= "Darkgreen")+theme_bw()+theme(aspect.ratio = 1)
fig_7<-fig_7+theme(text=element_text(size=20))
plot(fig_7)

fig_8<-ggplot(abund_table_mg,aes(x = Result, y = K04103))+geom_boxplot(width=0.5,lwd=1.5,color= "Chocolate")+theme_bw()+theme(aspect.ratio = 1)
fig_8<-fig_8+theme(text=element_text(size=20))
plot(fig_8)

fig_9<-ggplot(abund_table_mg,aes(x = Result, y = K07130))+geom_boxplot(width=0.5,lwd=1.5,color= "Maroon1")+theme_bw()+theme(aspect.ratio = 1)
fig_9<-fig_9+theme(text=element_text(size=20))
plot(fig_9)

fig_10<-ggplot(abund_table_mg,aes(x = Result, y = K00128))+geom_boxplot(width=0.5,lwd=1.5,color= "lightslateblue")+theme_bw()+theme(aspect.ratio = 1)
fig_10<-fig_10+theme(text=element_text(size=20))
plot(fig_10)

########Pathways differentially abundant###################

pathway_abund_table<-read.csv("C:/Akanksha/LAB WORK/fastq_metagenomeseq/path_abun_unstrat.csv")
pathway_abund_table<- data.frame(pathway_abund_table[,-1], row.names = pathway_abund_table[,1])

pathway_table_unnorm_filtered_matrix<-matrix(as.numeric(unlist(pathway_abund_table)),
                                           ncol=ncol(pathway_abund_table),
                                           nrow=nrow(pathway_abund_table))

rownames(pathway_table_unnorm_filtered_matrix)<-rownames(pathway_abund_table)
colnames(pathway_table_unnorm_filtered_matrix)<-colnames(pathway_abund_table)


pathway_table_norm_filtered_den<-matrix(rep(colSums(pathway_table_unnorm_filtered_matrix), nrow=nrow(pathway_table_unnorm_filtered_matrix)),
                                      ncol=ncol(pathway_table_unnorm_filtered_matrix), nrow=nrow(pathway_table_unnorm_filtered_matrix), 
                                      byrow = TRUE)

pathway_table_norm_filtered_matrix<-pathway_table_unnorm_filtered_matrix/pathway_table_norm_filtered_den

pathway_table_norm_filtered_matrix<-pathway_table_norm_filtered_matrix[,colSums(pathway_table_norm_filtered_matrix)>0]


#pathway_table_clr<-data.frame(clr(pathway_abund_table))
#pathway_table_clr<-pathway_abund_table[rowSums(pathway_abund_table)>0,]
#pathway_table_clr<-pathway_table_clr[,colSums(pathway_table_clr)>0]
#pathway_abund_table<-pathway_abund_table[rowSums(pathway_abund_table)>0,]
#pathway_sample<-colnames(pathway_abund_clr)
#metadata_path<-metadata[c(pathway_sample),]
set.seed(100)



pathway_ID<-c("TRPSYN-PWY","PWY-5022","CENTFERM-PWY")

pathway_abund_table_sp<-pathway_table_norm_filtered_matrix[c(pathway_ID),]

pathway_abund_table_sp<-pathway_abund_table_sp[,colSums(pathway_abund_table_sp)>0]



#anx_IDs_path<-anx_IDs[anx_IDs%in%colnames(pathway_abund_table_sp)]
#Noanx_IDs_path<-Noanx_IDs[Noanx_IDs%in%colnames(pathway_abund_table_sp)]

ko_wilcox.test_2<-NULL
for (iko_2 in 1:nrow(pathway_abund_table_sp)){
  print(iko_2)

  iko_wilcox_2<-wilcox.test(as.numeric(pathway_abund_table_sp[iko_2,match(anx_IDs,colnames(pathway_abund_table_sp))]), 
                            as.numeric(pathway_abund_table_sp[iko_2,match(Noanx_IDs,colnames(pathway_abund_table_sp))]),
                            conf.int = TRUE,conf.level = 0.9)

  ko_wilcox.test_2<-rbind(ko_wilcox.test_2,
                          cbind(rownames(pathway_abund_table_sp)[iko_2], iko_wilcox_2[[8]][1],
                                iko_wilcox_2[[8]][2],iko_wilcox_2[3]))
}

rownames(ko_wilcox.test_2)<-NULL
colnames(ko_wilcox.test_2)<-c("KO_ID", "LCI", "UCI","pvalue")
pvalue_fdrcorrected_2<-p.adjust(ko_wilcox.test_2[,"pvalue"], method="fdr")
ko_wilcox.test_2<-data.frame(ko_wilcox.test_2,pvalue_fdrcorrected_2)
ko_wilcox.test_2<-ko_wilcox.test_2[order(ko_wilcox.test_2$pvalue_fdrcorrected_2),]

pathway_table_mg<-data.frame(t(pathway_abund_table_sp))

SampleID<-rownames(pathway_table_mg)
rownames(pathway_table_mg)<-NULL
pathway_table_mg<-cbind.data.frame(SampleID,pathway_table_mg)


pathway_table_mg<-merge(pathway_table_mg,sample_data(ps2)[,c("SC5")], by.x="SampleID",
                      by.y=0, all.x=TRUE)
pathway_table_mg$Result[pathway_table_mg$SC5 > 9] <- "ANC"
pathway_table_mg$Result[pathway_table_mg$SC5 <= 9] <- "AHC"



fig_11<-ggplot(pathway_table_mg,aes(x = Result, y = TRPSYN.PWY))+geom_boxplot(width=0.5,lwd=1.5,color= "Darkgreen")+theme_bw()+theme(aspect.ratio = 1)
fig_11<-fig_11+xlab(label="Result")+ylab(label=paste("Tryptophan Synthesis pathway"))+theme(text=element_text(size=20))
plot(fig_11)


fig_12<-ggplot(pathway_table_mg,aes(x = Result, y = PWY.5022))+geom_boxplot(width=0.5,lwd=1.5,color= "chocolate")+theme_bw()+theme(aspect.ratio = 1)
fig_12<-fig_12+xlab(label="Result")+ylab(label=paste("GABA degradation pathway"))+theme(text=element_text(size=20))
plot(fig_12)


fig_13<-ggplot(pathway_table_mg,aes(x = Result, y = CENTFERM.PWY))+geom_boxplot(width=0.5,lwd=1.5,color= "Maroon1")+theme_bw()+theme(aspect.ratio = 1)
fig_13<-fig_13+xlab(label="Result")+ylab(label=paste("Butyrate Synthesis pathway"))+theme(text=element_text(size=20))
plot(fig_13)








################################################################################################################################################################


GAD7_scores<-metadata[,46:52]
GAD7_scores<-cbind(GAD7_scores,SC5)
#study_sample_names<-paste(metadata$Study,metadata$Subject)
#rownames(GAD7_scores)<-study_sample_names
GAD7_scores_scaled<-scale(GAD7_scores,center=TRUE,scale=TRUE)


heatmap.2(GAD7_scores_scaled, Colv= FALSE,dendogram="row",col =colorRampPalette(c("Darkgreen","Maroon1")),margins = c(5,20),key = TRUE,trace="none",key.xlab = "Z-score",key.ylab = "",
           cexCol = 0.5,cexRow = 1)


####################################################

###########

prevdf_3 = apply(X = otu_table(ps_v1),
               MARGIN = ifelse(taxa_are_rows(ps_v1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf_3 = data.frame(Prevalence = prevdf_3,
                    TotalAbundance = taxa_sums(ps_v1),
                    tax_table(ps_v1))

prevdf_3<-subset(prevdf_3, Phylum %in% get_taxa_unique(ps_v1, "Phylum"))
ggplot(prevdf_3, aes(TotalAbundance, Prevalence / nsamples(ps_v1),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 4, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(text = element_text(size = 10))+theme_bw()




prevdf_4 = apply(X = otu_table(ps_v2),
                 MARGIN = ifelse(taxa_are_rows(ps_v2), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

prevdf_4 = data.frame(Prevalence = prevdf_4,
                      TotalAbundance = taxa_sums(ps_v2),
                      tax_table(ps_v2))

prevdf_4<-subset(prevdf_4, Phylum %in% get_taxa_unique(ps_v2, "Phylum"))
ggplot(prevdf_4, aes(TotalAbundance, Prevalence / nsamples(ps_v2),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 4, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(text = element_text(size = 10))+theme_bw()



###############Maaslin2 method###########################################



ps3_otu<-data.frame(otu_table(ps3))

ps3_taxa<-data.frame(tax_table(ps3))

ps3_genus<-ps3_taxa$Genus

colnames(ps3_otu)<-ps3_genus

ps3_meta<-data.frame(sample_data(ps3))



write.table(ps3_otu, file = "matrixData_1.txt", sep = "\t")
write.table(ps3_meta, file="metadata_1.txt", sep ="\t")

df_input_data = read.table("C:/Akanksha/LAB WORK/fastq_metagenomeseq/matrixData_1.txt", sep = "\t",row.names = 1)

df_input_metadata = read.table("C:/Akanksha/LAB WORK/fastq_metagenomeseq/metadata_1.txt", sep = "\t", row.names = 1)

total<-rowSums(df_input_data)

rel_abund<-ave(df_input_data, FUN = function(x){x / sum(x)})

for (i in 1:ncol(df_input_data)){
  
}





fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata,
  output = "demo_output",
  fixed_effects = c("SC5","Normalized_Gestationalweeks"), random_effects = c("Subject"), normalization = c("CSS"))




slide_result <- sliding_spliner(data = df_input_metadata,
                                xvar = 'Normalized_Gestationalweeks', yvar = 'rel_abund',
                                category = 'Result', cases = 'Subject',
                                test_density = 10, cut_low = 7,
                                set_spar = 0.5)




slide_result <- sliding_spliner(data = df_input_data,
                                xvar = 'Normalized_Gestational weeks', yvar = 'relative_abundance',
                                category = 'antibiotics_y_n', cases = 'baby_id',
                                test_density = 10, cut_low = 7,
                                set_spar = 0.5)