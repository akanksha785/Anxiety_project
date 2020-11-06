library(readxl)
library(dada2)
library(DECIPHER)
library(phangorn)
library(phyloseq)
library(gridExtra)
library(dplyr)
library(plyr)
library(ggplot2)
library(tidyverse)
library(csv)
library(decontam)
library(metagenomeSeq)
library(splinectomeR)
library(reshape2)
library(tibble)
library(tidyr)
library(vegan)
library(igraph)
library(HMP)
library(dendextend)
library(rbin)
library(stats)
library(ggpubr)
library(nloptr)
library(coda)
library(gss)
library(remotes)
library(devtools)
library(SpiecEasi)


df<-'C:/Akanksha/LAB WORK/FastQ files'
list.files(df)

##Sorting Forward and Reverse reads

fwd_df<-sort(list.files(df, pattern = "_R1.fastq"))
rvs_df<-sort(list.files(df, pattern = "_R2.fastq"))

##Extracting sample names
sample.names <- sapply(strsplit(basename(fwd_df), "_"), `[`, 1)


#Inspecting profile reads

##Forward Reads
plotQualityProfile(fwd_df[1:2])
##plotQualityProfile(fwd_df[66:67])
##plotQualityProfile(fwd_df[131:132])

##Reverse Reads
plotQualityProfile(rvs_df[1:2])
##plotQualityProfile(rvs_df[66:67])
##plotQualityProfile(rvs_df[130:131])


# Place filtered files in filtered/ subdirectory
filtfwd_df <- file.path(df, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtrvs_df<- file.path(df, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtfwd_df) <- sample.names
names(filtrvs_df) <- sample.names


out<- filterAndTrim(fwd_df, filtfwd_df,rvs_df,filtrvs_df, truncLen=150,maxN=0,maxEE=2,truncQ = 2 ,rm.phix=TRUE,compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE



errFwd <- learnErrors(filtfwd_df, multithread=TRUE)
errRvs <- learnErrors(filtrvs_df, multithread=TRUE)
plotErrors(errFwd, nominalQ=TRUE)
plotErrors(errRvs, nominalQ=TRUE)

derepFs <- derepFastq(filtfwd_df)
derepRs <- derepFastq(filtrvs_df)

sam.names <- sapply(strsplit(basename(filtfwd_df), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names


dadaFwd <- dada(derepFs, err=errFwd, multithread = TRUE)
dadaRvs <- dada(derepRs, err=errRvs, multithread = TRUE)
dadaFwd[[1]]


plotErrors(dadaFwd)
plotErrors(dadaRvs)


mergers <- mergePairs(dadaFwd, filtfwd_df, dadaRvs, filtrvs_df, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
dim(seqtab)           ###this shows the number of ASV which is 6091
table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#rownames(seqtab.nochim)<-1:nrow(seqtab.nochim)
#id<-samdf$SampleID
#seqtab.nochim["SampleID"]<-id
dim(seqtab.nochim)   ######after removing chimeric ASV which is 4137
sum(seqtab.nochim)/sum(seqtab)

taxa <- assignTaxonomy(seqtab.nochim, "C:/Akanksha/LAB WORK/FastQ files/Silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

taxa_1<- addSpecies(taxa,"C:/Akanksha/LAB WORK/FastQ files/Silva/silva_species_assignment_v132.fa.gz")

######################################################################################################################################################
#rownames(taxa_1)<-1:nrow(taxa_1)
#taxa_1.print <- taxa_1 # Removing sequence rownames for display only

#rownames(taxa_1.print) <- NULL
#head(taxa_1.print)
#unname(head(taxa_1))
taxonomy<-array(taxa_1)
taxonomy0<-data.frame(taxonomy)

#Dividing taxonomy into kingdom,phylum,family,order,class,genus,species

Kingdom<-taxonomy0[1:4075, ]
Phylum<-taxonomy0[4076:8150, ]
Class<-taxonomy0[8151:12225, ]
Order<-taxonomy0[12226:16300, ]
Family<-taxonomy0[16301:20375, ]
Genus<-taxonomy0[20376:24450, ]
Species<-taxonomy0[24451:28525, ]


##library(msa)
mapfile0<-read.csv('C:/Akanksha/LAB WORK/Dataset/mapFile_UNC_09012020.csv', header = TRUE, stringsAsFactors = FALSE)
#sample_EPDS<-mapfile0$SC0
#sample_PSS<-mapfile0$SC4
#sample_GAD7<-mapfile0$SC5

samdf<-mapfile0[1:179, ]
samdf1 <- data.frame(samdf[,-1], row.names = samdf[,1])
samdf2<-sample_data(samdf1)


#a<-samdf %>%
#b<-remove_rownames() %>%
#c<-column_to_rownames(var = 'SampleID')


#seqtab.nochim<-cbind(SampleID = rownames(seqtab.nochim), seqtab.nochim)
#rownames(seqtab.nochim) <- 1:nrow(seqtab.nochim)


#mapfile1<-sort.int(samdf$SampleID,decreasing = FALSE, na.last = NA,index.return = FALSE)
#samdf1<-data.frame(sample_EPDS,sample_PSS,sample_GAD7)
#rownames(samdf)<-id
#mapfile3<-data.frame(mapfile1,samdf)
#mapfile4<-subset(samdf,select=samdf[ ,0])
#mapfile5<-rename(mapfile4,c("mapfile2" = "SampleID"))
#rownames(mapfile5)<-1:nrow(mapfile5)

seqs<-getSequences(seqtab.nochim)
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


ps<-phyloseq(tax_table(taxa_1),sample_data(samdf2),otu_table(seqtab.nochim,taxa_are_rows = FALSE), phy_tree(fitGTR_1$tree))
rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL) 

#df1 <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
#df1$LibrarySize <- sample_sums(ps)
#df1 <- df1[order(df1$LibrarySize),]
#df1$Index <- seq(nrow(df1))
#Sample_or_Control<-samdf2$Study
#ggplot(data=df1, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

###################################################################################################################3

negative<-samdf2$Sample_or_Control == "Control"
contaminant<-isContaminant(otu_table(ps), neg = negative, method = "prevalence", threshold = 0.05, normalize = TRUE, detailed = TRUE)


#contaminant_1<-cbind(contaminant, new_col = samdf$SampleID)

#rownames(contaminant)<-paste("Seq",1:ncol(seqtab.nochim), sep="_")

table(contaminant$contaminant)

#Filtering out true contaminant sequences
possibleCont<-filter(contaminant, contaminant == "TRUE")
possibleCont1<-rownames(possibleCont)

yy<-data.frame(otu_table(ps))
#colnames(yy)<-paste("Seq",1:ncol(yy), sep="_")

#colnames(seqtab.nochim)<-paste("Seq",1:ncol(seqtab.nochim), sep="_")

#Identify the manes of the control and cases samples

cases<-rownames(sample_data(ps)[which(sample_data(ps)[,"Study"]=="KMP" | sample_data(ps)[,"Study"]=="MPI"),])
controls<-rownames(sample_data(ps)[which(sample_data(ps)[,"Study"]!="KMP" & sample_data(ps)[,"Study"]!="MPI"),])

#Calculate the maximum and minimum for sequences

casesMax<-apply(yy[cases,possibleCont1],2,max)
casesMin<-apply(yy[cases,possibleCont1],2,min)
controlMax<-apply(yy[controls,possibleCont1],2,max)
controlMin<-apply(yy[controls,possibleCont1],2,min)

#Dataframe 

contaminationValues<-cbind(casesMax,casesMin,controlMax,controlMin)

#rownames(contaminationValues)<-possibleCont1

contamination_final<-data.frame(contaminationValues)
aa<-contamination_final$casesMax
bb<-contamination_final$controlMax
cc<-as.numeric(bb)
#dd<-vector()
xx<-filter(contamination_final, aa<3*cc) #if cases maximum is less than 3 times controls max, then it is a contaminant
zz<-rownames(xx)
#pp<-data.frame(zz)
allTaxa=taxa_names(ps)
allTaxa<-allTaxa[!(allTaxa %in% zz)]
ps.noncontam = prune_taxa(allTaxa,ps)

rank_names(ps.noncontam)
table(tax_table(ps.noncontam)[, "Phylum"], exclude = NULL) 

#noncontam_otu<-otu_table(ps.noncontam)

#table_otu<-otu_table(prune_taxa(taxa_sums(ps.noncontam)>0, ps.noncontam))
#p1<-phyloseq(tax_table(taxa_1),sample_data(samdf2),otu_table(table_otu,taxa_are_rows = FALSE), phy_tree(fitGTR_1$tree))
#ps11<- prune_taxa(table_otu,seqtab.nochim)

#table_otu_1<-data.frame(table_otu)

#colnames(yy)<-paste("Seq",1:ncol(yy), sep="_")
#colnames(seqtab.nochim)<-paste("Seq",1:ncol(seqtab.nochim), sep="_")
#colnames(yy)<-paste("Seq",1:ncol(yy), sep="_")
#rownames(taxa_1)<-paste("Seq",1:nrow(taxa_1), sep="_")
#noncontam<-seqtab.nochim[ ,-c(7,22,31,47,55,61,70,76,78,129,134,136,142,143,152,155,172,191,197,200,224,227,232,240,243,247,283,286,292,308,318,328,333,342,353,381,393,399,421,427,429,467,489,500,502,533,563,568,599,601,638,646,665,671,747,767,793,816,892,895,898,924,954,955,958,963,978,1011,1012,1030,1129,1140,1177,1193,1231,1289,1330,1400,1410,1431,1457,1466,1552,1566,1690,2081,2127,2205,2778)]
#rownames(noncontam)<-1:nrow(noncontam)
#rownames(taxa_1)<-paste("Seq",1:nrow(taxa_1), sep="_")
#colnames(taxa_2)<- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
#taxa_3<-t(taxa_2)
#a_tab = otu_table(as.matrix(noncontam), taxa_are_rows=TRUE) 
#t_tab = tax_table(as.matrix(taxa_2)) 

#######################################################################################################################
#keep.exprs <- filterByExpr(s1, min.count=1)
#filt1 <- data2[keep.exprs,]
#dim(filt1)
#g<-data.frame(data2)

samdf3<-samdf1[(samdf1$Sample_or_Control=="Sample"),]
samdf4<-sample_data(samdf3)
#samdf5<-samdf[(samdf$Sample_or_Control=="Sample"),]
#samdf6<-sample_data(samdf5)

ps0<-phyloseq(tax_table(taxa_1),sample_data(samdf4),otu_table(ps.noncontam,taxa_are_rows = FALSE), phy_tree(fitGTR_1$tree))


####ps1---->removed all NAs from phylum########################

ps1<- subset_taxa(ps0, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))  
table(phyloseq::tax_table(ps1)[, "Phylum"])
plot_bar(ps1, fill = "Phylum") + theme(legend.position="bottom")
tree<-plot_tree(ps1, nodelabf = nodeplotblank, color = NULL,
                shape = NULL, label.tips = "Phylum", sizebase = 5,
                text.size = 4, base.spacing = 0.01,
                ladderize = FALSE, plot.margin = 0.1)
plot(tree)


# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf_1 = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps1),
                      tax_table(ps1))

plyr::ddply(prevdf_1, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

filterPhyla = c("Deferribacteres","Dependentiae","Kiritimatiellaeota","Patescibacteria","Planctomycetes","Spirochaetes")              ##Removing these two phyla because their mean is too low and also its sum

###########ps2--------->Removing phyla because their mean and sum  is too low.############
ps2 = subset_taxa(ps1, !Phylum %in% filterPhyla)
tree1<-plot_tree(ps2, nodelabf = nodeplotblank, color = NULL,
                 shape = NULL, label.tips = "Genus",
                 text.size = 4, base.spacing = 0.02, sizebase = 5,
                 ladderize = FALSE, plot.margin = 0.1)
plot(tree1)

# Subset to the remaining phyla
prevdf_2 = subset(prevdf_1, Phylum %in% get_taxa_unique(ps2, "Phylum"))
ggplot(prevdf_2, aes(TotalAbundance, Prevalence / nsamples(ps2),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps2)
 
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf_2)[(prevdf_2$Prevalence >= prevalenceThreshold)]

##########ps3-------->keeping taxa which is more than 5% of total samples####################################
ps3 = prune_taxa(keepTaxa, ps2)

s1<-t(otu_table(ps3))
s2<-data.frame(s1)

for (i in 1:nrow(s2)) {
  for (j in 1:ncol(s2)) {
    if(s2[i,j]<=5){
      s2[i,j]<-0
    }
    
  }
  
}

data<-s2[rowSums(s2)>0,]
data1<-data[,colSums(data)>0]
data2<-t(data1)

ps4<-phyloseq(tax_table(taxa_1),sample_data(samdf4),otu_table(data2,taxa_are_rows = FALSE), phy_tree(fitGTR_1$tree))

r<-rownames(data2)
#pp<-data.frame(zz)
allTaxa=taxa_names(ps4)
allTaxa<-allTaxa[!(allTaxa %in% r)]
ps.rem = prune_taxa(allTaxa,ps4)

samdf5<-sample_data(ps.rem)

#rownames(taxa_1)<-paste("Seq",1:nrow(taxa_1), sep="_")
#colnames(data2)<-paste("Seq",1:ncol(data2), sep="_")

#data2<-t(data2)


#ps50<-phyloseq(tax_table(taxa_1),sample_data(samdf4),otu_table(data2,taxa_are_rows = TRUE), phy_tree(fitGTR_1$tree))



#Grouping by type 
samdf5$Type<-factor(samdf5$Type,
                    levels=c("Minimal", "Mild", "Moderate", "Severe"))
Anxiety_levels<-samdf5$Type

#Grouping by visits
samdf5$V<-factor(samdf5$V,levels=c("1","2","4"))
samdf5$visit<-samdf5$V

tab <- microbiome::alpha(ps.rem, index = "all")

#Add the diversity table to metadata

shannon<- tab$diversity_shannon


k1<-cbind.data.frame(samdf5,shannon)


#p1 <- ggviolin(samdf5, x = "Result", y = "chao1",add = "boxplot", fill = "Result") 
#print(p1)


##########################Tree plot for Phylum level###########################################################

length(get_taxa_unique(ps.rem, taxonomic.rank = "Phylum"))
ps5 = tax_glom(ps.rem, "Phylum", NArm = TRUE)
tree2<-plot_tree(ps5, nodelabf = nodeplotblank, color = NULL,
                 shape = NULL, label.tips = "Phylum",sizebase = 5,
                 text.size = 4, base.spacing = 0.02,
                 ladderize = FALSE, plot.margin = 0.2)
plot(tree2)

length(get_taxa_unique(ps.rem, taxonomic.rank = "Class"))
ps6= tax_glom(ps.rem, "Class", NArm = TRUE)
tree3<-plot_tree(ps6, nodelabf = nodeplotblank, color = NULL,
                 shape = NULL, label.tips = "Class", sizebase = 5,
                 text.size = 4, base.spacing = 0.02,
                 ladderize = FALSE, plot.margin = 0.2)
plot(tree3)

length(get_taxa_unique(ps.rem, taxonomic.rank = "Order"))
ps7 = tax_glom(ps.rem, "Order", NArm = TRUE)
tree4<-plot_tree(ps7, nodelabf = nodeplotblank, color = NULL,
                 shape = NULL, label.tips = "Order", sizebase = 5,
                 text.size = 4, base.spacing = 0.02,
                 ladderize = FALSE, plot.margin = 0.2)

plot(tree4)

length(get_taxa_unique(ps.rem, taxonomic.rank = "Family"))
ps8 = tax_glom(ps.rem, "Family", NArm = TRUE)
tree5<-plot_tree(ps8, nodelabf = nodeplotblank, color = NULL,
                 shape = NULL, label.tips = "Family", sizebase = 5,
                 text.size = 4, base.spacing = 0.02,
                 ladderize = FALSE, plot.margin = 0.2)
plot(tree5)

length(get_taxa_unique(ps.rem, taxonomic.rank = "Genus"))
ps9 = tax_glom(ps.rem, "Genus", NArm = TRUE)
tree6<-plot_tree(ps9, nodelabf = nodeplotblank, color = NULL,
                 shape = NULL, label.tips = "Genus", sizebase = 5,
                 text.size = 4, base.spacing = 0.02,
                 ladderize = FALSE, plot.margin = 0.2)
plot(tree6)

length(get_taxa_unique(ps.rem, taxonomic.rank = "Species"))
ps10 = tax_glom(ps.rem, "Species", NArm = TRUE)
tree7<-plot_tree(ps10, nodelabf = nodeplotblank, color = NULL, 
                 shape = NULL, label.tips = "Species", sizebase = 5,
                 text.size = 4, base.spacing = 0.02,
                 ladderize = FALSE, plot.margin = 0.2)
plot(tree7)

###############################################################################################################################################################################################################################################################

model<-lm(shannon ~ SC5, data=k1)
print(model)
summary(model) 

model<-lm(shannon ~ SC5+Age+as.numeric(GestationalWeeks_sinceConception)+factor(Subject), data=k1)
print(model)
summary(model) 

#Alpha diversity by Type
alpha_plot_1<-plot_richness(ps.rem, x = "SC5",color = "Type",measures=c("Shannon"))+geom_point(size = 4) +labs(title= "Shannon diversity in women using GAD-7 scores")+theme(aspect.ratio=1)                                                                                  
plot(alpha_plot_1)



model1<-lm(shannon ~ Visit , data=k1)
print(model1)
summary(model1) 

##Alpha diversity by visit
alpha_plot_2<-plot_richness(ps.rem,x = "Visit" ,color = "Type",measures=c("Shannon"))+geom_point(size = 4)+labs(title= "Shannon diversity in women for different trimester")+theme(aspect.ratio=1) 
plot(alpha_plot_2)



model2<-lm(shannon ~ GestationalWeeks_sinceConception , data=k1)
print(model2)
summary(model2) 

#Alpha diversity by Gestational weeks
alpha_plot_3<-plot_richness(ps.rem,x = "GestationalWeeks_sinceConception",color = "Type",measures = c("Shannon"))+ geom_point(size = 4)+labs(title= "Chao1 diversity in women for different Gestational weeks")
plot(alpha_plot_3)



#Generate a data.frame with adiv measures
adiv <- data.frame("Observed" = phyloseq::estimate_richness(ps.rem, measures = "Observed"),"Shannon" = phyloseq::estimate_richness(ps.rem, measures = "Shannon"),tree = phyloseq::phy_tree(ps.rem$tree),"Result" = phyloseq::sample_data(ps.rem)$Result)
head(adiv)

#Summarize
adiv %>%
  group_by(Result) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon))

w1<-wilcox.test(Observed ~ Result, data = adiv, exact = FALSE, conf.int = TRUE)
print(w1)


w2<-wilcox.test(Shannon ~ Result, data = adiv, conf.int = TRUE)              
print(w2)


#fail to reject the null hypothesis of no difference in location between groups.


#Alpha diversity boxplots

variable<- get_variable(ps.rem, "Result") %in% c("Anxious", "Non-Anxious")
#sample_data(ps3)$variable <- factor(variable)

fill <- "#4271AE"
line <- "#1F3552"

alpha_meas = c ("Shannon")
q<- plot_richness(ps.rem, x="variable", measures=alpha_meas)
qa<-q+ geom_boxplot(fill=fill,color=line,alpha=6)+facet_wrap(~samdf5$Type, scales = "fixed",ncol=4)+geom_point(size = 3)+ggtitle("Alpha diversity Boxplot of Samples by Anxiety Type")
qb<-qa+geom_jitter()
plot(qb)

#Plot adiv measures
alpha_plot_10<-adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
  ggplot(aes(x = Result, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Result), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none")

plot(alpha_plot_10)
#######################################################################################################################################################################################################################################################
#Beta Diversity

model3<-lm(evals ~ Type , data=k1)
print(model3)
summary(model3) 

out.bc<- ordinate(ps.rem, method = "MDS", distance = "bray")
evals<-out.bc$vectors
beta_plot_1<-plot_ordination(ps.rem,out.bc,color = "Type",shape = "Visit")+geom_point(size = 4)+ggtitle("Beta diversity plot of Samples by Anxiety Type")+stat_ellipse(aes(group = Type), linetype = 2)+theme(aspect.ratio=1)
plot(beta_plot_1)
#beta_plot_2<-plot_ordination(ps.rem,out.bc,color = "V",shape = "Result")+geom_point(size = 4)+ggtitle("Beta diversity plot of Samples by Visit")+stat_ellipse(aes(group = V), linetype = 2)
#plot(beta_plot_2)
#beta_plot_3<-plot_ordination(ps.rem,out.bc,color = "GestationalWeeks_sinceConception",shape = "Result")+geom_point(size = 4)+ggtitle("Beta diversity plot of Samples by Gestational weeks")
#plot(beta_plot_3)

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps.rem, method="unifrac", weighted=F)
ordination = ordinate(ps.rem, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rem, ordination, color="Type",shape = "Visit") + theme(aspect.ratio=1) + geom_point(size = 4)

#####################################################################################################################################


#############################################################################################################
#Using ADONIS function

#Null hypothesis(H0):-No significant differences between alpha diversity of anxious and non-anxious cases(equal to)
#Alternate hypothesis:-Significant differences between alpha diversity of anxious and non-anxious cases(not equal to)
#Significance level:- 0.05
clr_dist_matrix <- phyloseq::distance(ps.rem, method = "euclidean") 
permanova<-vegan::adonis(clr_dist_matrix ~ samdf5$Result)
print(permanova)

#p-value>0.05 , so fail to reject null hypothesis(cannot conclude a significant difference exists)
###############################################################################################################################

#t1<-permanova<-vegan::adonis(clr_dist_matrix ~ samdf5$Visit)
#print(t1)

#table_otu<-otu_table(prune_taxa(taxa_sums(ps.noncontam)>0, ps.noncontam))
#ps11<- prune_taxa(table_otu,seqtab.nochim)
#table_otu_1<-data.frame(table_otu)
#ps10<-phyloseq(tax_table(taxa_1),sample_data(samdf4),otu_table(table_otu,taxa_are_rows = FALSE), phy_tree(fitGTR_1$tree))
# community table from phyloseq object
#otu.table <- otu_table(ps.rem)
#otu.table_1<- otu.table[order(as.numeric(rownames(otu.table))),]
# Transform by sample ordered transposed matrix into data.frame
#df_otu_table <- as.data.frame(otu.table)
# attach metadata
#dataset<- as.data.frame(sample_data(ps.rem))
#dataset_1<- dataset[order(rownames(dataset)),]
# Add Sample type column to ordered OTU dataframe
#df_otu_table$Group <- dataset$Type
# create separate data frames for each sample type
#otu_minimal <- subset(df_otu_table, Group == "Minimal")
#otu_mild <- subset(df_otu_table, Group == "Mild")
#otu_moderate <- subset(df_otu_table, Group == "Moderate")
#otu_severe<- subset(df_otu_table, Group == "Severe")
#minimal_ancom.out <- ANCOM(otu_table(ps.rem),sample_data(ps.rem),sig = 0.05, multcorr = 2, tau = 0.02, theta = 0.1)

###########################################################################################################################################################
#metagenomeSeq is designed to determine features (be it Operational Taxanomic Unit (OTU), species, etc.) that are differentially abundant between two or more groups of multiple samples. metagenomeSeq is designed to address the effects of both normalization and under-sampling of microbial communities on disease association detection and the testing of feature correlations.
set.seed(123)
featureData =data.frame(tax_table(ps.rem))
rownames(featureData)<-paste("Seq",1:nrow(featureData), sep="_")
matrixData<-matrix(otu_table(ps.rem),ncol=ncol(otu_table(ps.rem)))
rownames(matrixData)<-rownames(otu_table(ps.rem))
colnames(matrixData)<-colnames(otu_table(ps.rem))
matrixData<-t(matrixData)
rownames(matrixData)<-paste("Seq",1:nrow(matrixData), sep="_")

metadata<-sample_data(ps.rem)[match(colnames(matrixData),rownames(sample_data(ps.rem))), ]
metadata<-metadata[match(colnames(matrixData),rownames(metadata)),]

otus_metagenomeSeq<-newMRexperiment(matrixData, 
                                    phenoData = AnnotatedDataFrame(metadata), 
                                    featureData = AnnotatedDataFrame(featureData))

otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq,present=5, depth = 5)
p = cumNormStatFast(otus_metagenomeSeq_filter)    #percentile by which to normalize counts
anxietyData = cumNorm(otus_metagenomeSeq, p = p)

anxiety_sample_data<-pData(anxietyData)
anxietyType<-anxiety_sample_data$Type
visitStatus<-anxiety_sample_data$Visit
gadscores<-anxiety_sample_data$SC5
anxietyResult<-anxiety_sample_data$Result
gestational_weeks<-anxiety_sample_data$GestationalWeeks_sinceConception
age<-anxiety_sample_data$Age
mod<-model.matrix(~anxietyResult+visitStatus)
settings<-zigControl(maxit = 100, verbose = TRUE)
fit<-fitZig(obj = anxietyData, mod = mod, control = settings, useCSSoffset = FALSE)

z1<-median(calculateEffectiveSamples(fit))
mean_eff_size = round(z1)

otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq, present = mean_eff_size, 
                                       depth = 5) 
p = cumNormStatFast(otus_metagenomeSeq_filter)    #percentile by which to normalize counts

anxietyData = cumNorm(otus_metagenomeSeq_filter, p = p)
anxiety_sample_data<-pData(anxietyData)
visitStatus<-anxiety_sample_data$Visit
anxietyType<-anxiety_sample_data$Type
anxietyResult<-anxiety_sample_data$Result
gadscores<-anxiety_sample_data$SC5
anxiety_study<-anxiety_sample_data$Study
#gestational_weeks<-as.factor(anxiety_sample_data$GestationalWeeks_sinceConception)
mod<-model.matrix(~anxietyResult+visitStatus)

settings<-zigControl(maxit = 100, verbose = TRUE)
fit<-fitZig(obj = anxietyData, mod = mod, control = settings, useCSSoffset = FALSE)

zigFit<-slot(fit,"fit")
finalMod<-slot(fit,"fit")$design

#contrast.matrix<-makeContrasts(anxietyStatusMinimal-Intercept,anxietyStatusModerate-Intercept,anxietyStatusSevere-Intercept,visitStatusV2-Intercept, visitStatusV4-Intercept,levels = finalMod)
#fit1<-contrasts.fit(zigFit,contrast.matrix)

fit1<-eBayes(zigFit, trend = TRUE, robust = TRUE)
topTable(fit1)

#res<-decideTests(fit1, method = "separate", adjust.method = "BH", p.value = 0.05)
#summary(res)

pvalues<-apply(fit1$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG<-cbind(fit1$coefficients,pvalues)
otus_mG<-data.frame(otus_mG)

otus_mG_filter_anxiety<-otus_mG[(otus_mG$anxietyResultANC.1<=0.05),]
otus_mG_filter_anxiety<-otus_mG_filter_anxiety[(otus_mG_filter_anxiety$visitStatusV2.1>=0.05),]
otus_mG_filter_anxiety<-otus_mG_filter_anxiety[(otus_mG_filter_anxiety$visitStatusV4.1>=0.05),]

h<-rownames(otus_mG_filter_anxiety)
#h1<-rownames(otus_mG_filter_v2)
#h2<-rownames(otus_mG_filter_v4)

m<-data.frame(matrixData)
n1<-data.frame(featureData)

#n2<-n[c(h),]

n3<-n1[c(h),]
#n4<-n1[c(h1),]
#n5<-n1[c(h2),]

m1<-m[c(h),]

#n4<-n3["Family"]
#n5<-data.frame(samdf5)
#n6<-filter(n5,Result =="Anxious")
#n7<-filter(n5, Result =="Non-Anxious")
#n8<-subset(n2, select = -c(S134,S145,S174,S261,S299,S427,S440,S93)) #for non-anxious
#n9<-subset(n2, select = c(S134,S145,S174,S261,S299,S427,S440,S93)) #for anxious
#n10<-apply(n8,1, function(x){x/sum(x)})
#n11<-t(n10)
#n12<-apply(n9,1,function(x){x/sum(x)})
#n13<-t(n12)
#rownames(n11)<-n4$Family
#rownames(n13)<-n4$Family
#n14<-data.frame(n11)
###########################################################################################################
#Using gad-7 scores

mod1<-model.matrix(~gadscores)
#colnames(mod1)<-c("Interaction","ANXIETY_TYPE", "VISITSTATUS")
settings1<-zigControl(maxit = 100, verbose = TRUE)
fit2<-fitZig(obj = anxietyData, mod = mod1, control = settings1, useCSSoffset = FALSE)

z2<-median(calculateEffectiveSamples(fit2))
mean_eff_size1 = round(z2)

otus_metagenomeSeq_filter1= filterData(otus_metagenomeSeq, present = mean_eff_size1, 
                                       depth = 5) 
p1 = cumNormStatFast(otus_metagenomeSeq_filter1)    #percentile by which to normalize counts
anxietyData1 = cumNorm(otus_metagenomeSeq_filter1, p = p1)

anxiety_sample_data1<-pData(anxietyData1)
visitStatus1<-anxiety_sample_data1$Visit
gadscores1<-anxiety_sample_data1$SC5

mod1<-model.matrix(~gadscores1)
settings1<-zigControl(maxit = 100, verbose = TRUE)
fit2<-fitZig(obj = anxietyData1, mod = mod1, control = settings1, useCSSoffset = FALSE)
zigFit1<-slot(fit2,"fit")
finalMod1<-slot(fit2,"fit")$design

fit3<-eBayes(zigFit1, trend = TRUE, robust = TRUE)
topTable(fit3)

pvalues1<-apply(fit3$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG1<-cbind(fit3$coefficients,pvalues1)
#res1<-decideTests(fit3, method = "separate", adjust.method = "BH", p.value = 0.05)
#summary(res1)
otus_mG1<-data.frame(otus_mG1)

otus_mG1_filter_gadscore<-otus_mG1[(otus_mG1$gadscores1.1<=0.05),]
#otus_mG1_filter_visit2<-otus_mG1_filter_gadscore[(otus_mG1_filter_gadscore$visitStatus1V2.1>=0.05),]
#otus_mG1_filter_visit4<-otus_mG1_filter_visit2[(otus_mG1_filter_visit2$visitStatus1V4.1>=0.05),]

#otus_mG1_filter_v2<-otus_mG1[(otus_mG1$visitStatus1V2.1<=0.05),]

#otus_mG1_filter_v4<-otus_mG1[(otus_mG1$visitStatus1V4.1<=0.05),]

h3<-rownames(otus_mG1_filter_gadscore)
#h4<-rownames(otus_mG1_filter_moderate)
#h5<-rownames(otus_mG1_filter_severe)

#h6<-rownames(otus_mG1_filter_v2)
#h7<-rownames(otus_mG1_filter_v4)

n6<-n1[c(h3),]
#n7<-n1[c(h4),]
#n8<-n1[c(h5),]

#n9<-n1[c(h6),]
#n10<-n1[c(h7),]
##############################################################################################################################################
# By anxiety Type

#mod2<-model.matrix(~anxietyType)
#settings2<-zigControl(maxit = 100, verbose = TRUE)
#fit4<-fitZig(obj = anxietyData, mod = mod2, control = settings1, useCSSoffset = FALSE)

#z3<-median(calculateEffectiveSamples(fit4))
#mean_eff_size2 = round(z3)

#otus_metagenomeSeq_filter2= filterData(otus_metagenomeSeq, present = mean_eff_size2, depth = 5) 
#p2 = cumNormStatFast(otus_metagenomeSeq_filter2)    #percentile by which to normalize counts
#anxietyData2 = cumNorm(otus_metagenomeSeq_filter2, p = p2)

#anxiety_sample_data2<-pData(anxietyData2)

#anxietyType2<-anxiety_sample_data2$Type

#mod2<-model.matrix(~anxietyType2)
#settings2<-zigControl(maxit = 100, verbose = TRUE)
#fit4<-fitZig(obj = anxietyData2, mod = mod2, control = settings2, useCSSoffset = FALSE)
#zigFit2<-slot(fit4,"fit")
#finalMod2<-slot(fit4,"fit")$design

#fit5<-eBayes(zigFit2, trend = TRUE, robust = TRUE)
#topTable(fit5)

#pvalues2<-apply(fit5$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
#otus_mG2<-cbind(fit5$coefficients,pvalues2)
#res1<-decideTests(fit3, method = "separate", adjust.method = "BH", p.value = 0.05)
#summary(res1)
#otus_mG2<-data.frame(otus_mG2)

#otus_mG2_filter_minimal<-otus_mG2[(otus_mG2$anxietyType2Minimal.1<=0.05),]
#otus_mG2_filter_moderate<-otus_mG2[(otus_mG2$anxietyType2Moderate.1<=0.05),]
#otus_mG2_filter_severe<-otus_mG2[(otus_mG2$anxietyType2Severe.1<=0.05),]

#h8<-rownames(otus_mG2_filter_minimal)
#h9<-rownames(otus_mG2_filter_moderate)
#h10<-rownames(otus_mG2_filter_severe)

#n11<-n1[c(h8),]
#n12<-n1[c(h9),]
#n13<-n1[c(h10),]

###############################################################################################################################
#By  1st visit only

metadata_V1<-metadata[which(metadata$Visit=="V1"),]
matrixData_V1<-matrixData[,match(rownames(metadata_V1),colnames(matrixData))]

otus_metagenomeSeq_V1<-newMRexperiment(matrixData_V1, 
                                       phenoData = AnnotatedDataFrame(metadata_V1), 
                                       featureData = AnnotatedDataFrame(featureData))

otus_metagenomeSeq_filter_V1 = filterData(otus_metagenomeSeq_V1,present=5, depth = 5)
p = cumNormStatFast(otus_metagenomeSeq_filter_V1)    #percentile by which to normalize counts
anxietyData_V1 = cumNorm(otus_metagenomeSeq_filter_V1, p = p)

anxiety_sample_data_V1<-pData(anxietyData_V1)

anxietyResult_V1<-anxiety_sample_data_V1$Result
visitStatus_V1<-anxiety_sample_data_V1$Visit
gadscores_V1<-anxiety_sample_data_V1$SC5
age_V1<-anxiety_sample_data_V1$Age

mod3<-model.matrix(~gadscores_V1)
settings3<-zigControl(maxit = 100, verbose = TRUE)
fit6<-fitZig(obj = anxietyData_V1, mod = mod3, control = settings3, useCSSoffset = FALSE)

z4<-median(calculateEffectiveSamples(fit6))
mean_eff_size3 = round(z4)

otus_metagenomeSeq_filter_V1= filterData(otus_metagenomeSeq_V1, present = mean_eff_size3, 
                                       depth = 5) 
p3 = cumNormStatFast(otus_metagenomeSeq_filter_V1)    #percentile by which to normalize counts
anxietyData_V1 = cumNorm(otus_metagenomeSeq_filter_V1, p = p3)

anxiety_sample_data_V1<-pData(anxietyData_V1)

anxietyResult_V1<-anxiety_sample_data_V1$Result
visitStatus_V1<-anxiety_sample_data_V1$Visit
gadscores_V1<-anxiety_sample_data_V1$SC5
age_V1<-anxiety_sample_data_V1$Age

mod3<-model.matrix(~gadscores_V1)
settings3<-zigControl(maxit = 100, verbose = TRUE)
fit6<-fitZig(obj = anxietyData_V1, mod = mod3, control = settings3, useCSSoffset = FALSE)
zigFit3<-slot(fit6,"fit")
finalMod3<-slot(fit6,"fit")$design

fit7<-eBayes(zigFit3, trend = TRUE, robust = TRUE)
topTable(fit7)

pvalues3<-apply(fit7$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG3<-cbind(fit7$coefficients,pvalues3)
#res1<-decideTests(fit3, method = "separate", adjust.method = "BH", p.value = 0.05)
#summary(res1)
otus_mG3<-data.frame(otus_mG3)

otus_mG3_filter_V1<-otus_mG3[(otus_mG3$gadscores_V1.1<=0.05),]
#otus_mG3_filter_v4<-otus_mG3[(otus_mG3$visitStatus2V4.1<=0.05),]
#otus_mG3_filter_V1<-otus_mG3_filter_V1[(otus_mG3_filter_V1$age_V1.1>=0.05),]


h11<-rownames(otus_mG3_filter_V1)
#h12<-rownames(otus_mG3_filter_v4)

n14<-n1[c(h11),]
#n15<-n1[c(h12),]
######################################################################################################################################################################################
#By visit2 only

metadata_V2<-metadata[which(metadata$Visit=="V2"),]
matrixData_V2<-matrixData[,match(rownames(metadata_V2),colnames(matrixData))]

otus_metagenomeSeq_V2<-newMRexperiment(matrixData_V2, 
                                       phenoData = AnnotatedDataFrame(metadata_V2), 
                                       featureData = AnnotatedDataFrame(featureData))

otus_metagenomeSeq_filter_V2 = filterData(otus_metagenomeSeq_V2,present=5, depth = 5)
p = cumNormStatFast(otus_metagenomeSeq_filter_V2)    #percentile by which to normalize counts
anxietyData_V2 = cumNorm(otus_metagenomeSeq_filter_V2, p = p)

anxiety_sample_data_V2<-pData(anxietyData_V2)

anxietyResult_V2<-anxiety_sample_data_V2$Result
visitStatus_V2<-anxiety_sample_data_V2$Visit
gadscores_V2<-anxiety_sample_data_V2$SC5
age_V2<-anxiety_sample_data_V2$Age

mod4<-model.matrix(~gadscores_V2)
settings4<-zigControl(maxit = 100, verbose = TRUE)
fit8<-fitZig(obj = anxietyData_V2, mod = mod4, control = settings4, useCSSoffset = FALSE)

z5<-median(calculateEffectiveSamples(fit8))
mean_eff_size4 = round(z5)

otus_metagenomeSeq_filter_V2= filterData(otus_metagenomeSeq_V2, present = mean_eff_size4, 
                                         depth = 5) 
p4 = cumNormStatFast(otus_metagenomeSeq_filter_V2)    #percentile by which to normalize counts
anxietyData_V2 = cumNorm(otus_metagenomeSeq_filter_V2, p = p4)

anxiety_sample_data_V2<-pData(anxietyData_V2)

anxietyResult_V2<-anxiety_sample_data_V2$Result
visitStatus_V2<-anxiety_sample_data_V2$Visit
gadscores_V2<-anxiety_sample_data_V2$SC5
age_V2<-anxiety_sample_data_V2$Age

mod4<-model.matrix(~gadscores_V2)
settings4<-zigControl(maxit = 100, verbose = TRUE)
fit8<-fitZig(obj = anxietyData_V2, mod = mod4, control = settings3, useCSSoffset = FALSE)
zigFit4<-slot(fit8,"fit")
finalMod4<-slot(fit8,"fit")$design

fit9<-eBayes(zigFit4, trend = TRUE, robust = TRUE)
topTable(fit9)

pvalues4<-apply(fit9$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG4<-cbind(fit9$coefficients,pvalues4)
#res1<-decideTests(fit3, method = "separate", adjust.method = "BH", p.value = 0.05)
#summary(res1)
otus_mG4<-data.frame(otus_mG4)

otus_mG4_filter_V2<-otus_mG4[(otus_mG4$gadscores_V2.1<=0.05),]
#otus_mG3_filter_v4<-otus_mG3[(otus_mG3$visitStatus2V4.1<=0.05),]
#otus_mG3_filter_V1<-otus_mG3_filter_V1[(otus_mG3_filter_V1$age_V1.1>=0.05),]


h12<-rownames(otus_mG4_filter_V2)
#h12<-rownames(otus_mG3_filter_v4)

n15<-n1[c(h12),]
#n15<-n1[c(h12),]


#################################################################################################################################################################################################################################################

#By visit4 only


metadata_V4<-metadata[which(metadata$Visit=="V4"),]
matrixData_V4<-matrixData[,match(rownames(metadata_V4),colnames(matrixData))]

otus_metagenomeSeq_V4<-newMRexperiment(matrixData_V4, 
                                       phenoData = AnnotatedDataFrame(metadata_V4), 
                                       featureData = AnnotatedDataFrame(featureData))

otus_metagenomeSeq_filter_V4 = filterData(otus_metagenomeSeq_V4,present=5, depth = 5)
p = cumNormStatFast(otus_metagenomeSeq_filter_V4)    #percentile by which to normalize counts
anxietyData_V4 = cumNorm(otus_metagenomeSeq_filter_V4, p = p)

anxiety_sample_data_V4<-pData(anxietyData_V4)

anxietyResult_V4<-anxiety_sample_data_V4$Result
visitStatus_V4<-anxiety_sample_data_V4$Visit
gadscores_V4<-anxiety_sample_data_V4$SC5
age_V4<-anxiety_sample_data_V4$Age


mod5<-model.matrix(~gadscores_V4)
settings5<-zigControl(maxit = 100, verbose = TRUE)
fit10<-fitZig(obj = anxietyData_V4, mod = mod5, control = settings5, useCSSoffset = FALSE)


z6<-median(calculateEffectiveSamples(fit10))
mean_eff_size5 = round(z6)


otus_metagenomeSeq_filter_V4= filterData(otus_metagenomeSeq_V4, present = mean_eff_size5, 
                                         depth = 5) 
p5 = cumNormStatFast(otus_metagenomeSeq_filter_V4)    #percentile by which to normalize counts
anxietyData_V4 = cumNorm(otus_metagenomeSeq_filter_V2, p = p5)

anxiety_sample_data_V4<-pData(anxietyData_V4)

anxietyResult_V4<-anxiety_sample_data_V4$Result
visitStatus_V4<-anxiety_sample_data_V4$Visit
gadscores_V4<-anxiety_sample_data_V4$SC5
age_V4<-anxiety_sample_data_V4$Age

mod5<-model.matrix(~gadscores_V4)
settings5<-zigControl(maxit = 100, verbose = TRUE)
fit10<-fitZig(obj = anxietyData_V4, mod = mod5, control = settings5, useCSSoffset = FALSE)
zigFit5<-slot(fit10,"fit")
finalMod5<-slot(fit10,"fit")$design

fit11<-eBayes(zigFit5, trend = TRUE, robust = TRUE)
topTable(fit11)

pvalues5<-apply(fit11$p.value,2,function(x){p.adjust(as.numeric(x),method="fdr")})
otus_mG5<-cbind(fit11$coefficients,pvalues5)
#res1<-decideTests(fit3, method = "separate", adjust.method = "BH", p.value = 0.05)
#summary(res1)
otus_mG5<-data.frame(otus_mG5)

otus_mG5_filter_V4<-otus_mG5[(otus_mG5$gadscores_V4.1<=0.05),]


h13<-rownames(otus_mG5_filter_V4)

n16<-n1[c(h13),]
#########################################################################################################################################################################3
