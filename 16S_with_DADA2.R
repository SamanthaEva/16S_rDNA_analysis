library(ggplot2)
library(phyloseq); packageVersion("phyloseq")
library(ShortRead)
library(dada2)
library(ape) #library for creating  tree

path<-"~/Documents/andy/fastq/" #set path
fns <- list.files(path)
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- paste0(path, fnFs)
fnRs <- paste0(path, fnRs)
plotQualityProfile(fnFs[[1]])
plotQualityProfile(fnFs[[2]])
plotQualityProfile(fnRs[[1]])
plotQualityProfile(fnRs[[2]])
filtFs <- paste0(path, sample.names, "_F_filt.fastq.gz")
filtRs <- paste0(path, sample.names, "_R_filt.fastq.gz")
#trim first 20 nt, filter based on quality, looking for q>30
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    trimLeft=c(20, 20), truncLen=c(240,200), 
                    maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
}

derepFs <- derepFastq(filtFs, verbose=TRUE) #dereplicate
derepRs <- derepFastq(filtRs, verbose=TRUE) #dereplicate
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
#Perform joint sample inference and error rate estimation
dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), selfConsist = TRUE)
dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), selfConsist = TRUE)
#merge paired ends
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#construct sequence table
seqtab <- makeSequenceTable(mergers)
#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
#assign taxonomy, choices are greengenes, rdp, and silva
taxa_gg <- assignTaxonomy(seqtab.nochim, paste0(path, "gg_13_8_train_set_97.fa"))
#please verify the number of columns before naming them, some reference taxonomy datasets only resolve to the Genus level
colnames(taxa_gg) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#transpose sequence table for external analyses
seqtab_nochim_T<-t(seqtab.nochim)
write.csv(seqtab_nochim_T, file = "~/Documents/andy/fastq/seqtab_nochim_T_gg.csv")
write.csv(taxa_gg, file = "~/Documents/andy/fastq/taxa_gg.csv")
write.csv(seqtab.nochim, file = "~/Documents/andy/fastq/seqtab_nochim_gg.csv")
#one way to add sample info is from the sample name.  In our case the sample names are AP-1, AP-2, etc.
#result is not very informative but could be with better sample nomenclature
# samples.out <- rownames(seqtab.nochim)
# sample <- (sapply(strsplit(samples.out, "-"), `[`, 2))
# samdf<-data.frame(Sample=sample)

#A better way to bring in sample names and metadata is with a metadata file.
meta <- read.csv("~/Documents/andy/fastq/meta_no19.csv", row.names=1) #csv with meta data
meta$cell_short<-substr(meta$cell,1,1)#add short version of cell
#invoke phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa_gg))

ps <- phyloseq(otu_table(seqtab_nochim_gg, taxa_are_rows=FALSE), 
               sample_data(meta_no19_sort), 
               tax_table(taxa))
#plot alpha diversity
p=plot_richness(ps, measures=c("Shannon", "Simpson"), color="Location") + theme_bw()
p+ geom_point(size=5, alpha = 0.7) + geom_text(mapping=aes(label=cell), vjust = 2.0)

#make a tree for unifrac
tree=rtree(ntaxa(ps), rooted = TRUE, tip.label = taxa_names(ps))
#add tree to phyloseq object
ps1=merge_phyloseq(ps,tree)
#remove sample 19, which only had 6000 sequences
ps1 <- prune_samples(sample_sums(ps1)>=7000, ps1)
#rarify for beta diversity
eso = rarefy_even_depth(ps1)
#perform ordination
ordu1 = ordinate(ps1, "PCoA", "unifrac", weighted=TRUE)
#p=plot_ordination(ps1, ordu, color = "Location", shape = "seed_or_sample", label="cell", title = "PCoA of unweighted unifrac distance matrix, untrimmed data")
#p+ geom_point(size=5, alpha = 0.7) 

#plot ordination
p=plot_ordination(eso, ordu1, color = "Location", shape = "seed_or_sample",  title = "PCoA of weighted unifrac distance matrix, rarified data")
p + geom_point(size=5, alpha = 0.7) + geom_text(mapping=aes(label=cell), vjust = 2.0)

#p + theme_bw() + geom_text(mapping = aes(label = cell), size = 10, vjust = 1.5) theme(text = element_text(size = 16)) + geom_point(size = 4)


#Previous results were with the raw output from dada2.  Next test is to trim data before ordination.  Not sure I need to do much besides rarify but we'll see.
#Remove OTUs that do not show appear more than 5 times in more than half the samples
wh0 = genefilter_sample(ps1, filterfun_sample(function(x) x > 5), A = 0.5 * nsamples(ps1))



ps2= prune_taxa(wh0, ps1)
#Transform sample counts to relative abundance.
ps2 = transform_sample_counts(ps2, function(x) 1e+06 * x/sum(x))
#ps2 = transform_sample_counts(ps1, function(x) 1e+06 * x/sum(x)) #run this insead of prev line if you don't want to remove rare taxa
#Keep only the most abundant nine phyla.
phylum.sum = tapply(taxa_sums(ps2), tax_table(ps2)[, "Phylum"], sum, na.rm = TRUE) #run to limit phyla
top9phyla = names(sort(phylum.sum, TRUE))[1:9] #run to limit phyla
ps2 = prune_taxa((tax_table(ps2)[, "Phylum"] %in% top9phyla), ps2) #run to limit phyla
#perform ordination and plot results
ordu = ordinate(ps2, "PCoA", "unifrac", weighted=TRUE)
q=plot_ordination(ps2, ordu, color = "Location", shape = "seed_or_sample",  title = "PCoA of weighted unifrac distance matrix, rarified/rare seq removed data")
q + geom_point(size=5, alpha = 0.7) + geom_text(mapping=aes(label=cell), vjust = 2.0)


#code below looks at multiple ordination methods.  Need to make output prettier.
dist = "wunifrac"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist) {
  ordi = ordinate(physeq, method = i, distance = dist)
  plot_ordination(physeq, ordi, "samples", color = "Location", shape = "Date_of_Sample", label = "cell")
}, eso2, dist)
names(plist) <- ord_meths

pdataframe = ldply(plist, function(x) {
  df = x$data[, 1:3]
  colnames(df) = c("Axis_1", "Axis_2", "Axis_3")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

p = ggplot(pdataframe, aes(Axis_1, Axis_2, color = Location, shape = Date_of_Sample, 
                           label = cell))
p = p + geom_point(size = 4)# + geom_polygon()
p = p + facet_wrap(~method, scales = "free")
p = p + scale_fill_brewer(type = "qual", palette = "Set1")
p = p + scale_colour_brewer(type = "qual", palette = "Set1")
p
#############################################################
#gap analysis

pam1 = function(x, k) {
  list(cluster = pam(x, k, cluster.only = TRUE))
}
x = phyloseq:::scores.pcoa(ordu, display = "sites")
# gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn = clusGap(x[, 1:2], FUN = pam1, K.max = 10, B = 100)
gskmn
##############################################################
#look at alpha diversity, MP samples only
ps3 = subset_samples(ps1, Location=="MP")

###############################################################
#adonis
#in excel I re-ordered the meta data to match the unifrac distance matrix
meta_no19_sort <- read.csv("~/meta_no19_sort.csv", row.names=1)eso = rarefy_even_depth(ps1)
meta_no19_sort$cell_short<-substr(meta_no19_sort$cell,1,1)#add short version of cell
ps1=merge_phyloseq(ps,tree)
eso = rarefy_even_depth(ps1)
eso2 = transform_sample_counts(eso, function(x) 1e+06 * x/sum(x))
ufrac<-UniFrac(eso2, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)

Loc_adonis<-adonis(q~meta_no19_sort$Location,permutations=999)
cell_adonis<-adonis(q~meta_no19_sort$cell,permutations=999)
date_adonis<-adonis(q~meta_no19_sort$Date_of_Sample,permutations=999)
cell_s_adonis<-adonis(q~meta_no19_sort$cell_short,permutations=999)

Loc_anosim<-anosim(q,meta_no19_sort$Location,permutations=999)
cell_anosim<-anosim(q,meta_no19_sort$cell,permutations=999)
date_anosim<-anosim(q,meta_no19_sort$Date_of_Sample,permutations=999)
cell_s_anosim<-anosim(q,meta_no19_sort$cell_short,permutations=999)

##############################################################
sample_sum_df <- data.frame(sum = sample_sums(ps))
sample_sum_meta<-merge(sample_sum_df, meta_no19_sort, by=0, all=TRUE)
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
smin <- min(sample_sums(ps))
smean <- mean(sample_sums(ps))
smax <- max(sample_sums(ps))
ggplot(data=sample_sum_meta, aes(x=Sample_Reason_Short, y=sum))+ggtitle("Sample sequencing depth")+
  geom_bar(stat="identity")+coord_flip()+ylab("Read counts")+xlab("Sample")+ theme_grey(base_size = 18) 
################################################################
#let's look at Bacteria only
ps_bac = subset_taxa(ps, Kingdom=="k__Bacteria")
p=plot_richness(ps, measures=c("Shannon", "Observed"), color="Location") + theme_bw()
p+ geom_point(size=5, alpha = 0.7) + geom_text(mapping=aes(label=cell), vjust = 2.0)

#make a tree for unifrac
tree=rtree(ntaxa(ps_bac), rooted = TRUE, tip.label = taxa_names(ps_bac))
#add tree to phyloseq object
ps_bac=merge_phyloseq(ps_bac,tree)
#remove sample 19, which only had 6000 sequences
ps3 <- prune_samples(sample_sums(ps2)>=7000, ps1)
#rarify for beta diversity
eso_bac = rarefy_even_depth(ps_bac)
#perform ordination
ordu = ordinate(eso_bac, "PCoA", "unifrac", weighted=TRUE)
#p=plot_ordination(ps1, ordu, color = "Location", shape = "seed_or_sample", label="cell", title = "PCoA of unweighted unifrac distance matrix, untrimmed data")
#p+ geom_point(size=5, alpha = 0.7) 

#plot ordination
p=plot_ordination(eso_bac, ordu, color = "Location", shape = "seed_or_sample",  title = "PCoA of weighted unifrac distance matrix, bacteria only, rarified data")
p + geom_point(size=5, alpha = 0.7) + geom_text(mapping=aes(label=cell), vjust = 2.0)

#p + theme_bw() + geom_text(mapping = aes(label = cell), size = 10, vjust = 1.5) theme(text = element_text(size = 16)) + geom_point(size = 4)


#Previous results were with the raw output from dada2.  Next test is to trim data before ordination.  Not sure I need to do much besides rarify but we'll see.
#Remove OTUs that do not show appear more than 5 times in more than half the samples
wh0 = genefilter_sample(eso_bac, filterfun_sample(function(x) x > 5), A = 0.5 * nsamples(eso_bac))



ps2= prune_taxa(wh0, eso_bac)
#Transform sample counts to relative abundance.
ps2 = transform_sample_counts(ps2, function(x) 1e+06 * x/sum(x))
#ps2 = transform_sample_counts(ps1, function(x) 1e+06 * x/sum(x)) #run this insead of prev line if you don't want to remove rare taxa
#Keep only the most abundant nine phyla.
phylum.sum = tapply(taxa_sums(ps2), tax_table(ps2)[, "Phylum"], sum, na.rm = TRUE) #run to limit phyla
top9phyla = names(sort(phylum.sum, TRUE))[1:9] #run to limit phyla
ps2 = prune_taxa((tax_table(ps2)[, "Phylum"] %in% top9phyla), ps2) #run to limit phyla
#perform ordination and plot results
ordu = ordinate(ps2, "PCoA", "unifrac", weighted=TRUE)
q=plot_ordination(ps2, ordu, color = "Location", shape = "seed_or_sample",  title = "PCoA of weighted unifrac distance matrix, rarified/rare seq removed data")
q + geom_point(size=5, alpha = 0.7) + geom_text(mapping=aes(label=cell), vjust = 2.0)


#code below looks at multiple ordination methods.  Need to make output prettier.
dist = "wunifrac"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist) {
  ordi = ordinate(physeq, method = i, distance = dist)
  plot_ordination(physeq, ordi, "samples", color = "Location", shape = "Date_of_Sample", label = "cell")
}, eso2, dist)
names(plist) <- ord_meths

pdataframe = ldply(plist, function(x) {
  df = x$data[, 1:3]
  colnames(df) = c("Axis_1", "Axis_2", "Axis_3")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

p = ggplot(pdataframe, aes(Axis_1, Axis_2, color = Location, shape = Date_of_Sample, 
                           label = cell))
p = p + geom_point(size = 4)# + geom_polygon()
p = p + facet_wrap(~method, scales = "free")
p = p + scale_fill_brewer(type = "qual", palette = "Set1")
p = p + scale_colour_brewer(type = "qual", palette = "Set1")
p

eso2 = transform_sample_counts(eso_bac, function(x) 1e+06 * x/sum(x))
ufrac<-UniFrac(eso2, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
q=as.matrix(ufrac)

Loc_adonis<-adonis(q~meta_no19_sort$Location,permutations=999)
cell_adonis<-adonis(q~meta_no19_sort$cell,permutations=999)
date_adonis<-adonis(q~meta_no19_sort$Date_of_Sample,permutations=999)
cell_s_adonis<-adonis(q~meta_no19_sort$cell_short,permutations=999)

Loc_anosim<-anosim(q,meta_no19_sort$Location,permutations=999)
cell_anosim<-anosim(q,meta_no19_sort$cell,permutations=999)
date_anosim<-anosim(q,meta_no19_sort$Date_of_Sample,permutations=999)
cell_s_anosim<-anosim(q,meta_no19_sort$cell_short,permutations=999)

