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

#one way to add sample info is from the sample name.  In our case the sample names are AP-1, AP-2, etc.
#result is not very informative but could be with better sample nomenclature
# samples.out <- rownames(seqtab.nochim)
# sample <- (sapply(strsplit(samples.out, "-"), `[`, 2))
# samdf<-data.frame(Sample=sample)

#A better way to bring in sample names and metadata is with a metadata file.
meta <- read.csv("~/Documents/andy/fastq/meta.csv", row.names=1) #csv with meta data
#invoke phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa_gg))
#plot alpha diversity
p=plot_richness(ps, measures=c("Shannon", "Simpson"), color="Reason_for_Sample") + theme_bw()
p+ geom_point(size=5, alpha = 0.7)

#make a tree for unifrac
tree=rtree(ntaxa(ps), rooted = TRUE, tip.label = taxa_names(ps))
#add tree to phyloseq object
ps1=merge_phyloseq(ps,tree)
#perform ordination
ordu = ordinate(ps1, "PCoA", "unifrac", weighted=TRUE)
#p=plot_ordination(ps1, ordu, color = "Location", shape = "seed_or_sample", label="cell", title = "PCoA of unweighted unifrac distance matrix, untrimmed data")
#p+ geom_point(size=5, alpha = 0.7) 

#plot ordination
p=plot_ordination(ps1, ordu, color = "Location", shape = "seed_or_sample",  title = "PCoA of weighted unifrac distance matrix, untrimmed data")
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
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist) {
  ordi = ordinate(physeq, method = i, distance = dist)
  plot_ordination(physeq, ordi, "samples", color = "Location", shape = "seed_or_sample", label = "cell")
}, ps2, dist)
names(plist) <- ord_meths

pdataframe = ldply(plist, function(x) {
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

p = ggplot(pdataframe, aes(Axis_1, Axis_2, color = Location, shape = seed_or_sample, 
                           label = cell))
p = p + geom_point(size = 4) + geom_polygon()
p = p + facet_wrap(~method, scales = "free")
p = p + scale_fill_brewer(type = "qual", palette = "Set1")
p = p + scale_colour_brewer(type = "qual", palette = "Set1")
p

