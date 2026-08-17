
#### script is designed to visualise the alignments of the ends of contigs in order to identify telomeric and Telomere-Linked-Repeats (TLRs)
suppressMessages(suppressWarnings(library(gggenomes)))

##get the regions of the HTR (+-TIPSIZE2kb) bed file
contigs=read.csv("contig_ends/SAMPLE.TIPSIZE2kb_ends.bed", sep='\t', header=T)
contigs$seq_id = contigs$contig
contigs$length = contigs$end


##get the minimap2 paf file for the links
links=suppressMessages(suppressWarnings(read_links("contig_ends_alignments/SAMPLE.TIPSIZE2kb_ends.nucmer.paf")))

##get the telomeric end positions
telomeres=read.csv("telomeric_repeats/SAMPLE.TIPSIZE2kb_ends.telomeres.bed", sep='\t' , header=T)
telomeres$seq_id = telomeres$contig
telomeres$length= telomeres$end-telomeres$start
telomeres$strand=telomeres$sense

##features being the candidate repeat regions
repeatbed=read.csv("subtelomeric_repeats/SAMPLE.repeat_rep.ends_blast.bed", sep='\t', header=T)
repeatbed$seq_id = repeatbed$contig
repeatbed$length = repeatbed$end-repeatbed$start


feat_list <- list(
  "subtelo" = repeatbed,
  "telo" = telomeres
)

##plot it...
##may need to rearrange the order of the regions manually
ends=suppressMessages(suppressWarnings(print(gggenomes(seqs=contigs, links=subset(links), feat=feat_list)%>%
                                               pick() %>%
                                               sync() %>%
                                               flip()+
                                               geom_link(aes(fill=((map_match/map_length)*100)) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
                                               scale_fill_gradientn(colours=c("grey100","grey75", "grey50"), name ="Identity (%)", labels=c(80,90,100), breaks=c(80,90,100), limits = c(80, 100))+
                                               geom_seq()+
                                               geom_feat(data=feats("subtelo"), color="red", linewidth=1)+
                                               geom_seq_label(hjust = 1.1, vjust = -0.1)+
                                               geom_variant(data=feats("telo"), shape=19, colour="black")+
                                               theme(legend.position = "top")+
                                               coord_cartesian(xlim = c(-10000,TIPSIZE)))))


widthFrac = max(contigs$length+100000) / 10000
#heightFrac = nrow(regionSeqs) / 2 # for default cases
heightFrac = nrow(contigs) / 2 #

ggsave("plotting_Rscripts/SAMPLE.end_alignments.svg", plot = ends, units = "in", height = heightFrac, width = widthFrac, limitsize = FALSE)


##same as above but now using only contig ends with at least one telomere or TLR
##first get a nonredundant list of all the contig ends with telomeres or TLRs
seq_list <- union(repeatbed$seq_id, telomeres$seq_id)
##now filter the contigs list by needing to contain one from the list
contigs_filtered <- contigs[contigs$seq_id %in% seq_list, ]

##now replot using the new contigs list
ends_filtered=suppressMessages(suppressWarnings(print(gggenomes(seqs=contigs_filtered, links=subset(links), feat=feat_list)%>%
                                                        pick() %>%
                                                        sync() %>%
                                                        flip()+
                                                        geom_link(aes(fill=((map_match/map_length)*100)) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
                                                        scale_fill_gradientn(colours=c("grey100","grey75", "grey50"), name ="Identity (%)", labels=c(80,90,100), breaks=c(80,90,100), limits = c(80, 100))+
                                                        geom_seq()+
                                                        geom_feat(data=feats("subtelo"), color="red", linewidth=1)+
                                                        geom_seq_label(hjust = 1.1, vjust = -0.1)+
                                                        geom_variant(data=feats("telo"), shape=19, colour="black")+
                                                        theme(legend.position = "top")+
                                                        coord_cartesian(xlim = c(-10000,TIPSIZE)))))

widthFrac = max(contigs_filtered$length+100000) / 10000
#heightFrac = nrow(regionSeqs) / 2 # for default cases
heightFrac = nrow(contigs_filtered) / 2 #

ggsave("plotting_Rscripts/SAMPLE.end_alignments.filtered.svg", plot = ends_filtered, units = "in", height = heightFrac, width = widthFrac, limitsize = FALSE)



####plotting the position of the repeats at a genome wide level

##just need two files, a bed file of the entire genome and coordinates of the repeats from alignment of a single repeat
genomebed=read.csv("assemblies/SAMPLE.genome.bed", sep='\t', header=T)
genomebed$seq_id = genomebed$contig
genomebed$length = genomebed$end


repeatbedWG=read.csv("subtelomeric_repeats/SAMPLE.repeat_rep.WG_blast.bed", sep='\t', header=T)
repeatbedWG$seq_id = repeatbedWG$contig
repeatbedWG$length= repeatbedWG$end-repeatbedWG$start


telomeresWG=read.csv("telomeric_repeats/SAMPLE.telomeres.bed", sep='\t' , header=T)
telomeresWG$seq_id = telomeresWG$contig
telomeresWG$length= telomeresWG$end-telomeresWG$start
telomeresWG$strand=telomeresWG$sense
telomeresWG2=subset(telomeresWG, length > 25)

length=max(genomebed$length+100000)

feat_list <- list(
  "subtelo" = repeatbedWG,
  "telo" = telomeresWG2
)


WG=suppressMessages(suppressWarnings(print(gggenomes(seqs=genomebed, feats=feat_list)%>%
                                             pick() +
                                             geom_seq()+
                                             geom_seq_label(hjust = 1.1, vjust = -0.1)+
                                             geom_variant(data=feats("subtelo"), color="red", position = position_nudge(y = .2), shape=6)+
                                             geom_variant(data=feats("telo"), color="black", size=1, shape=19)+
                                             xlim(-500000,length))))

WGwidthFrac = max(genomebed$length+100000) / 500000
#heightFrac = nrow(regionSeqs) / 2 # for default cases
WGheightFrac = nrow(genomebed) / 5 #

ggsave("plotting_Rscripts/SAMPLE.whole_genome.svg", plot = WG, units = "in", height = WGheightFrac, width = WGwidthFrac, limitsize = FALSE)
