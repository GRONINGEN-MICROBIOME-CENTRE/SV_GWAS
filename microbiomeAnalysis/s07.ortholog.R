source("functions.R")
library(UpSetR)
library(ComplexHeatmap)

uhgg_hits<-read.table("00.rawData/ortholog/Output.tsv",sep = "\t",header = T,check.names = F) %>% .[order(.$Evalue, decreasing = F),]
uhgg_info<-read.table("00.rawData/info/UHGG_286997assemblies.txt",sep = "\t", header = T, check.names = F)
gene_info <- read.table("../Documents/SV_anno/SV_gene.tsv", sep = "\t", header = T) %>% dplyr::filter(Step !="")

uhgg_info<-grep("Faecalicatena lactaris|Bifidobacterium bifidum|Collinsella aerofaciens|Faecalibacterium prausnitzii", uhgg_info$`Genome Taxonomy Database lineage`) %>%
  uhgg_info[.,]

uhgg_hits$Contig <- uhgg_hits$Genome
uhgg_hits$Genome <- str_extract_all(uhgg_hits$Contig, "GUT_.*_") %>% unlist %>% str_replace_all("_$", "")
uhgg_hits.rmdup <- (!duplicated(uhgg_hits[,c("Genome", "Gene")])) %>% uhgg_hits[.,]
uhgg_hits.long<-uhgg_hits.rmdup[,c("Gene", "Genome", "Evalue")]
uhgg_hits_mat<-spread(uhgg_hits.long, Gene, Evalue)

uhgg_hits_mat.anno<-full_join(uhgg_hits_mat, uhgg_info, by = "Genome")
uhgg_hits_mat.anno$Species <- uhgg_hits_mat.anno$`Genome Taxonomy Database lineage` %>% str_replace_all(".*s__", "")
colnames(uhgg_hits_mat.anno)[2:28]<-colnames(uhgg_hits_mat.anno)[2:28] %>% str_replace_all(".*;", "")
uhgg_hits_mat.anno<-uhgg_hits_mat.anno[,c(1,29:51,match(gene_info$GeneID, colnames(uhgg_hits_mat.anno)))]
colnames(uhgg_hits_mat.anno)[25:33]<-colnames(uhgg_hits_mat.anno)[25:33] %>% match(gene_info$GeneID) %>% gene_info$Step[.]

uhgg_hits_mat.anno$Step1_check <- (!is.na(uhgg_hits_mat.anno$Step_1A)) & 
  (!is.na(uhgg_hits_mat.anno$Step_1B)) &
  (!is.na(uhgg_hits_mat.anno$Step_1C)) &
  (!is.na(uhgg_hits_mat.anno$Step_1D))

uhgg_hits_mat.anno$Step2_check <- (!is.na(uhgg_hits_mat.anno$Step_2))
uhgg_hits_mat.anno$Step3_check <- (!is.na(uhgg_hits_mat.anno$Step_3))
uhgg_hits_mat.anno$Step4_check <- (!is.na(uhgg_hits_mat.anno$Step_4))
uhgg_hits_mat.anno$Step5_check <- (!is.na(uhgg_hits_mat.anno$Step_5.1)) | (!is.na(uhgg_hits_mat.anno$Step_5.2))

uhgg_hits_mat.anno$Complete_GalNAc<- uhgg_hits_mat.anno$Step1_check & 
  uhgg_hits_mat.anno$Step2_check &
  uhgg_hits_mat.anno$Step3_check &
  uhgg_hits_mat.anno$Step4_check &
  uhgg_hits_mat.anno$Step5_check

write.table(uhgg_hits_mat.anno, "07.orgholog/uhgg_hits.tsv",sep = "\t",col.names = T, row.names = F, quote = F)

uhgg_hits_mat.filter<-uhgg_hits_mat.anno %>% filter(Complete_GalNAc==T)
write.table(uhgg_hits_mat.filter, "07.orgholog/uhgg_hits.completeGalNAc.tsv",sep = "\t",col.names = T, row.names = F, quote = F)

step_list <- c("Step1", "Step2", "Step3", "Step4", "Step5")
GalNAc_check <- uhgg_hits_mat.anno[,34:38]
colnames(GalNAc_check) <- step_list
Species <- uhgg_hits_mat.anno$Species %>% str_replace_all("_.*", "")

if(!dir.exists("07.ortholog")){dir.create("07.orgholog")}
pdf("07.orgholog/UHGG.Pathway_completeness_upset.pdf", height = 3.5, width = 10)
ComplexUpset::upset(GalNAc_check,rev(step_list),
                    mode='exclusive_intersection',
                    sort_sets=FALSE,
                    width_ratio=0.2,
                    base_annotations=list(
                      'Intersection size'=ComplexUpset::intersection_size(
                        counts=F,mode='exclusive_intersection',
                        mapping=aes(fill=Species)
                      )+
                        scale_fill_manual(values = wesanderson::wes_palette("Darjeeling1", 4, type = "continuous"))+
                        theme_classic()
                    )
)
dev.off()

species_galnac<-table(Species,uhgg_hits_mat.anno$Complete_GalNAc) %>% as.data.frame

complete_prop <- data.frame(items = rep("Completeness", 2),
                           categories = c("TRUE", "FALSE"),
                           value = c(species_galnac[species_galnac$Var2==TRUE  & species_galnac$Species=="Faecalibacterium prausnitzii", 3],
                                     species_galnac[species_galnac$Var2==FALSE & species_galnac$Species=="Faecalibacterium prausnitzii", 3]))

p_complete_prop_fp<-my_pie(complete_prop, "Completeness",mycol = c("#F69100", "grey50"))

complete_prop <- data.frame(items = rep("Completeness", 2),
                            categories = c("TRUE", "FALSE"),
                            value = c(species_galnac[species_galnac$Var2==TRUE  & species_galnac$Species=="Collinsella aerofaciens", 3],
                                      species_galnac[species_galnac$Var2==FALSE & species_galnac$Species=="Collinsella aerofaciens", 3]))

p_complete_prop_ca<-my_pie(complete_prop, "Completeness",mycol = rev(c("#50A45C", "grey50")))

# Merge figures
p_complete_prop<-plot_grid(p_complete_prop_fp,p_complete_prop_ca,
                         rel_widths = c(1, 1),align = 'hv',
                         #labels = c("vSVs-BAs", "dSVs-BAs"),
                         ncol = 2,label_size	= 8,vjust = 0)

pdf("07.orgholog/GalNAc_completeness_proportion_species.pdf", height = 5, width = 10)
print(p_complete_prop)
dev.off()
