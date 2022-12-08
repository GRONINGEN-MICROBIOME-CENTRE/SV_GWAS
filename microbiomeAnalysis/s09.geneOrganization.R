library(treeio)
library(ggtree)
library(tidyverse)
library(gggenes)
library(ggthemes)
library(ggtext)
library(wesanderson)

tree <- read.tree("00.rawData/geneOrganization/ATCC.nwk")
strain_info <- read.delim('00.rawData/geneOrganization/group.txt', row.names = 1, sep = '\t')
example_genes <- read.table("00.rawData/geneOrganization/example_genes_577_579_direct.txt",fill=TRUE, header=TRUE,sep="\t")

group <- split(row.names(strain_info), strain_info$GalNAc_utilization)
SV_status <- split(row.names(strain_info), strain_info$SV_status)
example_genes$gene <- factor(example_genes$gene, levels = c("dinB","lacC","kbaZ","nagA","agaS","agaF", "kbaY","agaV","agaC", "agaD","pgmB","Unknown/others"))
example_genes.color <- c("grey50", wes_palette("Zissou1", 10, type = "continuous"),"grey80") 

p <- ggtree(groupOTU(tree,group), branch.length='none') + 
  geom_tiplab(aes(label = paste0("italic('  F. prausnitzii ", label, "')")), 
              parse = TRUE, size = 3.5) +
  geom_tippoint(aes(color = group),size = 2.5, show.legend = T) + 
  xlim_tree(9) + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),color = "white",
             data = example_genes, geom = geom_motif, panel = 'Alignment',
             on = 'dinB', align = 'centre') +
  scale_fill_manual("Genes", values = example_genes.color)+
  scale_color_manual("GalNAc utilization", values =  rev(c('orange','gray')))+
  scale_x_continuous(expand=c(0,4)) +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm'),
        legend.position="right")

facet_widths(p, widths=c(2,3))

if(!dir.exists("09.geneOrganization")){dir.create("09.geneOrganization")}
pdf(file="09.geneOrganization/snpPhyTree.pdf", width=9, height=3)
print(p)
dev.off()
