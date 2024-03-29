library(ggplot2)
library(tidyr)
library(dplyr)

#setwd('/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v3/plots')
d <- read.delim("dSVs_for_manhattan.txt", header = T, as.is = T, sep = "\t", check.names = F)
colnames(d) <- c("bac","SNP", "CHR", "BP", "P")
d$BP <- as.numeric(d$BP)

snps_to_highlight <- read.delim("dSV.snps_to_highlight.txt", header = F, as.is = T, sep = "\t", check.names = F)[,1]

don <- d %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(d, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight 
  mutate( is_highlight=ifelse(SNP %in% snps_to_highlight, "yes", "no"))
  

axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("dSV_manhattan_noannot.png", type = "cairo", width = 15, height = 5, units = 'in', res = 400)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)),  size=1) +
  scale_color_manual(values = rep(c("#495DA0", "#74AFDF"), 22 )) +
  geom_hline(yintercept=10.52288, color = "#EF3B2C", linewidth = 0.3) +
  # custom X axis:
  xlab("Chromosome") +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) , limits = c(0,46)) +
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=1) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(color="lightgrey", size = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
dev.off()