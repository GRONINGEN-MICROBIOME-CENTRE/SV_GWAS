library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)
library(ape)
library(microbiome)
library(corrplot)
library(Hmisc)
library(data.table)

library(ggmosaic)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(ggord)
library(gplots)
library(cowplot)

## color panels
mycolor2_blue_yellow <- c("#08B2E3","#FED766")
mycolor2_green_blue  <- c("#2EC4B6","#235789")
mycolor2_blue_red    <- c("#08B2E3","#ee6352")
mycolor2_green_blue <- c("#2EC4B6","#235789")

prepData <- function(inDf,getLevel="S",isHMP2=F,minBac=95,doDropZero=T,doRescale=T) {
  if (isHMP2) {
    taxRaw <- inDf
  } else {
    taxRaw <- as.data.frame(t(inDf))
  }
  if (doDropZero) {
    taxF <- taxRaw[rowSums(taxRaw) > 0,]
  } else {
    taxF <- taxRaw
  }
  
  #  taxF <- taxF[taxRaw$k__Bacteria >= minBac,]
  # take only bacteria
  taxF <- taxF[,grep('k__Bacteria',colnames(taxF))]
  # take only species
  taxF <- taxF[,grep('t__',colnames(taxF),invert = T)]
  taxF.spec <- taxF[,grep('s__',colnames(taxF))]
  # take only genera
  taxF <- taxF[,grep('s__',colnames(taxF),invert = T)]
  taxF.gen <- taxF[,grep('g__',colnames(taxF))]
  # take only families
  taxF <- taxF[,grep('g__',colnames(taxF),invert = T)]
  taxF.family <- taxF[,grep('f__',colnames(taxF))]
  # take only orders
  taxF <- taxF[,grep('f__',colnames(taxF),invert = T)]
  taxF.order <- taxF[,grep('o__',colnames(taxF))]
  # take only classes
  taxF <- taxF[,grep('o__',colnames(taxF),invert = T)]
  taxF.class <- taxF[,grep('c__',colnames(taxF))]
  # take only phyla
  taxF <- taxF[,grep('c__',colnames(taxF),invert = T)]
  taxF.phyla <- taxF[,grep('p__',colnames(taxF))]
  # rescale
  if (getLevel=="S") {taxF.final <- taxF.spec
  } else if (getLevel=="G") {taxF.final <- taxF.gen
  } else if (getLevel=="F") {taxF.final <- taxF.family
  } else if (getLevel=="O") {taxF.final <- taxF.order
  } else if (getLevel=="C") {taxF.final <- taxF.class
  } else if (getLevel=="P") {taxF.final <- taxF.phyla}
  if (doRescale) {
    taxF.final.rescaled <- taxF.final/rowSums(taxF.final)
  } else {
    taxF.final.rescaled <- taxF.final
  }
  # return
  colnames(taxF.final.rescaled) <- gsub('\\|','.',colnames(taxF.final.rescaled))
  taxF.final.rescaled
}

## shorten the metaphlan2 taxa name
shortenNames <- function(inDF,
                         sep = "\\|",
                         direction = 2) {
  if (direction == 2) {
    for (c in grep('__', rownames(inDF))) {
      cn <- rownames(inDF)[c]
      cns <- unlist(strsplit(cn, sep))
      cnsf <- cns[length(cns)]
      rownames(inDF)[c] <- cnsf
    }
    
    
  } else{
    for (c in grep('__', colnames(inDF))) {
      cn <- colnames(inDF)[c]
      cns <- unlist(strsplit(cn, sep))
      cnsf <- cns[length(cns)]
      colnames(inDF)[c] <- cnsf
    }
  }
  inDF
}

## The Normal Quantile Transformation
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}



## linear model
lm_btw_mats<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #y_mat<- y_mat_i #lld_intri[,c(1:3)]
  #x_mat<- x_mat #lld_vsv[,c(1:3)]
  #cov_mat<- cov_mat# lld_covar
  #covar<- covar_i #covar
  ## test block
  
  my_lm<-function(y,x){
    # y<-y_mat[,1]
    #  x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    try(lm_res <- summary(lm(Y~.,data = lm_input)), silent = T)
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 10, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}



##  mixed logistic model
lmm_btw_mats<-function(y_mat,x_mat,cov_mat,covar,random_factor,convert.x = T){
  require(reshape2)
  require(R.utils)
  require(lme4)
  library(lmerTest)
  
  ## test block
  #  y_mat<- hmp2_all_vsv
  #  x_mat<- hmp2_all_phen[,c("Ctrl_IBD", "Ctrl_UC", "Ctrl_CD", "UC_CD")]
  #  cov_mat<- hmp2_all_phen
  #  covar<-c("reads_filtered", "consent_age")
  #  random_factor<-c("Participant.ID","site_name")
  ## test block
  
  my_lm<-function(y,x){
    #    x<-x_mat[,1]
    #    y<-y_mat[,350]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar],random = cov_mat[,random_factor]) %>% na.omit
    #%>% sapply(as.numeric) 
    lm_input$X<-as.numeric(lm_input$X)
    lm_input$Y<-as.numeric(lm_input$Y)
    #    lm_input$random<-as.character(lm_input$random)
    for (i in grep("random", colnames(lm_input))) {
      lm_input[,i]<-as.factor(lm_input[,i])
    }
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    lm_input<-as.data.frame(lm_input)
    lm_input$Y<-qtrans(lm_input$Y)
    if(convert.x){
      lm_input$X<-qtrans(lm_input$X)
    }
    
    for (i in covar) {
      lm_input[,grep(covar[1], colnames(lm_input))]<-qtrans(lm_input[,grep(covar[1], colnames(lm_input))])
    }
    
    cov_formula<-colnames(lm_input)[-grep( "^Y|^X|^random",colnames(lm_input)) ]%>% paste(collapse = "+")
    rand_formula<-lm_input%>%colnames(.) %>%.[grep("random",.)]%>%paste("(1|",.,")",sep = "")%>%paste(collapse = "+")
    lmm_formula<-paste("Y~X",cov_formula,rand_formula,sep = "+")
    
    try(lmm_res <- lmer(lmm_formula,data = lm_input), silent = T)
    try(lmm_summ <- summary(lmm_res), silent = T)
    
    indv<-'X'
    
    try(beta    <- lmm_summ$coefficients[match(indv,rownames(lmm_summ$coefficients)),1],silent = T)
    try(se      <- lmm_summ$coefficients[match(indv,rownames(lmm_summ$coefficients)),2], silent = T)
    try(p.value <- lmm_summ$coefficients[match(indv,rownames(lmm_summ$coefficients)),5],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 10, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}


## correlation
cor_btw_mats<-function(y_mat,x_mat){
  require(reshape2)
  require(R.utils)
  
  ## test block
  # y_mat<- taxa.clr.filtered.dag3 #lld_intri[,c(1:3)]
  # x_mat<- gene.dag3 #lld_vsv[,c(1:3)]
  ## test block
  
  my_cor<-function(y,x){
    # y<-y_mat[,1]
    # x<-x_mat[,1]
    
    R    <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input <- data.frame(Y = y, X = x) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    # lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    cor_res <- cor.test(x, y, method = "spearman")
    
    R    <- cor_res$estimate
    p.value <- cor_res$p.value
    
    
    try(return(list(R = R,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_cor(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 9, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,2], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,2], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "R", "P","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}

## pie chart
my_pie<-function(in1,item,mycol = c("#F24D4D","#4D7EAE")){
  require(tibble)
  require(showtext)
  require(Cairo)
  require(ggsci)
  require(tibble)
  require(scales)
  require(ggrepel)
  require(forcats)
  require(scatterpie)
  
  showtext_auto()
  
  in1.tibble            <- as_tibble(subset(in1, in1$items%in%item))
  in1.tibble$categories <- fct_reorder(in1.tibble$categories, in1.tibble$value)
  in1.tibble            <- in1.tibble[order(in1.tibble$value, decreasing = TRUE), ]
  piepercent            <- round(100*in1.tibble$value/sum(in1.tibble$value), 2)
  my_labels             <- tibble(x.breaks = seq(1.3, 1.3, length.out = length(piepercent)),
                                  y.breaks = cumsum(in1.tibble$value) - in1.tibble$value/2,
                                  labels = paste(in1.tibble$categories, "\n",in1.tibble$value,", ",piepercent, "%", sep = ""),
                                  categories = in1.tibble$categories)
  
  pdfName    <- paste(item, ".pie.pdf", sep = "")
  ggplot(in1.tibble, aes(x = 1, y = value, fill = categories)) +
    ggtitle(paste(item)) +
    geom_bar(stat="identity", color='white') + 
    coord_polar(theta='y') + 
    theme(legend.position = "None",
          axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels.
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    #    scale_fill_brewer(palette = "Set3", direction = -1)+
    scale_fill_manual(values=mycol) + # set fill color manually
    geom_text_repel(data = my_labels, 
                    aes(x = x.breaks, y = y.breaks, label = labels),
                    size = 7,        # label text size
                    show.legend = FALSE,
                    inherit.aes = FALSE,
                    arrow = arrow(length=unit(0.01, "npc")),
                    force = 1,
                    segment.color = 'grey50'
    )
}

