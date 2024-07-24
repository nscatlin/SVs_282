#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(gridExtra)
library(tidyverse)
library(Dict)

#################
### FUNCTIONS ###
#################
########################################################################################################################################################################
# Calculate q-values for FDR multiple testing correction
.qvalue = function(p) {
  smooth.df = 3
  if (min(p) < 0 || max(p) > 1) {
    print("ERROR: p-values not in valid range.")
    return(0)
  }
  lambda = seq(0, 0.9, 0.05)
  m = length(p)
  pi0 = rep(0, length(lambda))
  for (i in 1:length(lambda)) {
    pi0[i] = mean(p >= lambda[i])/(1 - lambda[i])
  }
  spi0 = smooth.spline(lambda, pi0, df = smooth.df)
  pi0 = predict(spi0, x = max(lambda))$y
  pi0 = min(pi0, 1)
  if (pi0 <= 0) {
    print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
    return(0)
  }
  u = order(p)
  qvalue.rank = function(x) {
    idx = sort.list(x)
    fc = factor(x)
    nl = length(levels(fc))
    bin = as.integer(fc)
    tbl = tabulate(bin)
    cs = cumsum(tbl)
    tbl = rep(cs, tbl)
    tbl[idx] = tbl
    return(tbl)
  }
  v = qvalue.rank(p)
  qvalue = pi0 * m * p/v
  qvalue[u[m]] = min(qvalue[u[m]], 1)
  for (i in (m - 1):1) {
    qvalue[u[i]] = min(qvalue[u[i]], qvalue[u[i + 1]],
                       1)
  }
  return(qvalue)
}


# Calculates the number of effective markers from genotype correlations

Keff = function(r2,alpha) {
  m = nrow(r2)
  if (m > 1) {
    Q = sqrt(r2)
    Q[upper.tri(Q,diag=T)] = NA
    rmax = apply(Q[-1,],1,max,na.rm=T)
    kappa = sqrt(1-rmax^(-1.31*log10(alpha)))
    return(1+sum(kappa))
  } else {
    return(1)
  }
}

# Sets a threshold for results from GEMMA. You only need to supply a genotype file for the "M.eff" method, which estimates
# the effective number of markers and adjusts the significance threshold by dividing by the effective number of markers.
# "FDR" adjusts the significance threshold by the false discovery rate. "Bonferroni" adjusts the significance threshold
# by Bonferroni correction.


setThreshold = function(gemmaResults, genoBimbam = NULL, level = 0.05, method = c("M.eff", "Bonferroni", "FDR"), pvals = "p_lrt") {
  names(gemmaResults)[which(names(gemmaResults) == pvals)] <- "pval"
  gemmaResults$pval = as.numeric(gemmaResults$pval)
  iv = which(!is.na(-log10(gemmaResults$pval)))
  scores = as.vector(-log10(gemmaResults$pval)[iv])
  m = length(scores)
  if(method == "Bonferroni"){
    threshold = -log10(level/m)
  }
  if(method == "FDR"){
    tmp = cbind(10^(-scores), .qvalue(10^(-scores)))
    tmp = tmp[order(tmp[, 2]), ]
    if (tmp[1, 2] > level) {
      threshold = -log10(tmp[which.min(tmp[which(tmp[,2] == min(tmp[,2])),1])]) * 1.2
    } else{
      k = max(which(tmp[, 2] < level))
      threshold = -log10(mean(tmp[k:(k + 1), 1]))
    }}
  if(method == "M.eff"){
    if(is.null(genoBimbam)){
      stop("Genotype file is required to calculate effective number of markers")
    }
    
    gemmaResults = arrange(gemmaResults, chr, ps)
    
    if("ID" %in% colnames(genoBimbam)){
      genoBimbam = genoBimbam %>%
        rename(rs = ID)
    }
    
    genoOnly = genoBimbam %>%
      filter(rs %in% gemmaResults$rs) %>%
      right_join(gemmaResults[,1:3], by = "rs") %>%
      arrange(chr, ps) %>%
      select(-chr, -ps, -rs, -allele1, -allele0)
    
    chrom = levels(as.factor(gemmaResults$chr))
    n.chrom = length(chrom)
    r2 = vector("list", n.chrom)
    names(r2) = chrom
    for (i in chrom) {
      ix = which(gemmaResults$chr == i)
      genoNamed = t(genoOnly)
      colnames(genoNamed) = gemmaResults$rs
      r2[[i]] = cor(genoNamed[, ix], use = "pairwise.complete.obs")^2
    }
    
    me = 0
    
    for (chr in chrom) {
      if (length(r2[[chr]][,1]) > 1) {
        me = me + Keff(r2 = r2[[chr]], alpha = level)
      }
      else {
        me = me + 1
      }
    }
    
    threshold = -log10(level/me)
  }
  
  print(paste("Threshold: ", round(threshold, 2), sep = ""))
  fittedThreshold = gemmaResults %>%
    filter(!is.na(-log10(gemmaResults$pval))) %>%
    mutate(threshold = threshold, lod = scores)
  return(list(threshold, fittedThreshold))
}


manhattan.gemma = function(fittedThreshold, chr = NULL, threshold = NULL) {
  fittedThreshold = arrange(fittedThreshold, chr, ps)
  fittedThreshold$Chrom = as.factor(fittedThreshold$chr)
  
  if (is.null(chr)) {
    a = tapply(fittedThreshold$ps, fittedThreshold$chr, max)
    n = length(a)
    m = tapply(fittedThreshold$ps, fittedThreshold$chr, length)
    b = c(0,apply(array(1:(n-1)),1,function(k){sum(a[1:k])}))
    x = fittedThreshold$ps + rep(b,times=m)
    fittedThreshold$x = x
    ix = 1:nrow(fittedThreshold)
  }else {
    stopifnot(chr %in% levels(fittedThreshold$Chrom))
    ix = which(as.character(fittedThreshold$Chrom) == as.character(chr))
    x = fittedThreshold$ps[ix]/1e+06
    fittedThreshold = fittedThreshold[ix,]
    fittedThreshold$x = x
  }
  
  
  fittedThreshold$color = factor(ifelse(fittedThreshold$chr%%2 == 1, 1, 0))
  
  
  p = ggplot(data = fittedThreshold, aes(x = x, y = lod, color = color)) +
    ylab(expression(paste("-log"[10], "(p)"))) +
    theme_bw() + theme(text = element_text(size = 15), panel.grid = element_blank()) +
    geom_point() + guides(color = "none") +
    scale_shape(solid = FALSE)
  if (is.null(threshold)) {
    if (is.null(chr)) {
      allchr = factor(unique(fittedThreshold$chr))
      breaks = (tapply(x, fittedThreshold$chr, max) + tapply(x, fittedThreshold$chr, min))/2
      p = p + scale_x_continuous(name = "Chromosome", breaks = breaks,
                                 labels = allchr, expand = c(0,0))+ scale_y_continuous(expand=c(0,0.1)) + scale_colour_manual(values = c("#7BBD00", "#18453B")) +
        geom_hline(yintercept = unique(fittedThreshold$threshold), linetype = 2, colour = "grey50") +
        theme(strip.background = element_rect(linewidth = 1, color = "black"))
    }else {
      p = p + scale_x_continuous(name = "Position (Mb)", expand=c(0,0))+ scale_y_continuous(expand=c(0,0.1)) +
        scale_colour_manual(values = "#18453B") +
        geom_hline(yintercept = unique(fittedThreshold$threshold), linetype = 2, colour = "grey50") +
        theme(strip.background = element_rect(linewidth = 1, color = "black"))
    }}else {
      if (is.null(chr)) {
        allchr = factor(unique(fittedThreshold$chr))
        breaks = (tapply(x, fittedThreshold$chr, max) + tapply(x, fittedThreshold$chr, min))/2
        p = p + scale_x_continuous(name = "Chromosome", breaks = breaks,
                                   labels = allchr, expand = c(0,0))+ scale_y_continuous(expand=c(0,0.1)) + scale_colour_manual(values = c("#7BBD00", "#18453B")) +
          geom_hline(yintercept = threshold, linetype = 2, colour = "grey50") +
          theme(strip.background = element_rect(linewidth = 1, color = "black"))
      }else {
        p = p + scale_x_continuous(name = "Position (Mb)", expand = c(0,0))+ scale_y_continuous(expand=c(0,0.1)) +
          scale_colour_manual(values = "#18453B") +
          geom_hline(yintercept = threshold, linetype = 2, colour = "grey50") +
          theme(strip.background = element_rect(linewidth = 1, color = "black"))
      }
    }
  return(p)
}
########################
### END OF FUNCTIONS ###
########################
########################################################################################################################################################################

# Set working dir as the dir that the file exists in
setwd('.')

# setwd(dirname(getActiveDocumentContext()$path))

pheno_dict <- Dict$new("1" = "Growing Degree Days to Silking",
                       "2" = "Growing Degree Days to Anthesis",
                       "3" = "Growing Degree Days to Anthesis-Silking Interval",
                       "4" = "Days to Silking",
                       "5" = "Days to Anthesis",
                       "6" = "Anthesis-Silking Interval",
                       "7" = "Plant Height",
                       "8" = "Ear Height",
                       "9" = "Difference between Plant Height and Ear Height",
                       "10" = "Ratio of Ear Height to Plant Height",
                       "11" = "Ratio of Plant Height and Days to Anthesis"
                       
)



# args[1]  =  association file from gemma00

df <- read.table("PHENONUM_all_lmm_options_miss0.1.assoc.txt", header = TRUE, colClasses=c("allele1"="character"))

df$p_lrt_adjust_fdr <- p.adjust(df$p_lrt, method='fdr')
sig <- df %>% 
  filter(p_lrt_adjust_fdr <= 0.05)
if(nrow(sig) > 0){
  write.table(sig,file="PHENONUM_all_lmm_options_miss0.1.assoc.txt_sig_results.txt",row.names=FALSE, quote=FALSE)
  df$qvalue <- .qvalue(df$p_lrt)
  if (!file.exists("man_plots")){
    dir.create("man_plots")
  }
  if (!file.exists("qq_plots")){
    dir.create("qq_plots")
  } 
  df_fdr <- as.data.frame(setThreshold(df,method = 'FDR'))
  manhattan.gemma(df_fdr)+theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, vjust = -5, size  = 15, face = "bold"),
          plot.tag = element_text(size = 16),
          axis.title = element_text(size = 15),
          legend.position="NONE",
          axis.text.x = element_text(color = "black", size = 16),
          axis.text.y = element_text(color = "black", size = 16),
          axis.ticks = element_line(color = "black"))
  ggtitle(paste(pheno_dict$get(PHENONUM), "Manhattan Plot", sep = ' '))+
  ggsave("./man_plots/_man_plot.png", df_man_plot, device = "png", width = 10, height = 8)

  
  qq_plot <- function(df){
    df <-  df %>% 
      # colnames(c("SNP", "chr", "pos", "pval"))
      rename("chr"="CHR",
             "pval"="P",
             "pos"="BP")
    # df$chr <- df$CHR
    # df$pval <- df$P
    # df$pos <- df$BP
    
    df_nrow=nrow(df)
    
    exp.pval=(1:df_nrow-0.5)/df_nrow
    
    exp.pval.log=as.data.frame(-log10(exp.pval))
    
    var.pval=df$pval[order(df$pval)]
    
    var.pval.log=as.data.frame(-log10(var.pval))
    
    N=df_nrow
    
    cupper=-log10(qbeta(0.95,1:N,N-1:N+1))
    
    clower=-log10(qbeta(1-0.95,1:N,N-1:N+1))
    
    df2=cbind(exp.pval.log,var.pval.log,cupper,clower)
    
    colnames(df2)=c("expected","var","cup","clow")
    
    g2=ggplot(df2)+geom_point_rast(aes(x=expected,y=var),colour="black",size=0.5)+geom_abline(slope=1, intercept=0,colour="grey")+geom_line(aes(expected,cup),linetype=2)+geom_line(aes(expected,clow),linetype=2)+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+xlab(expression(paste(-log[10],"(expected ", italic(p))))+ylab(expression(paste(-log[10],"(observed ", italic(p),")")))+theme_classic()
    
    # line below is plotting both the man plot and qq plot
    # g=g1+g2+plot_layout(widths=c(2,1))+plot_annotation(tag_levels='a')+theme(plot.tag=element_text(size=18))
    return(g2)
  }
  df_qq_plot <- qq_plot(df_fdr)  +
    ggtitle(paste(pheno_dict$get(PHENONUM), "QQ Plot", sep = ' '))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave("./qq_plots/PHENONUM_qq_plot.png", df_qq_plot, device = "png", width = 10, height = 8)
  
}
