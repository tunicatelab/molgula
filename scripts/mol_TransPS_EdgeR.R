
# Libraries ---------------------------------------------------------------

library(edgeR)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(ggrepel)

# Load_count_files --------------------------------------------------------
setwd("/Users/elijahlowe/Desktop/tails/scaffolded_gene_models/")


Mocu.3hpf <- read_tsv("SRR1197965_mocu_GG.xprs") %>% 
  dplyr::rename(Mocu_f3 = eff_counts, mocu_gene = target_id)
Mocu.3hpf #test to see if loaded properly
Mocu.4hpf <- read_tsv("SRR1197522_mocu_GG.xprs") %>% 
  dplyr::rename(Mocu_f4 = eff_counts, mocu_gene = target_id)
Mocu.6hpf <- read_tsv("SRR1197972_mocu_GG.xprs") %>% 
  dplyr::rename(Mocu_f6 = eff_counts, mocu_gene = target_id)

nrow(Mocu.3hpf) == nrow(Mocu.4hpf)

#Now let's combine them all into one dataset

Mocu.data <- Mocu.3hpf %>% 
  inner_join(Mocu.4hpf, by = "mocu_gene") %>% 
  inner_join(Mocu.6hpf, by = "mocu_gene") %>% 
  select(mocu_gene, Mocu_f3, Mocu_f4, Mocu_f6) # %>% 
  # mutate(mocu_gene=gsub("\\|m.+", "", mocu_gene))


# Load M. occulta samples -------------------------------------------------

Mocc.3hpf   <- read_tsv("SRR1197985_mocc_GG.xprs") %>% 
  dplyr::rename(Mocc_f3 = eff_counts, mocc_gene = target_id)
Mocc.3.2hpf <- read_tsv("SRR1197986_mocc_GG.xprs") %>% 
  dplyr::rename(Mocc_f3.2 = eff_counts, mocc_gene = target_id)  
Mocc.4hpf   <- read_tsv("SRR1199464_mocc_GG.xprs") %>% 
  dplyr::rename(Mocc_f4 = eff_counts, mocc_gene = target_id)
Mocc.5hpf   <- read_tsv("SRR1199259_mocc_GG.xprs") %>% 
  dplyr::rename(Mocc_f5 = eff_counts, mocc_gene = target_id)
Mocc.6hpf   <- read_tsv("SRR1199268_mocc_GG.xprs") %>% 
  dplyr::rename(Mocc_f6 = eff_counts, mocc_gene = target_id)

Mocc.data <- Mocc.3hpf %>% 
  inner_join(Mocc.3.2hpf, by = "mocc_gene") %>% 
  inner_join(Mocc.4hpf, by = "mocc_gene") %>% 
  inner_join(Mocc.5hpf, by = "mocc_gene") %>% 
  inner_join(Mocc.6hpf, by = "mocc_gene") %>%
  select(mocc_gene, Mocc_f3, Mocc_f3.2, Mocc_f4, Mocc_f5, Mocc_f6)# %>% 
  #mutate(mocc_gene=gsub("\\|m.+", "", mocc_gene))


# Load hybrid samples -----------------------------------------------------

hyb_ocu.3hpf <- read_tsv("SRR1198337_mocu_GG.xprs") %>% 
  dplyr::rename(hyb_ocu_f3 = eff_counts, mocu_gene = target_id)
hyb_ocu.4hpf <- read_tsv("SRR1198321_mocu_GG.xprs") %>% 
  dplyr::rename(hyb_ocu_f4 = eff_counts, mocu_gene = target_id)
hyb_ocu.6hpf <- read_tsv("SRR1198346_mocu_GG.xprs") %>% 
  dplyr::rename(hyb_ocu_f6 = eff_counts, mocu_gene = target_id)

hyb_ocu.data <- hyb_ocu.3hpf %>% 
  inner_join(hyb_ocu.4hpf, by = "mocu_gene") %>% 
  inner_join(hyb_ocu.6hpf, by = "mocu_gene") %>% 
  select(mocu_gene, hyb_ocu_f3, hyb_ocu_f4, hyb_ocu_f6) #%>% 
 # mutate(mocu_gene=gsub("\\|m.+", "", mocu_gene))

hyb_occ.3hpf <- read_tsv("SRR1198337_mocc_GG.xprs") %>% 
  dplyr::rename(hyb_occ_f3 = eff_counts, mocc_gene = target_id)
hyb_occ.4hpf <- read_tsv("SRR1198321_mocc_GG.xprs") %>% 
  dplyr::rename(hyb_occ_f4 = eff_counts, mocc_gene = target_id)
hyb_occ.6hpf <- read_tsv("SRR1198346_mocc_GG.xprs") %>% 
  dplyr::rename(hyb_occ_f6 = eff_counts, mocc_gene = target_id)

hyb_occ.data <- hyb_occ.3hpf %>% 
  inner_join(hyb_occ.4hpf, by = "mocc_gene") %>% 
  inner_join(hyb_occ.6hpf, by = "mocc_gene") %>% 
  select(mocc_gene, hyb_occ_f3, hyb_occ_f4, hyb_occ_f6)# %>% 
 # mutate(mocc_gene=gsub("\\|m.+", "", mocc_gene))

# Join all datasets, and filter with reciprocal  --------------------------

recip_table <- read_tsv('may_ortho.proteinortho', col_types = list(col_double(),col_double(),col_double(),col_character(),col_character(),col_character()), col_names = c("# Species","Genes", "Alg.-Conn.", "KH", "mocc_gene","mocu_gene")) %>% 
  filter(`# Species` >= 2, !grepl('\\*|,', mocc_gene),!grepl('\\*|,', mocu_gene)) %>% 
  select(-KH, -`# Species`, -Genes, -`Alg.-Conn.`) %>% 
  mutate(mocc_gene=gsub("\\|m.+", "", mocc_gene)) %>% 
  mutate(mocu_gene=gsub("\\|m.+", "", mocu_gene))

#recip <- read_csv('all-mocu_GG_TransPS.x.mocc_GG_TransPS.recip', col_types = list(col_character(),col_character()), col_names = c("mocu_gene","mocc_gene")) #nucleotide orthology
recip <- read_csv('../mocu.x.mocc.TransPS.tblastx.recip', col_types = list(col_character(),col_character()), col_names = c("mocu_gene","mocc_gene")) #translated protein orthology
Mocu.data %>% 
  mutate(mocu_gene=gsub("\\|m.+", "", mocu_gene)) %>% 
  inner_join(recip_table)

all.data <- inner_join(Mocu.data, recip) %>% 
  inner_join(Mocc.data) %>% 
  inner_join(hyb_ocu.data) %>% 
  inner_join(hyb_occ.data) %>% 
  mutate(hyb_f3 = hyb_ocu_f3 + hyb_occ_f3, hyb_f4 = hyb_ocu_f4 + hyb_occ_f4, hyb_f6 = hyb_ocu_f6 + hyb_occ_f6,
         mol_f3 = Mocu_f3 + Mocc_f3.2, mol_f4 = Mocu_f4 + Mocc_f4, mol_f6 = Mocu_f6 + Mocc_f6) %>% 
  column_to_rownames('mocu_gene') %>% 
  select(-mocc_gene)

#####################
# filter <- apply(all.data,1,function(x) mean(x)>10)
# table(filter)
# 
# common <- intersect(names(Mocu.data),
#                     rownames(all.data[filter,]))
# length(common)
# 
# feature <- data.frame(gc=yeastGC,length=yeastLength)
# data <- newSeqExpressionSet(counts=as.matrix(all.data[common,]),
#                             featureData=feature[common,],
#                             phenoData=data.frame(
#                               conditions=c(rep("mut",2),rep("wt",2)),
#                               row.names=colnames(all.data)))
# data
# 
# boxplot(data,col=colors[1:4])

##############
trans_lens <- inner_join(Mocu.3hpf, recip) %>% 
  inner_join(Mocc.3hpf, by = "mocc_gene") %>% 
  select("length.x", "length.y")
plot(x=trans_lens$length.x, y=trans_lens$length.y)
TransPS_fit <- lm(trans_lens$length.x ~ trans_lens$length.y)
summary(TransPS_fit)

# all.data.h <- inner_join(Mocu.data, recip) %>% 
#   inner_join(Mocc.data) %>% 
#   inner_join(hyb_ocu.data) %>% 
#   inner_join(hyb_occ.data) %>% 
#   mutate(hyb_f3 = hyb_ocu_f3 + hyb_occ_f3, hyb_f4 = hyb_ocu_f4 + hyb_occ_f4, hyb_f6 = hyb_ocu_f6 + hyb_occ_f6) %>% 
#   select(-mocc_gene)

#We now have a data frame with all the data in it:
#View(all.data)


# Experimental setup ------------------------------------------------------

Group <- factor(paste(species=c("ocu","ocu","ocu","occ","occ","occ","occ","occ","hyb_ocu","hyb_ocu","hyb_ocu","hyb_occ","hyb_occ","hyb_occ","hyb","hyb","hyb","mol","mol","mol"),
                      stage=c("f3","f4","f6","f3","f3","f4","f5","f6","f3","f4","f6","f3","f4","f6","f3","f4","f6","f3","f4","f6"),sep="."))
edger.data<-DGEList(counts=all.data, group = Group)

# housekeeping <- read.table("all-mocu_GG_TransPS.x.mocc_GG_TransPS.housekeeping")
# 
# edger.data1 <- edger.data
# edger.data1$samples$group <- 1
# edger.data0 <- estimateDisp(edger.data1[housekeeping$V1,], trend="none", tagwise=FALSE)
# edger.data$common.dispersion <- edger.data0$common.dispersion
# fit <- glmFit(edger.data, design)
# lrt <- glmLRT(fit, coef = ocuf4vsocuf6)
# 
keep <- rowSums(cpm(edger.data)>0.5) >=2
table(keep) 
edger.data<-edger.data[keep,, keep.lib.sizes=FALSE]

design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)


ocuf3vsoccf3 = c("occ.f3","ocu.f3")
ocuf4vsoccf4 = c("occ.f4","ocu.f4")
ocuf6vsoccf6 = c("occ.f6","ocu.f6")
hocuf3vshoccf3 = c("hyb_occ.f3","hyb_ocu.f3")
hocuf4vshoccf4 = c("hyb_occ.f4","hyb_ocu.f4")
hocuf6vshoccf6 = c("hyb_occ.f6","hyb_ocu.f6")
ocuf3vshybf3 = c("ocu.f3","hyb.f3")
occf3vshybf3 = c("occ.f3","hyb.f3")
ocuf4vshybf4 = c("ocu.f4","hyb.f4")
occf4vshybf4 = c("occ.f4","hyb.f4")
ocuf6vshybf6 = c("ocu.f6","hyb.f6")
occf6vshybf6 = c("occ.f6","hyb.f6")
ocuf3vsocuf4 = c("ocu.f3","ocu.f4")
occf3vsoccf4 = c("occ.f3","occ.f4")
ocuf4vsocuf6 = c("ocu.f4","ocu.f6")
occf4vsoccf6 = c("occ.f4","occ.f6")
ocuf3vsocuf6 = c("ocu.f3","ocu.f6")
occf3vsoccf6 = c("occ.f3","occ.f6")
hybf3vshybf4 = c("hyb.f3","hyb.f4")
hybf4vshybf6 = c("hyb.f4","hyb.f6")
hybf3vshybf6 = c("hyb.f3","hyb.f6")
molf3vshybf3 = c("mol.f3","hyb.f3")
molf4vshybf4 = c("mol.f4","hyb.f4")
molf6vshybf6 = c("mol.f6","hyb.f6")


#housekeeping <- read_csv("mocu_housekeeping.list", col_names = "housekeeping")
#housekeeping1 <- read_csv("mocu_gimme.x.KH2012_housekeeping.csv", col_names = c("housekeeping","KH_ID"))
#d1 <- edger.data
#d1<-estimateCommonDisp(edger.data)
#d4 <- estimateGLMTagwiseDisp(d3,design)
edger.data <- calcNormFactors(edger.data)
d2 <- estimateGLMCommonDisp(edger.data,design,verbose=TRUE)
d3 <- estimateGLMTagwiseDisp(d2,design)

plotMDS(d3, main="edgeR MDS Plot")

edgeRvol <- function(edgeR_results, plot_title){
  mutateddf <- topTags(edgeR_results, n=nrow(edgeR_results$table))$table %>% 
    rownames_to_column(var = "gene") %>% 
    dplyr::mutate(sig=ifelse(FDR < 0.05 & abs(logFC) > 1.5, "padj<0.05", "Not Sig"))
  volc = ggplot(mutateddf, aes(logCPM,logFC)) + #volcanoplot with log2Foldchange versus pvalue
#  MA <- ggplot(mutateddf, aes(FDR, logFC))  +
    geom_point(aes(col=sig, alpha = sig)) + #add points colored by significance
    scale_color_manual(values=c("black", "red")) +
    scale_alpha_manual(guide='none', values = list('Not Sig' = 0.3, 'padj<0.05' = 0.5)) +
#    xlim(-10.25,10.25) +
    ylim(-10.5, 10.5) +
    ggtitle(plot_title) + #e.g. 'Volcanoplot DESeq2'
    theme_minimal() # simplify plot for publications
  volc+geom_text_repel(data=head(mutateddf, 20), aes(label=gene)) #adding text for the top 20 genes
  #ggsave("Volcanoplot.pdf", device="pdf") #In case you want to easily save to disk
  volc
}

et.mocu.3v4<-exactTest(d3,pair=ocuf3vsocuf4)
summary(de.mocu.3v4 <- decideTestsDGE(et.mocu.3v4, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu.3v4)]
pdf('mocu.3v4.MA.pdf')
edgeRvol(et.mocu.3v4, "M. oculata neurula vs gastrula")
#plotSmear(et.mocu.3v4, ylim=c(-10,10), de.tags=detags, main="M. oculata neurula vs gastrula")
dev.off()

et.mocc.3v4<-exactTest(d3,pair=occf3vsoccf4)
summary(de.mocc.3v4 <- decideTestsDGE(et.mocc.3v4, p=0.05))
detags <- rownames(d3)[as.logical(de.mocc.3v4)]
pdf('mocc.3v4.MA.pdf')
edgeRvol(et.mocc.3v4, "M. occulta neurula vs gastrula")
#plotSmear(et.mocc.3v4, ylim=c(-10,10), de.tags=detags, main="M. occulta neurula vs gastrula")
dev.off()

et.hyb.3v4<-exactTest(d3,pair=hybf3vshybf4)
summary(de.hyb.3v4 <- decideTestsDGE(et.hyb.3v4, p=0.05))
detags <- rownames(d3)[as.logical(de.hyb.3v4)]
pdf('hyb.3v4.MA.pdf')
edgeRvol(et.hyb.3v4, "Hybrid neurula vs gastrula")
#plotSmear(et.hyb.3v4, ylim=c(-10,10), de.tags=detags, main="Hybrid neurula vs gastrula")
dev.off()

et.mocu.4v6<-exactTest(d3,pair=ocuf4vsocuf6)
summary(de.mocu.4v6 <- decideTestsDGE(et.mocu.4v6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu.4v6)]
pdf('mocu.4v6.MA.pdf')
edgeRvol(et.mocu.4v6, "M. oculata tailbud vs neurula")
#plotSmear(et.mocu.4v6, ylim=c(-10,10), de.tags=detags, main="M. oculata tailbud vs neurula")
dev.off()

et.mocc.4v6<-exactTest(d3,pair=occf4vsoccf6)
summary(de.mocc.4v6 <- decideTestsDGE(et.mocc.4v6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocc.4v6)]
pdf('mocc.4v6.MA.pdf')
edgeRvol(et.mocc.4v6, "M. occulta tailbud v neurula")
#plotSmear(et.mocc.4v6, ylim=c(-10,10), de.tags=detags, main="M. occulta tailbud vs neurula")
dev.off()

et.hyb.4v6<-exactTest(d3,pair=hybf4vshybf6)
summary(de.hyb.4v6 <- decideTestsDGE(et.hyb.4v6, p=0.05))
detags <- rownames(d3)[as.logical(de.hyb.4v6)]
pdf('hyb.4v6.MA.pdf')
edgeRvol(et.hyb.4v6, "Hybrid tailbud vs neurula")
#plotSmear(et.hyb.4v6, ylim=c(-10,10), de.tags=detags, main="Hybrid tailbud vs neurula")
dev.off()

et.mocu.3v6<-exactTest(d3,pair=ocuf3vsocuf6)
summary(de.mocu.3v6 <- decideTestsDGE(et.mocu.3v6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu.3v6)]
pdf('mocu.3v6.MA.pdf')
edgeRvol(et.mocu.3v6, "M. oculata tailbud vs gastrula")
#plotSmear(et.mocu.4v6, ylim=c(-10,10), de.tags=detags, main="M. oculata tailbud vs neurula")
dev.off()

et.mocc.3v6<-exactTest(d3,pair=occf3vsoccf6)
summary(de.mocc.3v6 <- decideTestsDGE(et.mocc.3v6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocc.3v6)]
pdf('mocc.3v6.MA.pdf')
edgeRvol(et.mocc.3v6, "M. occulta tailbud v gastrula")
#plotSmear(et.mocc.4v6, ylim=c(-10,10), de.tags=detags, main="M. occulta tailbud vs neurula")
dev.off()

et.hyb.3v6<-exactTest(d3,pair=hybf3vshybf6)
summary(de.hyb.3v6 <- decideTestsDGE(et.hyb.3v6, p=0.05))
detags <- rownames(d3)[as.logical(de.hyb.3v6)]
pdf('hyb.3v6.MA.pdf')
edgeRvol(et.hyb.3v6, "Hybrid tailbud vs gastrula")
#plotSmear(et.hyb.v6, ylim=c(-10,10), de.tags=detags, main="Hybrid tailbud vs neurula")
dev.off()

et.mocuhyb3vmocchyb3<-exactTest(d3,pair=hocuf3vshoccf3)
summary(de.mocuhyb3vmocchyb3 <- decideTestsDGE(et.mocuhyb3vmocchyb3, p=0.05))
detags <- rownames(d3)[as.logical(de.mocuhyb3vmocchyb3)]
edgeRvol(et.mocuhyb3vmocchyb3, "Hybrid 3hpf")
plotSmear(et.mocuhyb3vmocchyb3, ylim=c(-10,10), de.tags=detags, main="Hybrid 3hpf")

et.mocuhyb4vmocchyb4<-exactTest(d3,pair=hocuf4vshoccf4)
summary(de.mocuhyb4vmocchyb4 <- decideTestsDGE(et.mocuhyb4vmocchyb4, p=0.05))
detags <- rownames(d3)[as.logical(de.mocuhyb4vmocchyb4)]
edgeRvol(et.mocuhyb4vmocchyb4, "Hybrid 4hpf")
plotSmear(et.mocuhyb4vmocchyb4, ylim=c(-10,10), de.tags=detags, main="Hybrid 4hpf")

et.mocuhyb6vmocchyb6<-exactTest(d3,pair=hocuf6vshoccf6)
summary(de.mocuhyb6vmocchyb6 <- decideTestsDGE(et.mocuhyb6vmocchyb6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocuhyb6vmocchyb6)]
plotSmear(et.mocuhyb6vmocchyb6, ylim=c(-10,10), de.tags=detags, main="Hybrid 6hpf")

et.mocu3vmocc3<-exactTest(d3,pair=ocuf3vsoccf3)
summary(de.mocu3vmocc3 <- decideTestsDGE(et.mocu3vmocc3, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu3vmocc3)]
plotSmear(et.mocu3vmocc3, ylim=c(-10,10), de.tags=detags, main="M. oculata vs M. occulta 3hpf")

et.mocu4vmocc4<-exactTest(d3,pair=ocuf4vsoccf4)
summary(de.mocu4vmocc4 <- decideTestsDGE(et.mocu4vmocc4, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu4vmocc4)]
plotSmear(et.mocu4vmocc4, ylim=c(-10,10), de.tags=detags, main="M. oculata vs M. occulta 4hpf")

et.mocu6vmocc6<-exactTest(d3,pair=ocuf6vsoccf6)
summary(de.mocu6vmocc6 <- decideTestsDGE(et.mocu6vmocc6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu6vmocc6)]
plotSmear(et.mocu6vmocc6, ylim=c(-10,10), de.tags=detags, main="M. oculata vs M. occulta 6hpf")

et.molf3vshybf3<-exactTest(d3,pair=molf3vshybf3)
summary(de.molf3vshybf3 <- decideTestsDGE(et.molf3vshybf3, p=0.05))
detags <- rownames(d3)[as.logical(de.molf3vshybf3)]
plotSmear(et.molf3vshybf3, ylim=c(-10,10), de.tags=detags, main="M. oculata vs M. occulta 3hpf")

et.molf4vshybf4<-exactTest(d3,pair=molf4vshybf4)
summary(de.molf4vshybf4 <- decideTestsDGE(et.molf4vshybf4, p=0.05))
detags <- rownames(d3)[as.logical(de.molf4vshybf4)]
plotSmear(et.molf4vshybf4, ylim=c(-10,10), de.tags=detags, main="M. oculata vs M. occulta 4hpf")

et.molf6vshybf6<-exactTest(d3,pair=molf6vshybf6)
summary(de.molf6vshybf6 <- decideTestsDGE(et.molf6vshybf6, p=0.05))
detags <- rownames(d3)[as.logical(de.molf6vshybf6)]
plotSmear(et.molf6vshybf6, ylim=c(-10,10), de.tags=detags, main="M. oculata vs M. occulta 6hpf")

# top genes ---------------------------------------------------------------
cond_list <- list(et.mocu.3v4, et.mocu.4v6,et.mocu.3v4,et.mocc.4v6,et.hyb.3v4,et.hyb.4v6,
                  et.mocuhyb3vmocchyb3,et.mocuhyb4vmocchyb4,et.mocuhyb6vmocchyb6)
lapply(cond_list, function(x){
  rownames_to_column(topTags(as.name(x), n=nrow(as.name(x)$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% get(paste0(as.name(x), ".up.csv",sep = ""))
  #write_csv(get(paste0(x,".up.csv",sep = "")))
})

for (i in cond_list) {
  rownames_to_column(topTags(i, n=nrow(i$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% get(paste0(i, ".up.csv",sep = ""))
  write_csv(get(paste0(i,".up.csv",sep = "")))
}

mocu.3v4.up <- rownames_to_column(topTags(et.mocu.3v4, n=nrow(et.mocu.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.3v4.up.csv")
mocu.3v4.down <- rownames_to_column(topTags(et.mocu.3v4, n=nrow(et.mocu.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.3v4.down.csv")
mocc.3v4.up <- rownames_to_column(topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.3v4.up.csv")
mocc.3v4.down <- rownames_to_column(topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.3v4.down.csv")
hyb.3v4.up <- rownames_to_column(topTags(et.hyb.3v4, n=nrow(et.hyb.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.3v4.up.csv")
hyb.3v4.down <- rownames_to_column(topTags(et.hyb.3v4, n=nrow(et.hyb.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.3v4.down.csv")

mocu.4v6.up <- rownames_to_column(topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.4v6.up.csv")
mocu.4v6.down <- rownames_to_column(topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.4v6.down.csv")
mocc.4v6.up <- rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.mocc.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.4v6.up.csv")
mocc.4v6.down <- rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.mocc.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.4v6.down.csv")
hyb.4v6.up <- rownames_to_column(topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.4v6.up.csv")
hyb.4v6.down <- rownames_to_column(topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.4v6.down.csv")

mocu.3v6.up <- rownames_to_column(topTags(et.mocu.3v6, n=nrow(et.mocu.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.3v6.up.csv")
mocu.3v6.down <- rownames_to_column(topTags(et.mocu.3v6, n=nrow(et.mocu.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.3v6.down.csv")
mocc.3v6.up <- rownames_to_column(topTags(et.mocc.3v6, n=nrow(et.mocc.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.3v6.up.csv")
mocc.3v6.down <- rownames_to_column(topTags(et.mocc.3v6, n=nrow(et.mocc.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.3v6.down.csv")
hyb.3v6.up <- rownames_to_column(topTags(et.hyb.3v6, n=nrow(et.hyb.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.3v6.up.csv")
hyb.3v6.down <- rownames_to_column(topTags(et.hyb.3v6, n=nrow(et.hyb.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.3v6.down.csv")

# write all results -------------------------------------------------------


mocc_4v6 <- topTags(et.mocc.4v6, n=nrow(et.mocc.4v6$table))$table %>% 
  rownames_to_column(var = "gene")
hyb_4v6 <- topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb4vmocchyb4 <- topTags(et.mocuhyb4vmocchyb4, n=nrow(et.mocuhyb4vmocchyb4$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb6vmocchyb6 <- topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table %>% 
  rownames_to_column(var = "gene")


all_4v6 <- topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table %>% 
  rownames_to_column(var = "gene") %>%
  inner_join(recip, by = c('gene' ='mocu_gene')) %>%
  inner_join(mocc_4v6, by = 'gene') %>% 
  inner_join(hyb_4v6, by = 'gene') %>%
  inner_join(mocuhyb4vmocchyb4, by = 'gene') %>%
  inner_join(mocuhyb6vmocchyb6, by = 'gene') %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id))  %>% 
  rename(logFC.mocu_tb_nr = logFC.x, logCPM.mocu_tb_nr = logCPM.x, PValue.mocu_tb_nr = PValue.x, FDR.mocu_tb_nr = FDR.x,
         logFC.mocc_tb_nr = logFC.y, logCPM.mocc_tb_nr = logCPM.y, PValue.mocc_tb_nr = PValue.y, FDR.mocc_tb_nr = FDR.y,
         logFC.hyb_tb_nr = logFC.x.x, logCPM.hyb_tb_nr = logCPM.x.x, PValue.hyb_tb_nr = PValue.x.x, FDR.hyb_tb_nr = FDR.x.x,
         logFC.hyb_nr = logFC.y.y, logCPM.hyb_nr = logCPM.y.y, PValue.hyb_nr = PValue.y.y, FDR.hyb_nr = FDR.y.y,
         logFC.hyb_tb = logFC, logCPM.hyb_tb = logCPM, PValue.hyb_tb = PValue, FDR.hyb_tb = FDR) %>% 
  mutate(mocc_gene=gsub("\\@KH.+", "", mocc_gene))

all_4v6 %>% View()
  write_csv("4v6_Mol_GG_edgeR_results.csv")

mocc_3v4 <- topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table %>% 
  rownames_to_column(var = "gene")
hyb_3v4 <- topTags(et.hyb.3v4, n=nrow(et.hyb.3v4$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb3vmocchyb3 <- topTags(et.mocuhyb3vmocchyb3, n=nrow(et.mocuhyb3vmocchyb3$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb4vmocchyb4 <- topTags(et.mocuhyb4vmocchyb4, n=nrow(et.mocuhyb4vmocchyb4$table))$table %>% 
  rownames_to_column(var = "gene")


all_3v4 <- topTags(et.mocu.3v4, n=nrow(et.mocu.3v4$table))$table %>% 
  rownames_to_column(var = "gene") %>%
  inner_join(recip, by = c('gene' ='mocu_gene')) %>%
  inner_join(mocc_3v4, by = 'gene') %>% 
  inner_join(hyb_3v4, by = 'gene') %>%
  inner_join(mocuhyb3vmocchyb3, by = 'gene') %>%
  inner_join(mocuhyb4vmocchyb4, by = 'gene') %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  rename(logFC.mocu_nr_gt = logFC.x, logCPM.mocu_nr_gt = logCPM.x, PValue.mocu_nr_gt = PValue.x, FDR.mocu_nr_gt = FDR.x,
         logFC.mocc_nr_gt = logFC.y, logCPM.mocc_nr_gt = logCPM.y, PValue.mocc_nr_gt = PValue.y, FDR.mocc_nr_gt = FDR.y,
         logFC.hyb_nr_gt = logFC.x.x, logCPM.hyb_nr_gt = logCPM.x.x, PValue.hyb_nr_gt = PValue.x.x, FDR.hyb_nr_gt = FDR.x.x,
         logFC.hyb_gt = logFC.y.y, logCPM.hyb_gt = logCPM.y.y, PValue.hyb_gt = PValue.y.y, FDR.hyb_gt = FDR.y.y,
         logFC.hyb_nr = logFC, logCPM.hyb_nr = logCPM, PValue.hyb_nr = PValue, FDR.hyb_nr = FDR) %>% 
  mutate(mocc_gene=gsub("\\@KH.+", "", mocc_gene))
  
all_3v4 %>% 
  write_csv("3v4_Mol_GG_edgeR_results.csv")

mocc_3v6 <- topTags(et.mocc.3v6, n=nrow(et.mocc.3v6$table))$table %>% 
  rownames_to_column(var = "gene")
hyb_3v6 <- topTags(et.hyb.3v6, n=nrow(et.hyb.3v6$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb3vmocchyb3 <- topTags(et.mocuhyb3vmocchyb3, n=nrow(et.mocuhyb3vmocchyb3$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb6vmocchyb6 <- topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table %>% 
  rownames_to_column(var = "gene")

all_3v6 <- topTags(et.mocu.3v6, n=nrow(et.mocu.3v6$table))$table %>% 
  rownames_to_column(var = "gene") %>%
  inner_join(recip, by = c('gene' ='mocu_gene')) %>%
  inner_join(mocc_3v6, by = 'gene') %>% 
  inner_join(hyb_3v6, by = 'gene') %>%
  inner_join(mocuhyb3vmocchyb3, by = 'gene') %>%
  inner_join(mocuhyb6vmocchyb6, by = 'gene') %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  rename(logFC.mocu_tb_gt = logFC.x, logCPM.mocu_tb_gt = logCPM.x, PValue.mocu_tb_gt = PValue.x, FDR.mocu_tb_gt = FDR.x,
         logFC.mocc_tb_gt = logFC.y, logCPM.mocc_tb_gt = logCPM.y, PValue.mocc_tb_gt = PValue.y, FDR.mocc_tb_gt = FDR.y,
         logFC.hyb_tb_gt = logFC.x.x, logCPM.hyb_tb_gt = logCPM.x.x, PValue.hyb_tb_gt = PValue.x.x, FDR.hyb_tb_gt = FDR.x.x,
         logFC.hyb_gt = logFC.y.y, logCPM.hyb_gt = logCPM.y.y, PValue.hyb_gt = PValue.y.y, FDR.hyb_gt = FDR.y.y,
         logFC.hyb_tb = logFC, logCPM.hyb_tb = logCPM, PValue.hyb_tb = PValue, FDR.hyb_tb = FDR) %>% 
  mutate(mocc_gene=gsub("\\@KH.+", "", mocc_gene))
all_3v6 %>% 
  write_csv("3v6_Mol_GG_edgeR_results.csv")

library(readxl)
insitu_data <- read_excel("KHID-UniqueName-URLs-InSitu-COMPLETE.xlsx") #readxl library

all_all <- inner_join(all_3v4, all_4v6) %>%
  inner_join(all_3v6) %>% 
  right_join(insitu_data, ., by = c("KHID" = "KH.id")) %>% 
  select(-everything(),"gene","mocc_gene",everything()) #%>% 
  #mutate(mocc_gene=gsub("\\@KH.+", "", mocc_gene))
  


dim(all_all)
all_all %>% 
  write_csv("all_Mol_GG_edgeR_results_tblastx.csv")
#Loading the rvest package
library('rvest')
library(tibble)
library(data.table)
library(readxl)

#Scrap txt file from url and skip metadata



# Allele specific analysis  -----------------------------------------------


mol3vhyb3_genes <- rownames_to_column(topTags(et.molf3vshybf3, n=nrow(et.molf3vshybf3$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < 0.1) %>% 
  select(gene)
mol4vhyb4_up_genes <- rownames_to_column(topTags(et.molf4vshybf4, n=nrow(et.molf4vshybf4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.1) %>% 
  select(gene)
mol6vhyb6_genes <- rownames_to_column(topTags(et.molf6vshybf6, n=nrow(et.molf6vshybf6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < 0.1) %>% 
  select(gene)

mocuhyb3vmocchyb3_FDR_genes <- rownames_to_column(topTags(et.mocuhyb3vmocchyb3, n=nrow(et.mocuhyb3vmocchyb3$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < 0.1)
mocuhyb4vmocchyb4_up_genes <- rownames_to_column(topTags(et.mocuhyb4vmocchyb4, n=nrow(et.mocuhyb4vmocchyb4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.1) %>% 
  select(gene)
mocuhyb6vmocchyb6_up_genes <- rownames_to_column(topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.1) %>% 
  select(gene)

mocu3vmocc3_FDR_genes <- rownames_to_column(topTags(et.mocu3vmocc3, n=nrow(et.mocu3vmocc3$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < 0.1)
mocu4vmocc4_up_genes <- rownames_to_column(topTags(et.mocu4vmocc4, n=nrow(et.mocu4vmocc4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.1) %>% 
  select(gene)
mocu6vmocc6_genes <- rownames_to_column(topTags(et.mocu4vmocc4, n=nrow(et.mocu6vmocc6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5) %>% 
  select(gene)


mocu_4v6_up_genes <- rownames_to_column(topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.1) %>% 
  select(gene)
hyb_4v6_up_genes <- rownames_to_column(topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.1) %>% 
  select(gene)
mocc_4v6_down_genes <- rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < 0) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id))  %>%
  inner_join(mocc_pseudo_filtered, by = c("KH.id"="KH")) %>% filter(cate == "high") %>%  View()
  select(gene)
mocc_4v6_up_genes <- rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% 
  select(gene)

mocu_3v4_up_genes <- rownames_to_column(topTags(et.mocu.3v4, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% 
  select(gene)
hyb_3v4_up_genes <- rownames_to_column(topTags(et.hyb.3v4, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% 
  select(gene)
mocc_3v4_up_genes <- rownames_to_column(topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% 
  select(gene)
rownames_to_column(topTags(et.mocuhyb4vmocchyb4, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  filter((gene %in% mocu_4v6_up_genes$gene) & !(gene %in% hyb_4v6_up_genes$gene) & (gene %in% mocc_4v6_up_genes$gene)) %>%
  mutate(allele = logFC > 0) %>% count()

overlap_allele.05 <- rownames_to_column(topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < 0.05, gene %in% hyb_4v6_up_genes$gene & gene %in% mocu_4v6_up_genes$gene) %>% 
  mutate(allele = logFC > 0)
overlap_allele.all <- rownames_to_column(topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(gene %in% hyb_4v6_up_genes$gene & gene %in% mocu_4v6_up_genes$gene) %>% 
  mutate(allele = logFC > 0)
overlap_allele <-rbind(overlap_allele.05,overlap_allele.all)

hybUp_allele.05 <- rownames_to_column(topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < 0.05, gene %in% hyb_4v6_up_genes$gene) %>% 
  mutate(allele = logFC > 0)
hybUp_allele.all <- rownames_to_column(topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < gene %in% hyb_4v6_up_genes$gene) %>% 
  mutate(allele = logFC > 0)
hybUp_allele <-rbind(hybUp_allele.05,hybUp_allele.all)

hybAll_allele.05 <- rownames_to_column(topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(FDR < 0.05) %>% 
  mutate(allele = logFC > 0)
hybAll_allele.all <- rownames_to_column(topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table, "gene") %>%
  as_tibble() %>%
  mutate(allele = logFC > 0)
hybAll_allele <-rbind(hybAll_allele.05,hybAll_allele.all)  

overlap_allele.p <- ggplot(overlap_allele, aes(x=logFC, fill = allele, color = allele)) + 
  geom_histogram(data = overlap_allele.all, fill = "black", color = "gray", alpha = 0.2) +
  geom_histogram(data = overlap_allele.05) +
  xlim(-12,12) + 
  scale_fill_manual(values=c("FALSE"="lightpink1","TRUE"="slategray1"), labels = c("M. occulta", "M. oculata")) + 
  scale_color_manual(values=c("FALSE"="red","TRUE"="blue"), labels = c("M. occulta", "M. oculata")) +
  labs(title = "Allele specific logFC of overlapping hybrid & M. oculata up regulated genes at 6hpf") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = 2)

hybUp_allele.p <- ggplot(hybUp_allele, aes(x=logFC, fill = allele, color = allele)) + 
  geom_histogram(data = hybUp_allele.all, fill = "black", color = "gray", alpha = 0.2) +
  geom_histogram(data = hybUp_allele.05) +
  xlim(-12,12) + 
  scale_fill_manual(values=c("FALSE"="lightpink1","TRUE"="slategray1"), labels = c("M. occulta", "M. oculata")) + 
  scale_color_manual(values=c("FALSE"="red","TRUE"="blue"), labels = c("M. occulta", "M. oculata"))+
  labs(title = "Allele specific logFC of hybrid up regulated genes at 6hpf") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = 2)

hybAll_allele.p <- ggplot(hybAll_allele, aes(x=logFC, fill = allele, color = allele)) +
  geom_histogram(data = hybAll_allele.all, fill = "black", color = "gray", alpha = 0.2) +
  geom_histogram(data = hybAll_allele.05) +
  xlim(-12,12) + 
  scale_fill_manual(values=c("FALSE"="lightpink1","TRUE"="slategray1"), labels = c("M. occulta", "M. oculata")) + 
  scale_color_manual(values=c("FALSE"="red","TRUE"="blue"), labels = c("M. occulta", "M. oculata")) +
  labs(title = "Allele specific expression of all hybrid DEG at 6hpf") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = 2)
pdf("allele_hyb_tb_DE.pdf")
  hybAll_allele.p + theme_minimal()
dev.off()

pdf("allele_hyb_tb_all.pdf")
  grid.arrange(hybAll_allele.p, hybUp_allele.p, overlap_allele.p)
dev.off()
# Print count files -------------------------------------------------------

print_counts <- function(file_name, species_data, gene_sp){
  recip_file <- read_csv(file_name, col_names = c("ref",gene_sp))
  left_join(recip_file, species_data)%>% 
    write_excel_csv(paste0(file_name, "counts.csv"))
}

mocu_files <- list.files(pattern = "\\.mocu_GG.recip$")
for (i in mocu_files){
  print_counts(i, Mocu.data, "mocu_gene")
}

mocu.x.meta.recip <- mocu_files %>%
  map(read_csv, col_names = FALSE) %>%
  reduce(bind_rows)
mocu.x.meta.recip

mocc.x.meta.recip <- mocc_files %>%
  map(read_csv, col_names = FALSE) %>%
  reduce(bind_rows)
mocc.x.meta.recip


# ebf counts ---------

cpm_counts <- cpm(d3,normalized.lib.sizes=TRUE) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble()

mocu_ebf <- c("MOCU.TRINITY_GG_11616_c0_g2@KH.L24.10.v1.A.SL1-1")
mocu_islet <- c("MOCU.TRINITY_GG_4354_c0_g1@KH.L152.2.v1.A.nonSL5-1")
mocu_nova <- c("MOCU.TRINITY_GG_6526_c4_g1@KH.C11.417.v1.B.ND1-1")
mocu_ngn <- c("MOCU.TRINITY_GG_5983_c6_g1@KH.C6.129.v1.R.nonSL4-1")
mocu_ChaT <- c("MOCU.TRINITY_GG_1324_c2_g1@KH.C1.498.v1.A.SL1-1")
mocu_onecut <- c("MOCU.TRINITY_GG_4392_c0_g1@KH.C6.185.v1.A.SL1-1")
mocc_ebf <- c("MOCC.TRINITY_GG_11361_c1_g1@KH.L24.10.v2.A.ND2-2")
mocc_ChaT <- c("MOCC.TRINITY_GG_17058_c13_g1@KH.C1.498.v1.A.SL1-1")
mocc_ngn <- c("MOCC.TRINITY_GG_5757_c1_g1@KH.C6.129.v1.R.nonSL4-1")


hyb_ocu.ChaT<- hyb_ocu.data %>% 
  filter(mocu_gene %in% mocu_ChaT) %>% 
  summarise_if(is.numeric, sum) %>% 
  gather(key = time, value = counts)

plot_parents_gene_cpm <- function(gene, plot_title, cpm_count = cpm_counts){
cpm_count %>% 
  dplyr::filter(rowname == gene) %>% 
  select(Mocu_f3,Mocu_f4,Mocu_f6,Mocc_f3,Mocc_f4,Mocc_f6) %>% 
  gather(key = time, value = counts)%>% 
  separate(time, sep = "_f(?=[:digit:])", into = c("speices", "hour")) %>% 
  ggplot(aes(x=hour, y = counts, fill = speices)) + 
  geom_bar(stat = 'identity', position=position_dodge(), color = rep(c("blue", "red"),3))  + labs(y = "CPM*", x = "hours post fertilization", caption = "* normalized by library size", title = paste(plot_title,"Gene expression")) +
  scale_fill_manual(values = c("lightpink1","slategray1"), labels = c("M. occulta", "M. oculata"), name = "Species") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
}

plot_allele_gene_cpm <- function(gene, plot_title, cpm_count = cpm_counts){
cpm_count %>% 
  dplyr::filter(rowname == gene) %>% 
  select(hyb_ocu_f3,hyb_ocu_f4,hyb_ocu_f6,hyb_occ_f3,hyb_occ_f4,hyb_occ_f6) %>% 
  gather(key = time, value = counts)%>% 
  separate(time, sep = "_f(?=[:digit:])", into = c("speices", "hour")) %>% 
  ggplot(aes(x=hour, y = counts, fill = speices)) + 
  geom_bar(stat = 'identity', position=position_dodge(), color = rep(c("blue", "red"),3))  + labs(y = "CPM*", x = "hours post fertilization", caption = "* normalized by library size", title = paste(plot_title,"Allelic expression")) +
  scale_fill_manual(values = c("lightpink1","slategray1"), labels = c("M. occulta", "M. oculata"), name = "Parent allele") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
}

mg_genes <- c(mocu_ebf,mocu_islet,mocu_nova,mocu_ngn,mocu_ChaT,mocu_onecut)

plot_parents_gene_cpm("MOCU.TRINITY_GG_2107_c14_g1|m.4751", "unknown")
svg('ebf_cpm.svg')
plot_parents_gene_cpm(mocu_ebf, "Ebf")
dev.off()
svg('ebf_allele_cpm.svg')
plot_allele_gene_cpm(mocu_ebf, "Ebf")
dev.off()
svg('islet_cpm.svg')
plot_parents_gene_cpm(mocu_islet, "islet")
dev.off()
svg('islet_allele_cpm.svg')
plot_allele_gene_cpm(mocu_islet, "islet")
dev.off()
svg('nova_cpm.svg')
plot_parents_gene_cpm(mocu_nova, "nova")
dev.off()
svg('nova_allele_cpm.svg')
plot_allele_gene_cpm(mocu_nova, "nova")
dev.off()
svg('ngn_cpm.svg')
plot_parents_gene_cpm(mocu_ngn, "ngn")
dev.off()
svg('ngn_allele_cpm.svg')
plot_allele_gene_cpm(mocu_ngn, "ngn")
dev.off()
svg('ChaT_cpm.svg')
plot_parents_gene_cpm(mocu_ChaT, "ChaT")
dev.off()
svg('ChaT_allele_cpm.svg')
plot_allele_gene_cpm(mocu_ChaT, "ChaT")
dev.off()
svg('onecut_cpm.svg')
plot_parents_gene_cpm(mocu_onecut, "Onecut")
dev.off()
svg('onecut_allele_cpm.svg')
plot_allele_gene_cpm(mocu_onecut, "Onecut")
dev.off()

ggplot(cpm_counts, aes(x=log2(hyb_f6)-log2(Mocu_f6), y=log2(hyb_f6)-log2(Mocc_f6))) + 
  geom_point(alpha = 0.5)

ggplot(cpm_counts, aes(x=log2(Mocu_f6/Mocc_f6), y=log2(hyb_ocu_f6/hyb_occ_f6))) + 
  geom_point(alpha = 0.5)

write_tsv(cpm_counts, "Mol_cpm.txt")
et.mocuhyb3vmocchyb3["MOCU.TRINITY_GG_1324_c2_g1@KH.C1.498.v1.A.SL1-1",]
et.mocuhyb4vmocchyb4["MOCU.TRINITY_GG_1324_c2_g1@KH.C1.498.v1.A.SL1-1",]
et.mocuhyb6vmocchyb6["MOCU.TRINITY_GG_1324_c2_g1@KH.C1.498.v1.A.SL1-1",]
et.mocc.4v6["MOCU.TRINITY_GG_1324_c2_g1@KH.C1.498.v1.A.SL1-1",]
cpm(d3,normalized.lib.sizes=TRUE) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() #['MOCU.TRINITY_GG_1324_c2_g1@KH.C1.498.v1.A.SL1-1',]

mocc_non_immune <- read_tsv("../mocc.x.FASTA_Metamorphosis_-_not_Immunity.output", col_names = c("meta","mocc")) %>% 
  unique() %>% 
  mutate(type = "non_immune")
mocc_innate_imm <- read_tsv("../mocc.x.FASTA_Innate_Immunity_Metamorphosis.output", col_names = c("meta","mocc")) %>% 
  unique() %>% 
  mutate(type = "innate_immune")
mocc_housekeeping <- read_tsv("../mocc.x.FASTA_Housekeeping_Proteins.output", col_names = c("meta","mocc")) %>% 
  unique() %>% 
  mutate(type = "housekeeping")
mocc_other_meta <- read_tsv("../mocc.x.FATSA_other.output", col_names = c("meta","mocc")) %>% 
  unique() %>% 
  mutate(type = "other")
mocc_mata <- bind_rows(mocc_non_immune,mocc_innate_imm,mocc_housekeeping,mocc_other_meta)

mocc_meta <- inner_join(mocc_mata, Mocc.data, by = c("mocc" ="mocc_gene")) %>% 
  group_by(type) %>% 
  summarise_at(vars(Mocc_f3:Mocc_f6), mean, na.rm = TRUE)

mocu_non_immune <- read_tsv("../mocu.x.FASTA_Metamorphosis_-_not_Immunity.output", col_names = c("meta","mocu")) %>% 
  unique() %>% 
  mutate(type = "non_immune")
mocu_innate_imm <- read_tsv("../mocu.x.FASTA_Innate_Immunity_Metamorphosis.output", col_names = c("meta","mocu")) %>% 
  unique() %>% 
  mutate(type = "innate_immune")
mocu_housekeeping <- read_tsv("../mocu.x.FASTA_Housekeeping_Proteins.output", col_names = c("meta","mocu")) %>% 
  unique() %>% 
  mutate(type = "housekeeping")
mocu_other_meta <- read_tsv("../mocu.x.FATSA_other.output", col_names = c("meta","mocu")) %>% 
  unique() %>% 
  mutate(type = "other")

cpm_counts <- as.data.frame(cpm(d3)) %>% 
  rownames_to_column(var = "mocu_gene")

mocu_mata <- bind_rows(mocu_non_immune,mocu_innate_imm,mocu_housekeeping,mocu_other_meta)
inner_join(mocu_mata, Mocu.data, by = c("mocu" ="mocu_gene")) %>%
  group_by(type) %>% 
  summarise_at(vars(Mocu_f3:Mocu_f6), mean, na.rm = TRUE)
mocu_meta <- inner_join(mocu_mata, cpm_counts, by = c("mocu" ="mocu_gene"))
mocu_meta %>% 
  gather(time, counts,Mocu_f3:Mocu_f6) %>% 
  ggplot(aes(x = type, y = counts, fill = time)) + geom_boxplot() + ylim(c(0.1, 1000)) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)

manx <- c("MOCU.TRINITY_GG_4725_c0_g1@KH.C2.957.v1.A.nonSL1-1")



# Transcription Factor genes ----------------------------------------------

tf_url <- 'http://ghost.zool.kyoto-u.ac.jp/TF_KH.html'
transcription_factors <- tf_url %>%
  read_html() %>%
  html_nodes('.Table1') %>%
  html_table(fill = TRUE)
transcription_factors <- transcription_factors[[1]][1:3] %>% 
  as.tibble()

bHLH <- slice(transcription_factors, (grep('bHLH', transcription_factors$X3, value = FALSE)[1]+3):grep('bZIP', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'bHLH')
bZIP <- slice(transcription_factors, (grep('bZIP', transcription_factors$X3, value = FALSE)[1]+3):grep('Ets', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'bZip')
Ets <- slice(transcription_factors, (grep('Ets', transcription_factors$X3, value = FALSE)[1]+3):grep('Fox', transcription_factors$X3, value = FALSE)[1]-1)%>% mutate(type = 'Ets')
Fox <- slice(transcription_factors, (grep('Fox', transcription_factors$X3, value = FALSE)[1]+3):grep('HMG', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'Fox')
HMG <- slice(transcription_factors, (grep('HMG', transcription_factors$X3, value = FALSE)[1]+3):grep('homeobox', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'HMG')
homeobox <- slice(transcription_factors, (grep('homeobox', transcription_factors$X3, value = FALSE)[1]+3):grep('Nuclear Receptor', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'homeobox')
Nuclear_receptor <- slice(transcription_factors, (grep('Nuclear Receptor', transcription_factors$X3, value = FALSE)[1]+3):grep('T-box', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'Nuclear Receptor')
T_box <- slice(transcription_factors, (grep('T-box', transcription_factors$X3, value = FALSE)[1]+3):grep('Other TFs', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'T-box')
Other <- slice(transcription_factors, (grep('Other TFs', transcription_factors$X3, value = FALSE)[1]+3):grep('Genes that may or may not encode a seqence-specific transcription factor', transcription_factors$X3, value = FALSE)[1]-1) %>% mutate(type = 'other')
possibles <- slice(transcription_factors, (grep('Genes that may or may not encode a seqence-specific transcription factor', transcription_factors$X3, value = FALSE)[1]+2):n()) %>% mutate(type = 'possibles')

TF <- bHLH %>% 
  bind_rows(bZIP) %>% 
  bind_rows(Ets) %>% 
  bind_rows(Fox) %>% 
  bind_rows(HMG) %>% 
  bind_rows(homeobox) %>% 
  bind_rows(Nuclear_receptor) %>% 
  bind_rows(T_box) %>% 
  bind_rows(Other) %>% 
  bind_rows(possibles) %>% 
  mutate(cat = 'TF')

# Signaling Molecule genes ------------------------------------------------

signaling_url <- 'http://ghost.zool.kyoto-u.ac.jp/ST_KH.html'
signaling_molecules <- signaling_url %>%
  read_html() %>%
  html_nodes('.Table1') %>% 
  html_table(fill = TRUE)
signaling_molecules <- signaling_molecules[[1]][1:3] %>% 
  as.tibble()

headers <- c('RTK', 'TGF', 'Wnt','hedgehog', 'JAK', 'NFkB', 'Notch')

RTK <- slice(signaling_molecules, (grep('RTK', signaling_molecules$X3, value = FALSE)[1]+1):grep('TGF', signaling_molecules$X3, value = FALSE)[1]-1) %>% mutate(type = 'RTK')
TGF <- slice(signaling_molecules, (grep('TGF', signaling_molecules$X3, value = FALSE)[1]+1):grep('Wnt', signaling_molecules$X3, value = FALSE)[1]-1) %>% mutate(type = 'TGF')
Wnt <- slice(signaling_molecules, (grep('Wnt', signaling_molecules$X3, value = FALSE)[1]+3):grep('hedgehog', signaling_molecules$X3, value = FALSE)[1]-1)%>% mutate(type = 'Wnt')
hedgehog <- slice(signaling_molecules, (grep('hedgehog', signaling_molecules$X3, value = FALSE)[1]+3):grep('JAK', signaling_molecules$X3, value = FALSE)[1]-1) %>% mutate(type = 'hedgehog')
JAK <- slice(signaling_molecules, (grep('JAK', signaling_molecules$X3, value = FALSE)[1]+3):grep('NFkB', signaling_molecules$X3, value = FALSE)[1]-1) %>% mutate(type = 'JAK')
NFkB <- slice(signaling_molecules, (grep('NFkB', signaling_molecules$X3, value = FALSE)[1]+3):grep('Notch', signaling_molecules$X3, value = FALSE)[1]-1) %>% mutate(type = 'NFkB')
Notch <- slice(signaling_molecules, (grep('Notch', signaling_molecules$X3, value = FALSE)[1]+2):n()) %>% mutate(type = 'Notch')

signaling <- RTK %>% 
  bind_rows(TGF) %>% 
  bind_rows(Wnt) %>% 
  bind_rows(hedgehog) %>% 
  bind_rows(JAK) %>% 
  bind_rows(NFkB) %>% 
  bind_rows(Notch) %>% 
  mutate(cat = 'signaling')

# Zinc finger genes -------------------------------------------------------

zF_url <- 'http://ghost.zool.kyoto-u.ac.jp/ZF_KH.html'
zinc_finger <- zF_url %>%
  read_html() %>%
  html_nodes('.Table1') %>% 
  html_table(fill = TRUE)
zinc_finger <- zinc_finger[[1]][1:3] %>% 
  as.tibble()

zf_TF <- slice(zinc_finger, (grep('Genes that encode a possible transcription factor with Zinc Finger motifs', zinc_finger$X3, value = FALSE)[1]+3):grep('Genes for Zinc Finger proteins that are NOT transcription factors', zinc_finger$X3, value = FALSE)[1]-1) %>% mutate(type = 'zf_TF')
zf_nonTF <- slice(zinc_finger, (grep('Genes for Zinc Finger proteins that are NOT transcription factors', zinc_finger$X3, value = FALSE)[1]+2):n()) %>% mutate(type = 'zf_nonTF')


zf <- zf_TF %>% 
  bind_rows(zf_nonTF) %>% 
  mutate(cat = 'zf')


# All categories  ---------------------------------------------------------

annotated_genes <- TF %>% 
  bind_rows(signaling) %>% 
  bind_rows(zf) %>% 
  rename(KH.id = X3, ghost.id = X2, gene_name = X1)

# Volcano plots -----------------------------------------------------------

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

install.packages("devtools")
devtools::install_github("kevinblighe/EnhancedVolcano")
library("EnhancedVolcano")

EnhancedVolcano(et.mocu.4v6,
                lab = rownames(et.mocu.4v6),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8))


# UpSetR plots ------------------------------------------------------------
#install.packages("UpSetR")
library(UpSetR)
listInput <- list('M. oculata neur v gast up' = mocu.3v4.up$gene, 'M. oculata neur v gast down' = mocu.3v4.down$gene, 'M. occulta neur v gast up' = mocc.3v4.up$gene, 'M. occulta neur v gast down' = mocc.3v4.down$gene, 'M. hybrid neur v gast up' = hyb.3v4.up$gene, 'M. hybrid neur v gast down' = hyb.3v4.down$gene,
                  'M. oculata tail v neur up' = mocu.4v6.up$gene, 'M. oculata tail v neur down' = mocu.4v6.down$gene, 'M. occulta tail v neur up' = mocc.4v6.up$gene, 'M. occulta tail v neur down' = mocc.4v6.down$gene, 'M. hybrid tail v neur up' = hyb.4v6.up$gene, 'M. hybrid tail v neur down' = hyb.4v6.down$gene,
                  'M. oculata tail v gast up' = mocu.3v6.up$gene, 'M. oculata tail v gast down' = mocu.3v6.down$gene, 'M. occulta tail v gast up' = mocc.3v6.up$gene, 'M. occulta tail v gast down' = mocc.3v6.down$gene, 'M. hybrid tail v gast up' = hyb.3v6.up$gene, 'M. hybrid tail v gast down' = hyb.3v6.down$gene)
pdf('molgula_de_upsetplot.pdf')
upset(fromList(listInput), order.by = "freq", nsets = 18)
dev.off()
