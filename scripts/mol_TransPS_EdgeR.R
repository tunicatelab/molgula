
# Libraries ---------------------------------------------------------------

library(edgeR)
library(tidyverse)
library(gridExtra)

# Load_count_files --------------------------------------------------------
setwd("/Users/elijahlowe/Desktop/tails/scaffolded_gene_models/")


Mocu.3hpf <- read_tsv("SRR1197522_mocu_GG.xprs") %>% 
  dplyr::rename(Mocu_f3 = eff_counts, mocu_gene = target_id)
Mocu.3hpf #test to see if loaded properly
Mocu.4hpf <- read_tsv("SRR1197965_mocu_GG.xprs") %>% 
  dplyr::rename(Mocu_f4 = eff_counts, mocu_gene = target_id)
Mocu.6hpf <- read_tsv("SRR1197972_mocu_GG.xprs") %>% 
  dplyr::rename(Mocu_f6 = eff_counts, mocu_gene = target_id)

nrow(Mocu.3hpf) == nrow(Mocu.4hpf)

#Now let's combine them all into one dataset

Mocu.data <- Mocu.3hpf %>% 
  inner_join(Mocu.4hpf, by = "mocu_gene") %>% 
  inner_join(Mocu.6hpf, by = "mocu_gene") %>% 
  select(mocu_gene, Mocu_f3, Mocu_f4, Mocu_f6)


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
  select(mocc_gene, Mocc_f3, Mocc_f3.2, Mocc_f4, Mocc_f5, Mocc_f6)


# Load hybrid samples -----------------------------------------------------

hyb_ocu.3hpf <- read_tsv("SRR1198321_mocu_GG.xprs") %>% 
  dplyr::rename(hyb_ocu_f3 = eff_counts, mocu_gene = target_id)
hyb_ocu.4hpf <- read_tsv("SRR1198337_mocu_GG.xprs") %>% 
  dplyr::rename(hyb_ocu_f4 = eff_counts, mocu_gene = target_id)
hyb_ocu.6hpf <- read_tsv("SRR1198346_mocu_GG.xprs") %>% 
  dplyr::rename(hyb_ocu_f6 = eff_counts, mocu_gene = target_id)

hyb_ocu.data <- hyb_ocu.3hpf %>% 
  inner_join(hyb_ocu.4hpf, by = "mocu_gene") %>% 
  inner_join(hyb_ocu.6hpf, by = "mocu_gene") %>% 
  select(mocu_gene, hyb_ocu_f3, hyb_ocu_f4, hyb_ocu_f6)

hyb_occ.3hpf <- read_tsv("SRR1198321_mocc_GG.xprs") %>% 
  dplyr::rename(hyb_occ_f3 = eff_counts, mocc_gene = target_id)
hyb_occ.4hpf <- read_tsv("SRR1198337_mocc_GG.xprs") %>% 
  dplyr::rename(hyb_occ_f4 = eff_counts, mocc_gene = target_id)
hyb_occ.6hpf <- read_tsv("SRR1198346_mocc_GG.xprs") %>% 
  dplyr::rename(hyb_occ_f6 = eff_counts, mocc_gene = target_id)

hyb_occ.data <- hyb_occ.3hpf %>% 
  inner_join(hyb_occ.4hpf, by = "mocc_gene") %>% 
  inner_join(hyb_occ.6hpf, by = "mocc_gene") %>% 
  select(mocc_gene, hyb_occ_f3, hyb_occ_f4, hyb_occ_f6)

# Join all datasets, and filter with reciprocal  --------------------------

recip <- read_csv('all-mocu_GG_TransPS.x.mocc_GG_TransPS.recip', col_types = list(col_character(),col_character()), col_names = c("mocu_gene","mocc_gene"))
all.data <- inner_join(Mocu.data, recip) %>% 
  inner_join(Mocc.data) %>% 
  inner_join(hyb_ocu.data) %>% 
  inner_join(hyb_occ.data) %>% 
  mutate(hyb_f3 = hyb_ocu_f3 + hyb_occ_f3, hyb_f4 = hyb_ocu_f4 + hyb_occ_f4, hyb_f6 = hyb_ocu_f6 + hyb_occ_f6) %>% 
  column_to_rownames('mocu_gene') %>% 
  select(-mocc_gene)

# all.data.h <- inner_join(Mocu.data, recip) %>% 
#   inner_join(Mocc.data) %>% 
#   inner_join(hyb_ocu.data) %>% 
#   inner_join(hyb_occ.data) %>% 
#   mutate(hyb_f3 = hyb_ocu_f3 + hyb_occ_f3, hyb_f4 = hyb_ocu_f4 + hyb_occ_f4, hyb_f6 = hyb_ocu_f6 + hyb_occ_f6) %>% 
#   select(-mocc_gene)

#We now have a data frame with all the data in it:
#View(all.data)


# Experimental setup ------------------------------------------------------

Group <- factor(paste(species=c("ocu","ocu","ocu","occ","occ","occ","occ","occ","hyb_ocu","hyb_ocu","hyb_ocu","hyb_occ","hyb_occ","hyb_occ","hyb","hyb","hyb"),
                      stage=c("f3","f4","f6","f3","f3","f4","f5","f6","f3","f4","f6","f3","f4","f6","f3","f4","f6"),sep="."))
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
#  keep <- rowSums(cpm(edger.data)>1) >=2
#  keep 
#  edger.data<-edger.data[keep,]

design <- model.matrix(~Group)
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

#housekeeping <- read_csv("mocu_housekeeping.list", col_names = "housekeeping")
#housekeeping1 <- read_csv("mocu_gimme.x.KH2012_housekeeping.csv", col_names = c("housekeeping","KH_ID"))
#d1 <- edger.data
#d1<-estimateCommonDisp(edger.data)
#d4 <- estimateGLMTagwiseDisp(d3,design)
d2 <- estimateGLMCommonDisp(edger.data,design,verbose=TRUE)
d3 <- estimateGLMTagwiseDisp(d2,design)

plotMDS(d3, main="edgeR MDS Plot")

et.mocu.3v4<-exactTest(d3,pair=ocuf3vsocuf4)
summary(de.mocu.3v4 <- decideTestsDGE(et.mocu.3v4, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu.3v4)]
jpeg('et.mocu.3v4.jpg')
plotSmear(et.mocu.3v4, ylim=c(-10,10), de.tags=detags, main="M. oculata neurula vs gastrula")
dev.off()

et.mocc.3v4<-exactTest(d3,pair=occf3vsoccf4)
summary(de.mocc.3v4 <- decideTestsDGE(et.mocc.3v4, p=0.05))
detags <- rownames(d3)[as.logical(de.mocc.3v4)]
jpeg('et.mocc.3v4.jpg')
plotSmear(et.mocc.3v4, ylim=c(-10,10), de.tags=detags, main="M. occulta neurula vs gastrula")
dev.off()

et.mocu.4v6<-exactTest(d3,pair=ocuf4vsocuf6)
summary(de.mocu.4v6 <- decideTestsDGE(et.mocu.4v6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocu.4v6)]
jpeg('et.mocu.4v6.jpg')
plotSmear(et.mocu.4v6, ylim=c(-10,10), de.tags=detags, main="M. oculata tailbud vs neurula")
dev.off()

et.mocc.4v6<-exactTest(d3,pair=occf4vsoccf6)
summary(de.mocc.4v6 <- decideTestsDGE(et.mocc.4v6, p=0.05))
detags <- rownames(d3)[as.logical(de.mocc.4v6)]
jpeg('et.mocc.4v6.jpg')
plotSmear(et.mocc.4v6, ylim=c(-10,10), de.tags=detags, main="M. occulta tailbud vs neurula")
dev.off()

et.hyb.3v4<-exactTest(d3,pair=hybf3vshybf4)
summary(de.hyb.3v4 <- decideTestsDGE(et.hyb.3v4, p=0.05))
detags <- rownames(d3)[as.logical(de.hyb.3v4)]
jpeg('et.hyb.3v4.jpg')
plotSmear(et.hyb.3v4, ylim=c(-10,10), de.tags=detags, main="Hybrid neurula vs gastrula")
dev.off()

et.hyb.4v6<-exactTest(d3,pair=hybf4vshybf6)
summary(de.hyb.4v6 <- decideTestsDGE(et.hyb.4v6, p=0.05))
detags <- rownames(d3)[as.logical(de.hyb.4v6)]
jpeg('et.hyb.4v6.jpg')
plotSmear(et.hyb.4v6, ylim=c(-10,10), de.tags=detags, main="Hybrid tailbud vs neurula")
dev.off()

et.mocuhyb3vmocchyb3<-exactTest(d3,pair=hocuf3vshoccf3)
summary(de.mocuhyb3vmocchyb3 <- decideTestsDGE(et.mocuhyb3vmocchyb3, p=0.05))
detags <- rownames(d3)[as.logical(de.mocuhyb3vmocchyb3)]
plotSmear(et.mocuhyb3vmocchyb3, ylim=c(-10,10), de.tags=detags, main="Hybrid 3hpf")

et.mocuhyb4vmocchyb4<-exactTest(d3,pair=hocuf4vshoccf4)
summary(de.mocuhyb4vmocchyb4 <- decideTestsDGE(et.mocuhyb4vmocchyb4, p=0.05))
detags <- rownames(d3)[as.logical(de.mocuhyb4vmocchyb4)]
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

rownames_to_column(topTags(et.mocu.3v4, n=nrow(et.mocu.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.3v4.up.csv")
rownames_to_column(topTags(et.mocu.3v4, n=nrow(et.mocu.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.3v4.down.csv")
rownames_to_column(topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.3v4.up.csv")
rownames_to_column(topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.3v4.down.csv")
rownames_to_column(topTags(et.hyb.3v4, n=nrow(et.hyb.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.3v4.up.csv")
rownames_to_column(topTags(et.hyb.3v4, n=nrow(et.hyb.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.3v4.down.csv")

rownames_to_column(topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.4v6.up.csv")
rownames_to_column(topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocu.4v6.down.csv")
rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.mocc.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.4v6.up.csv")
rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.mocc.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  inner_join(recip, by = c("gene" = "mocu_gene")) %>%
  mutate(gene = mocc_gene) %>% 
  select(-mocc_gene) %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.mocc.4v6.down.csv")
rownames_to_column(topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.4v6.up.csv")
rownames_to_column(topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < -1.5, FDR < 0.05) %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("et.hyb.4v6.down.csv")



# Allele specific analysis  -----------------------------------------------


mocu_4v6_up_genes <- rownames_to_column(topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% 
  select(gene)
hyb_4v6_up_genes <- rownames_to_column(topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC > 1.5, FDR < 0.05) %>% 
  select(gene)
mocc_4v6_down_genes <- rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  filter(logFC < 0) %>% 
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
hybAll_allele.p + theme_minimal()

grid.arrange(hybAll_allele.p, hybUp_allele.p, overlap_allele.p)
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

#install.packages("readxl")
library(readxl)

metamorph_blast <- read_excel("blast_results_table_metamorphosis_genes.xlsx", skip = 1)

metamorph_blast %>% 
  group_by(`Species studied`) %>% 
  count()

left_join(metamorph_blast, mocc.x.meta.recip, by = c("Genbank Accession" = "X1")) %>% 
  select("Species studied", "X2") %>% 
  View()
group_by(`Species studied`) %>% 
  summarise()

# ebf counts ---------
mocc_ebf <- c("MOCC.TRINITY_GG_11361_c1_g1@KH.L24.10.v2.A.ND2-2")
mocu_ebf <- c("MOCU.TRINITY_GG_11616_c0_g2@KH.L24.10.v1.A.SL1-1")
hyb_ocu.ebf <- hyb_ocu.data %>% 
  filter(mocu_gene %in% mocu_ebf) %>% 
  summarise_if(is.numeric, sum) %>% 
  gather(key = time, value = counts)
hyb_occ.ebf <- hyb_occ.data %>% 
  filter(mocc_gene %in% mocc_ebf) %>% 
  summarise_if(is.numeric, sum) %>% 
  gather(key = time, value = counts)
bind_rows(hyb_ocu.ebf,hyb_occ.ebf) %>% 
  separate(time, sep = "_f(?=[:digit:])", into = c("speices", "hour")) %>% 
  ggplot(aes(x=hour, y = counts, fill = speices)) + 
  geom_bar(stat = 'identity', position=position_dodge(), color = rep(c("blue", "red"),3)) + coord_flip() + labs(y = "effective counts", x = "hours post fertilization") +
  scale_fill_manual(values = c("lightpink1","slategray1"), labels = c("M. occulta", "M. oculata"), name = "Parent allele") +
  theme_minimal()

et.mocuhyb3vmocchyb3["MOCU.TRINITY_GG_11616_c0_g2@KH.L24.10.v1.A.SL1-1",]
et.mocuhyb4vmocchyb4["MOCU.TRINITY_GG_11616_c0_g2@KH.L24.10.v1.A.SL1-1",]
et.mocuhyb6vmocchyb6["MOCU.TRINITY_GG_11616_c0_g2@KH.L24.10.v1.A.SL1-1",]

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


# write all results -------------------------------------------------------


mocc_4v6 <- topTags(et.mocc.4v6, n=nrow(et.mocc.4v6$table))$table %>% 
  rownames_to_column(var = "gene")
hyb_4v6 <- topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb4vmocchyb4 <- topTags(et.mocuhyb4vmocchyb4, n=nrow(et.mocuhyb4vmocchyb4$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb6vmocchyb6 <- topTags(et.mocuhyb6vmocchyb6, n=nrow(et.mocuhyb6vmocchyb6$table))$table %>% 
  rownames_to_column(var = "gene")


topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table %>% 
  rownames_to_column(var = "gene") %>% 
  inner_join(mocc_4v6, by = 'gene') %>% 
  inner_join(hyb_4v6, by = 'gene') %>%
  inner_join(mocuhyb4vmocchyb4, by = 'gene') %>%
  inner_join(mocuhyb6vmocchyb6, by = 'gene') %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("4v6_Mol_GG_edgeR_results.csv")

mocc_3v4 <- topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table %>% 
  rownames_to_column(var = "gene")
hyb_3v4 <- topTags(et.hyb.3v4, n=nrow(et.hyb.3v4$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb3vmocchyb3 <- topTags(et.mocuhyb3vmocchyb3, n=nrow(et.mocuhyb3vmocchyb3$table))$table %>% 
  rownames_to_column(var = "gene")
mocuhyb4vmocchyb4 <- topTags(et.mocuhyb4vmocchyb4, n=nrow(et.mocuhyb4vmocchyb4$table))$table %>% 
  rownames_to_column(var = "gene")


topTags(et.mocu.3v4, n=nrow(et.mocu.3v4$table))$table %>% 
  rownames_to_column(var = "gene") %>% 
  inner_join(mocc_3v4, by = 'gene') %>% 
  inner_join(hyb_3v4, by = 'gene') %>%
  inner_join(mocuhyb3vmocchyb3, by = 'gene') %>%
  inner_join(mocuhyb4vmocchyb4, by = 'gene') %>% 
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>%
  write_csv("3v4_Mol_GG_edgeR_results.csv")

#Loading the rvest package
library('rvest')
library(tibble)
library(data.table)
library(readxl)

#Scrap txt file from url and skip metadata


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
