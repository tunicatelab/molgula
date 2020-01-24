#Loading the rvest package
library('rvest')
library(tibble)
library(tidyverse)
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
  as_tibble()

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
  as_tibble()

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
  as_tibble()

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

View(annotated_genes)


# get expression patterns -------------------------------------------------

expression_url <- ("http://ghost.zool.kyoto-u.ac.jp/cgi-bin/txtgetkh.cgi?inkey=citb008o08&source=kh2013")
read_html(expression_url)





# Label microarray probes -------------------------------------------------

data_to_annotate <- read_excel("~/Desktop/Alberto/AlbertoRMA1-9b-LOSTPROBESETS-plusSequences.xlsx")
data_to_annotate

read_jgi_to_KH <- fread("~/Desktop/Alberto/KH-JGI.blast", select = c(1,2), col.names = c("JGI", "KH"), header = FALSE)
jgi_to_KH <- read_jgi_to_KH %>% 
  mutate(KH=gsub("\\.v.+", "", KH)) %>% 
  mutate(KH=gsub("KH2012:", "", KH)) %>% 
  unique() %>% 
  as_tibble()

all_probes_to_KH <- right_join(jgi_to_KH, probe_data, by = c("JGI" = "Composite Element Database Entry[aniseed]")) %>% 
  select('KH', 'Composite Element Name', 'JGI')

data_to_annotate %>% 
  left_join(all_probes_to_KH, by = c('Probeset' = 'Composite Element Name'), fill = "no annotation") %>%
  mutate('KH ID' = replace_na(KH,"No annotation"), other = ifelse(JGI %in% "", other, JGI)) %>% 
  select(-KH, -JGI) %>% 
  write_excel_csv('~/Desktop/Alberto/AlbertoRMA1-9b-LOSTPROBESETS-plusSequences+KH.csv')
