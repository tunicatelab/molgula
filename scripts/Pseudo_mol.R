library(readxl)
library(tidyverse)

data_to_annotate <- read_excel("~/Desktop/Alberto/GPCR-list.xls", skip = 5) %>% 
  mutate(`Accession No`=gsub('>', '',`Accession No`))
data_to_annotate

lookup_table <- read_tsv("~/Desktop/Alberto/KH-JGI.blast", col_names = FALSE) %>% 
  select(X1, X2) %>% 
  write_csv("~/Desktop/Alberto/KH.xJGI.lookup.csv")
lookup_table

data_to_annotate %>% 
  left_join(lookup_table, by = c(`Accession No`="X1")) %>% 
  write_csv("~/Desktop/Alberto/GPCR-list.csv")

mocu_pseudo <- read_delim("~/Desktop/tails/scaffolded_gene_models/mocu_pseudo.txt", delim = ' ',col_names = c("KH", "subject", "match", "eval", "stop_high","stop_low","FS_high", "FS_low"))
mocc_pseudo <- read_delim("~/Desktop/tails/scaffolded_gene_models/mocc_pseudo.txt", delim = ' ',col_names = c("KH", "subject", "match", "eval", "stop_high","stop_low","FS_high", "FS_low"))

mocu_stop_high <- mocu_pseudo %>% filter(stop_high > 0) %>% select(KH) %>% unique()
mocu_no_stop_high <- mocu_pseudo %>% filter(stop_high == 0) %>% select(KH) %>% unique()

mocc_stop_high <- mocc_pseudo %>% filter(stop_high > 0) %>% select(KH) %>% unique()
mocc_no_stop_high <- mocc_pseudo %>% filter(stop_high == 0) %>% select(KH) %>% unique()

sum(mocu_stop_high$KH %in% mocc_stop_high$KH)/length(mocu_stop_high$KH)

mocc_pseudo <- mocc_pseudo %>%
  mutate(cate = case_when(stop_high > 0 | FS_high > 0 ~ 'high',
                          stop_low > 0 | FS_low > 0 ~ 'low',
                          TRUE ~ 'coding')) %>% 
  group_by(cate) %>% 
  select(-subject, -match, -eval) %>% 
  mutate(KH=gsub("\\.v.+", "", KH)) %>%
  unique() 

mocc_coding <- mocc_pseudo %>%
  filter(cate == 'coding') %>% 
  select(KH,cate)#summarise(n())

mocc_pseudo_filtered <- mocc_pseudo %>% 
  filter(!KH %in% mocc_coding$KH)

mocc_pseudo_filtered %>%
  group_by(cate) %>%
  summarise(n())


mocu_pseudo <- mocu_pseudo %>%
  mutate(cate = case_when(stop_high > 0 | FS_high > 0 ~ 'high',
                          stop_low > 0 | FS_low > 0 ~ 'low',
                          TRUE ~ 'coding')) %>%
  select(-subject, -match, -eval) %>% 
  mutate(KH=gsub("\\.v.+", "", KH)) %>%
  unique() 
  
  mocu_coding <- mocu_pseudo %>%
    filter(cate == 'coding') %>% 
    select(KH,cate)#summarise(n())
  
mocu_pseudo_filtered <-   mocu_pseudo %>% 
    filter(!KH %in% mocu_coding$KH) 

mocu_pseudo_filtered %>%
    group_by(cate) %>%
    summarise(n())

mocc_pseudo_filtered %>% 
  filter(!KH %in% mocu_pseudo_filtered$KH) %>% 
  inner_join(mocc_4v6, by =c("KH" = "KH.id")) %>% View()
length(mocc_pseudo_filtered)
  
group_by(cate) %>%
  summarise(n())
         