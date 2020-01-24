install.packages("httr")
install.packages("jsonlite")

library(httr)
library(jsonlite)

library(RCurl)
library("tidyverse")
library(dplyr)





path <- "https://www.aniseed.cnrs.fr/api/all_genes:?"
request <- GET(url = path, 
               query = list(
                 organism_id = 464)
)

request$status_code
response <- content(request, as = "text", encoding = "UTF-8")
df <- fromJSON(response, flatten = TRUE) %>% 
  data.frame()
as_tibble(df) %>% 
  mutate(gene_model = gsub("KH2012:", "", gene_model)) -> all_genes

  #df$genes[[1]] %>% as_tibble()


get_api_genes()
path <- "https://www.aniseed.cnrs.fr/api/all_genes:?organism_id=464"
request <- GET(url = path, 
               query = list(
                 cell = cell,
                 organism_id = 464)
)

get_territories <- function(KH_id){
  #KH_id %>% as_tibble() %>% mutate(KH_id = gsub("\\.v.+", "", value))
  gene_names <- all_genes %>% dplyr::filter(gene_model == KH_id) %>% 
    dplyr::select(gene_name) %>% 
    mutate(gene_name = gsub("\\;.*","",gene_name))
  
  path <- "https://www.aniseed.cnrs.fr/api/all_territories_by_gene:?"
  request <- GET(url = path, 
                 query = list(
                   gene = gene_names,
                   organism_id = 464)
  )
  
  request$status_code
  response <- content(request, as = "text", encoding = "UTF-8")
  df <- fromJSON(response, flatten = TRUE) %>% 
    data.frame()
  as_tibble(df) %>% 
    mutate(KHid = KH_id)
}

hyb_3v4 <- rownames_to_column(topTags(et.hyb.3v4, n=nrow(et.hyb.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)
mocu_3v4 <- rownames_to_column(topTags(et.mocu.3v4, n=nrow(et.mocu.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)
mocc_3v4 <- rownames_to_column(topTags(et.mocc.3v4, n=nrow(et.mocc.3v4$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)

hyb_4v6 <- rownames_to_column(topTags(et.hyb.4v6, n=nrow(et.hyb.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)
mocu_4v6 <- rownames_to_column(topTags(et.mocu.4v6, n=nrow(et.mocu.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)
mocc_4v6 <- rownames_to_column(topTags(et.mocc.4v6, n=nrow(et.mocc.4v6$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)

hyb_3v6 <- rownames_to_column(topTags(et.hyb.3v6, n=nrow(et.hyb.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)
mocu_3v6 <- rownames_to_column(topTags(et.mocu.3v6, n=nrow(et.mocu.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)
mocc_3v6 <- rownames_to_column(topTags(et.mocc.3v6, n=nrow(et.mocc.3v6$table))$table, "gene") %>%
  as_tibble() %>%
  separate(gene, c("gene", "KH.id"), sep = "@") %>%
  mutate(KH.id=gsub("\\.v.+", "", KH.id)) %>% 
  filter(FDR < 0.1)


#test_KHid <- test_KHid %>% as_tibble() %>% mutate(test_KHid = gsub("\\.v.+", "", value))
#res_removeRep2.g_0.05a_territories <- lapply(res_removeRep2.g_0.05a$KHid, get_territories)
territories <- hyb_3v6 %>%
  filter(logFC < 0) %>%
  filter(KH.id != 'NA')
territories <- lapply(territories$KH.id, get_territories)
#mocu_3v6_territories <- territories
hyb_3v6_Dterritories <- territories
#mocu3v4upgenes_territories <- lapply(mocu3v4upgenes$KH.id, get_territories)


territories %>% {
  tibble(
    KHid = map(., "KHid"),
    anatomical_territory = map(., "anatomical_territory"),
    stage = map(., "stage")
  )
} -> territories_

svg("hyb_3v6down_territories.svg")
territories_ %>% unnest(cols = c(anatomical_territory, stage)) %>% 
  #filter(stage > "Stage 22") %>% 
  mutate(anatomical_territory = strsplit(as.character(anatomical_territory), ",")) %>% 
  unnest(anatomical_territory) %>%
  group_by(anatomical_territory) %>% 
  count() %>%
  filter(n > 0) %>% 
  ggplot(aes(x=reorder(anatomical_territory, n), y = n, fill = anatomical_territory)) + 
  geom_bar(stat = "identity") + theme_minimal()+coord_flip() + theme(legend.position='none')
dev.off()


test_territories_counts <- map(test_territories, 1) %>% base::unlist() %>% as_tibble() %>% 
  mutate(value = strsplit(as.character(value), ",")) %>% 
  unnest(value) %>% 
  group_by(value) %>% 
  count(value) %>% View()

ggplot(test_territories_counts, aes(x=value, fill = value)) + geom_bar() + theme_minimal()

test_territories_ <- test_territories %>% {
  tibble(
    anatomical_territory = map(., "anatomical_territory"),
    stage = map(., "stage")
  )
}