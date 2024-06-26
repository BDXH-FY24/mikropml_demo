#### CC107
library(tidyverse)
library(readxl)
library(ggtext)
library(glue)
library(RColorBrewer)


set.seed(1234)


metadata <- read_excel(path="raw_data/schubert.metadata.xlsx", na="NA") %>% 
  select(sample_id, disease_stat) %>% 
  drop_na(disease_stat)

otu_counts <- read_tsv("raw_data/schubert.subsample.shared") %>% 
  select(Group, starts_with("Otu")) %>% 
  rename(sample_id = Group) %>% 
  pivot_longer(-sample_id, names_to = "otu", values_to = "count")

taxonomy <- read_tsv("raw_data/schubert.cons.taxonomy") %>% 
  select("OTU", "Taxonomy") %>% 
  rename_all(tolower) %>% 
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace_all(taxonomy, ";$", "")) %>% 
  separate(taxonomy,
           into = c("kingdom", "phylum", "class", "order",
                    "family", "genus"),
           sep = ";")

otu_rel_abund <- inner_join(metadata, otu_counts, by = "sample_id") %>% 
  inner_join(., taxonomy, by = "otu") %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = count / sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  pivot_longer(c("kingdom","phylum", "class", "order",
                 "family", "genus", "otu"),
               names_to = "level", 
               values_to = "taxon") %>% 
  mutate(disease_stat = factor(disease_stat,
                               levels = c("Case",
                                          "DiarrhealControl",
                                          "NonDiarrhealControl"
                               )))

  
taxon_rel_abund <- otu_rel_abund %>% 
    filter(level == "genus") %>% 
    group_by(disease_stat, sample_id, taxon) %>% 
    summarize(rel_abund = sum(rel_abund), .groups = "drop") %>% 
    group_by(disease_stat, taxon) %>% 
    summarize(mean_rel_abund  = 100 * mean(rel_abund), .groups = "drop") %>% 
    mutate(taxon = str_replace_all(taxon,
                                   "(.*)_unclassified", "Unclassified<br>*\\1*"),
           taxon = str_replace(taxon, "^([^<]*)$", "*\\1*"),
           taxon = str_replace_all(taxon, "_"," ")
    )
  
taxon_pool <- taxon_rel_abund %>%
    group_by(taxon) %>% 
    summarize(pool = max(mean_rel_abund) < 1,
              mean =mean(mean_rel_abund) ,
              .groups = "drop")
  
inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>% 
    mutate(taxon = if_else(pool, "Other", taxon)) %>% 
    group_by(disease_stat, taxon) %>% 
    summarize(mean_rel_abund = sum(mean_rel_abund),
              mean = min(mean),
              .groups = "drop") %>% 
    mutate(taxon = factor(taxon),
           taxon = fct_reorder(taxon, mean, .desc = F)) %>% 
    ggplot(aes(y =taxon, x = mean_rel_abund, 
               fill = disease_stat))+
    geom_col(position = position_dodge(), width = 0.6)+
    scale_fill_manual(name=NULL,
                      breaks = c("NonDiarrhealControl",
                                "DiarrhealControl",
                                "Case"),
                     labels = c("Healthy",
                                "Diarreal,<br>*C. difficile* negative",
                                "Diarreal,<br>*C. difficile* positive"),
                     values = c("gray", "dodgerblue", "red")
                     )+
  scale_x_continuous(expand = c(0, NA))+
    labs(y = NULL,
           x = " Mean Relative Abundance (%)"
         )+
    theme_classic() +
    theme(
      axis.text.y = element_markdown(hjust = 1, vjust =1, size = 6),
      legend.text = element_markdown(size = 5),
      legend.position = c(0.8,0.6), 
      legend.key.height = unit(10, "pt"),
      panel.grid.major.x = element_line(color = "lightgray")
    )

 
ggsave("figures/cc107_dodgermap.png", width = 5, height = 6)


















