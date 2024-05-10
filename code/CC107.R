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
                               levels = c("NonDiarrhealControl",
                                          "DiarrhealControl",
                                          "Case"
                               )))

  
  taxon_rel_abund <- otu_rel_abund %>% 
    filter(level == "phylum") %>% 
    group_by(disease_stat, sample_id, taxon) %>% 
    summarize(rel_abund  = 100 * sum(rel_abund), .groups = "drop") %>% 
    mutate(taxon = str_replace_all(taxon,
                                   "(.*)_unclassified", "Unclassified<br> *\\1*"),
           taxon = str_replace(taxon, "^(\\S*)$", "*\\1*")
    )
  
  taxon_pool <- taxon_rel_abund %>%
    group_by(disease_stat,taxon) %>% 
    summarize(mean = mean(rel_abund), .groups = "drop") %>% 
    group_by(taxon) %>% 
    summarize(pool = max(mean) < 3,
              mean =mean(mean),
              .groups = "drop")

sample_order <- taxon_rel_abund %>% 
  filter(taxon == "*Firmicutes*") %>% 
  arrange(desc(rel_abund)) %>% 
  mutate(order = 1:nrow(.)) %>% 
  select(sample_id, order)
  

pretty <-  c("NonDiarrhealControl" = "Healthy",
             "DiarrhealControl" = "Diarreal,<br>*C. difficile*<br>negative",
             "Case" = "Diarreal,<br>*C. difficile*<br>positive")

inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>% 
    mutate(taxon = if_else(pool, "Other", taxon)) %>% 
    group_by(sample_id, disease_stat, taxon) %>% 
    summarize(rel_abund = sum(rel_abund),
              mean = min(mean),
              .groups = "drop") %>% 
    mutate(taxon = factor(taxon),
           taxon = fct_reorder(taxon, mean, .desc = TRUE),
           taxon =fct_shift(taxon,n=1)) %>%
  inner_join(., sample_order, by = "sample_id") %>% 
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>% 
    ggplot(aes(x =sample_id, y = rel_abund, 
               fill = taxon))+

    geom_col(width = 1)+
    
    scale_fill_manual(name = NULL,
                      breaks = c("*Bacteroidetes*", "*Firmicutes*",
                                 "*Proteobacteria*", "*Verrucomicrobia*",
                                 "Other"),
                      values = c(brewer.pal(4, "Dark2"), "gray")
                      )+
    # scale_x_discrete(breaks = c("NonDiarrhealControl",
    #                             "DiarrhealControl",
    #                             "Case"),
    #                  labels = c("Healthy",
    #                             "Diarreal,<br>*C. difficile*<br>negative",
    #                             "Diarreal,<br>*C. difficile*<br>positive"),
    #                  position = "top")+
    # scale_fill_gradient(name = "Mean<br>Relative Abundance",
    #                     low = "#FFFFFF",
    #                     high = "#FF0000",
    #                     expand = c(0,0),
    #                     limits = c(0,NA))+
  scale_y_continuous(expand = c(0,0))+
  # facet_wrap(~disease_stat,nrow = 1, scales = "free_x")+
  facet_grid(~disease_stat, 
             scales = "free_x", 
             space = "free", 
             switch = "x",
             labeller = labeller(disease_stat = pretty) )+
    labs(x = NULL,
           y = "Mean Relative Abundance")+
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_markdown(),
      legend.title = element_markdown(size = 8, hjust = 0.5),
      legend.title.position = "top",
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.height = unit(12, "pt")
      
      
      
      
      # legend.text = element_markdown(),
      # legend.key.size = unit(10, "pt")
    ) #+
    #coord_fixed(ratio = 0.3)

 
  

  


ggsave("figures/cc106_facet_wrap.png", width = 5, height = 5)


















