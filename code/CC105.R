#### CC105
library(tidyverse)
library(readxl)
library(ggtext)
library(glue)
library(brew)

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
    summarize(rel_abund = sum(rel_abund), .groups = "drop") %>% 
    group_by(disease_stat, taxon) %>% 
    summarize(mean_rel_abund  = 100 * mean(rel_abund), .groups = "drop") %>% 
    mutate(taxon = str_replace_all(taxon,
                                   "(.*)_unclassified", "Unclassified<br> *\\1*"),
           taxon = str_replace(taxon, "^(\\S*)$", "*\\1*")
    )
  
  taxon_pool <- taxon_rel_abund %>%
    group_by(taxon) %>% 
    summarize(pool = max(mean_rel_abund) < 3,
              mean =mean(mean_rel_abund),
              .groups = "drop")
  
  htmap <- inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>% 
    mutate(taxon = if_else(pool, "Other", taxon)) %>% 
    group_by(disease_stat, taxon) %>% 
    summarize(mean_rel_abund = sum(mean_rel_abund),
              mean = min(mean),
              .groups = "drop") %>% 
    mutate(taxon = factor(taxon),
           taxon = fct_reorder(taxon, mean, .desc = FALSE) )
  
 htmap %>% 
           
    ggplot(aes(x =disease_stat, fill = mean_rel_abund, 
               y = taxon))+
   geom_tile()
    # geom_col()+
    geom_tile(show.legend = FALSE)+
    geom_text(aes(label = format(round( mean_rel_abund, 1), nsamll = 1)) )+
    # scale_fill_manual(name = NULL,
    #                   breaks = c("*Bacteroidetes*", "*Firmicutes*",
    #                              "*Proteobacteria*", "*Verrucomicrobia*",
    #                              "Other"), 
    #                   values = my36colors)+
    scale_x_discrete(breaks = c("NonDiarrhealControl",
                                "DiarrhealControl",
                                "Case"),
                     labels = c("Healthy",
                                "Diarreal,<br>*C. difficile*<br>negative",
                                "Diarreal,<br>*C. difficile*<br>positive"),
                     position = "top")+
    scale_fill_gradient(name = "Mean<br>Relative Abundance",
                        low = "#FFFFFF",
                        high = "#FF0000",
                        expand = c(0,0),
                        limits = c(0,NA))+
    labs(x = NULL,
           y = NULL)+
    theme_classic() +
    theme(
      axis.text.x.top = element_markdown(vjust = 0.5),
      axis.text.y = element_markdown(),
      legend.title = element_markdown(size = 8, hjust = 0.5),
      legend.title.position = "top",
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.height = unit(12, "pt")
      
      
      
      
      # legend.text = element_markdown(),
      # legend.key.size = unit(10, "pt")
    )+
    coord_fixed(ratio = 0.3)

 
  

  


ggsave("figures/cc105_heatmap.png", width = 5, height = 5)


















