library(mikropml)
library(tidyverse)
library(ggtext)
library(broom)
library(glue)
### using stringR to replace and clean the data
#### code 121: https://www.youtube.com/watch?v=tBxGVfvx-Gc&list=PLmNrK_nkqBpIIRdQTS2aOs5OD7vVMKWAi&index=41&t=53s


shared <- read_tsv("raw_data/baxter.subsample.shared",
         col_types = cols(Group = col_character(),
                          .default = col_double())) %>% 
  rename_all(tolower) %>% 
  select(group, starts_with("otu")) %>% 
  pivot_longer(-group, names_to = "otu", values_to = "count")   
  

taxonomy <- read_tsv("raw_data/baxter.cons.taxonomy") %>% 
  rename_all(tolower) %>%
  select(otu, taxonomy) %>% 
  mutate(otu = tolower(otu),
         taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";unclassified", 
                                "_unclassified"),
         taxonomy = str_replace_all(taxonomy, 
                                    ";unclassified", ""),
         taxonomy = str_replace_all(taxonomy, ";$", ""),
         taxonomy = str_replace_all(taxonomy, ".*;", ""))  
 

metadata <-  read_tsv("raw_data/baxter.metadata.tsv",
         col_types = cols(sample = col_character())) %>% 
  rename_all(tolower)  %>% 
  rename(group = sample) %>% 
  mutate(srn = dx_bin == "Adv Adenoma" | dx_bin == "Cancer",
         lesion = dx_bin == "Adv Adenoma" | dx_bin == "Cancer" | dx_bin == "Adenoma")   
        
 
composite <- inner_join(shared, taxonomy, by = c("otu")) %>% 
  group_by(group, taxonomy) %>% 
  summarize(count= sum(count), .groups = "drop") %>% 
  group_by(group) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  inner_join(., metadata, by = "group")




sig_genera <- composite%>% 
  nest(data = -taxonomy) %>%
  mutate(test = map(.x = data, 
                    ~wilcox.test(rel_abund~srn,data = .x) %>% tidy)) %>% 
  unnest(test) %>% 
  filter(p.value < 0.05) %>% 
  mutate(p_adujst = p.adjust(p.value, method = "BH")) %>% 
  filter(p_adujst < 0.05) %>% 
  select(taxonomy, p_adujst)
  
composite %>% 
  inner_join(sig_genera, by ="taxonomy") %>% 
  ggplot(aes(x = rel_abund, y = taxonomy, color = srn))+
  geom_jitter()





 
ggsave("figures/significant_genera.tiff", width = 6, height = 4)
 
  















  
  