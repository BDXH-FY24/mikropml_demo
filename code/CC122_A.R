library(mikropml)
library(tidyverse)
library(ggtext)
library(broom)
rm(list = ls())
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
 

# metadata <-  read_tsv("raw_data/baxter.metadata.tsv",
#          col_types = cols(sample = col_character())) %>% 
#   rename_all(tolower)  %>% 
#   rename(group = sample) %>% 
#   
#   mutate(srn = dx_bin == "Adv Adenoma" | dx_bin == "Cancer",
#          lesion = dx_bin == "Adv Adenoma" | dx_bin == "Cancer" | dx_bin == "Adenoma")   
#         
 
metadata <-  read_tsv("raw_data/baxter.metadata.tsv",
                      col_types = cols(sample = col_character())) %>% 
  mutate(group = sample) %>% 
  rename_all(tolower)  %>% 
  # rename(group = sample) %>% 
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




###

sig_genera <- composite %>% 
  select(-sample) %>% 
  nest(data = -taxonomy) %>% 
  mutate(test = map(.x = data, ~wilcox.test(rel_abund ~ srn, data =.x) %>% tidy)) %>% 
  unnest(test) %>% 
  mutate(p.adjust = p.adjust(p.value,method = "BH")) %>% 
  filter(p.adjust < 0.05) %>% 
  select(taxonomy, p.adjust)


composite %>% 
  inner_join(sig_genera, by = "taxonomy") %>% 
  mutate(rel_abund = 100* rel_abund + 1/20000,
         taxonomy = str_replace(taxonomy, "(.*)", "*\\1*"),
         taxonomy = str_replace_all(taxonomy, "\\*(.*)_unclassified", "Unclassified<br>*\\1"),
         srn = factor(srn, levels =c (T,F))) %>%
ggplot(aes(x = rel_abund, y=taxonomy, color = srn, fill = srn))+
  
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                        jitter.width = 0.3),
              shape = 21)+
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange",
               position = position_dodge(width = 0.8),
               color = "black", show.legend = F)+
  scale_x_log10()+
  theme_minimal()+
  scale_color_manual(NULL,
                     breaks = c(F, T),
                     values = c("gray", "dodgerblue"),
                     labels = c("healthy", "SRN"))+
  scale_fill_manual(NULL,
                     breaks = c(F, T),
                     values = c("gray", "dodgerblue"),
                     labels = c("healthy", "SRN"))+
  labs(x = "relative abundance (%)", y = NULL)+
  theme(
    axis.text.y = element_markdown()
  )
  










ggsave("figures/jitter.png")


















  
  