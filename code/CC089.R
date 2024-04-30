#### CC08

library(tidyverse)
library(readxl)
library(ggtext)
library(glue)

metadata <- read_excel(path="raw_data/schubert.metadata.xlsx", na="NA")

alpha_diversity <- read_tsv("raw_data/schubert.groups.ave-std.summary") %>%
  filter(method == "ave") %>%
  select(-label, -method)

metadata_alpha <- inner_join(metadata, alpha_diversity, by=c('sample_id'='group')) %>%
  mutate(disease_stat = factor(disease_stat,
                               levels=c("NonDiarrhealControl",
                                        "DiarrhealControl",
                                        "Case"))
  )

healthy_color <- "#BEBEBE"
diarrhea_color <- "#0000FF"
case_color <- "#FF0000"



invsimpson_summary <- metadata_alpha %>% 
  group_by(disease_stat) %>% 
  summarize(mean = mean(invsimpson),
            se = sd(invsimpson) / sqrt(n()),
            max =mean+se,
            min=mean-se,
            N=n(),
            .groups = "drop") 

healthy_n <- invsimpson_summary %>% 
  filter(disease_stat == "NonDiarrhealControl") %>% 
  pull(N)



diarrhea_n <- invsimpson_summary %>% 
  filter(disease_stat == "DiarrhealControl") %>% 
  pull(N)


case_n <- invsimpson_summary %>% 
  filter(disease_stat == "Case") %>% 
  pull(N)


metadata_alpha %>% 
  # ggplot(aes(x = disease_stat, y = mean, ymin = min, ymax = max, fill = disease_stat))+
  # geom_errorbar(width = 0.3)+
  # geom_col(show.legend = F)+
  ggplot(aes(x = disease_stat, y =invsimpson, color = disease_stat)) +
  stat_summary(fun.data = median_hilow, 
               geom = 'pointrange', 
               show.legend = F,
               color = "gray",
               group = 1)+
  stat_summary(fun.data = median_hilow,
               geom = 'point',
               show.legend = F,
               size = 3)+
  # stat_summary(fun.data = median_hilow,
  #              geom = 'line',
  #              show.legend = F,
  #              group = 1,
  #              color = "black")+
  
  
  
  labs(x = NULL, y = "Inverse simpson index" )+
  scale_x_discrete(breaks = c("NonDiarrhealControl",
                              "DiarrhealControl",
                              "Case"),
                   labels = c(glue("Health<br>(N={healthy_n})"), 
                   glue("Diarrhea and<br>*C.difficile negative*<br>(N={diarrhea_n})"),
                   glue("Diarrhea and<br>*C.difficile positive*<br>\\
                              (N={case_n})")))+
  scale_color_manual(name = NULL,
                     breaks=c("NonDiarrhealControl",
                              "DiarrhealControl",
                              "Case"),
                     labels = c("Health", "Diarrhea and<br>*C.difficile negative*",
                                "Diarrhea and<br>*C.difficile positive*"),
                     values = c(healthy_color, diarrhea_color,case_color))+
  scale_fill_manual(name = NULL,
                     breaks=c("NonDiarrhealControl",
                              "DiarrhealControl",
                              "Case"),
                     labels = c("Health", "Diarrhea and<br>*C.difficile negative*",
                                "Diarrhea and<br>*C.difficile positive*"),
                    values = c(healthy_color, diarrhea_color,case_color))+
  theme_minimal()+
  theme(axis.text.x = element_markdown())

ggsave("figures/diversity.tiff", width = 4.5, height = 3.5)


















