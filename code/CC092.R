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

# set.seed(1234)
metadata_alpha %>% 
  # ggplot(aes(x = disease_stat, y = mean, ymin = min, ymax = max, fill = disease_stat))+
  # geom_errorbar(width = 0.3)+
  # geom_col(show.legend = F)+
  ggplot(aes(x = disease_stat, 
             y =invsimpson, 
             fill = disease_stat)) +
  
  # geom_violin(show.legend = F,alpha = 0.55)+
  geom_dotplot(binaxis = "y",
               binwidth = 0.35,
               stackdir = "center" ,
               show.legend = F)+
  stat_summary(fun = median, show.legend = F,
               geom = "crossbar",width = 0.3,
               color = "black")+
  geom_boxplot(show.legend = F,outlier.shape = NA,
               alpha = 0.25, width = .6,
               coef = 1)+ # Coef: by default 1.5 of IQR
  # geom_jitter(show.legend = F,
  #             color = "black",
  #             shape = 21,
  #               position = position_jitter(width = 0.25, seed = 1234))+
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

ggsave("figures/StatSummary-cc092.png", width = 4.5, height = 3.5)


















