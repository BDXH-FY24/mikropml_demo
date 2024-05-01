library(tidyverse)
library(readxl)
library(ggtext)
library(glue)

set.seed(1234)
metadata <- read_excel(path="raw_data/schubert.metadata.xlsx", na="NA")

alpha_diversity <- read_tsv("raw_data/schubert.groups.ave-std.summary") %>%
  filter(method == "ave") %>%
  select(-label, -method)

metadata_alpha <- inner_join(metadata, alpha_diversity, by=c('sample_id'='group')) %>%
  mutate(disease_stat = factor(disease_stat,
                               levels=rev(c("NonDiarrhealControl",
                                        "DiarrhealControl",
                                        "Case")))
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

kt <- kruskal.test(invsimpson ~ disease_stat, data = metadata_alpha)
pt <- pairwise.wilcox.test(metadata_alpha$invsimpson, 
                     g=metadata_alpha$disease_stat,
                     p.adjust.method = "BH")


if(kt$p.value < 0.05){
  pt <- pairwise.wilcox.test(metadata_alpha$invsimpson, 
                             g=metadata_alpha$disease_stat,
                             p.adjust.method = "BH")
} #else{}



strip_chart <- metadata_alpha %>% 
  ggplot(aes(x = disease_stat, 
             y =invsimpson, 
             fill = disease_stat)) +
  stat_summary(fun = median, show.legend = F,
               geom = "crossbar")+
  geom_jitter(show.legend = F,
              color = "black",
              width = 0.25,
              shape = 21)+
  labs(x = NULL, y = "Inverse simpson index" )+
  scale_x_discrete(#guide = guide_axis(angle = 40),
                   breaks = c("NonDiarrhealControl",
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
  theme(axis.text.y = element_markdown(hjust = 1)) +
  coord_flip()

strip_chart + 
  geom_line(data = tibble(x =c(2,3), y=c(23,23)), 
            aes(x = x, y =y),
            inherit.aes = F)+
  geom_line(data = tibble(x =c(1,2.5), y=c(33,33)), 
            aes(x = x, y =y),
            inherit.aes = F)+
  geom_text(data = tibble(x = 2.5, y=23.5), 
            aes(x = x, y =y, label = "n.s"),
            size = 6,
            inherit.aes = F)+
  geom_text(data = tibble(x = 1.75, y=33), 
            aes(x = x, y =y, label = "*"),
            size = 10,
            inherit.aes = F)
  
  

ggsave("figures/StatSummary-cc093.png", width = 4.5, height = 3.5)


















