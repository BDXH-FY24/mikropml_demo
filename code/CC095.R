#### CC095
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

kt <- kruskal.test(sobs ~ disease_stat, data = metadata_alpha)
pt <- pairwise.wilcox.test(metadata_alpha$sobs, 
                     g=metadata_alpha$disease_stat,
                     p.adjust.method = "BH")


if(kt$p.value < 0.05){
  pt <- pairwise.wilcox.test(metadata_alpha$invsimpson, 
                             g=metadata_alpha$disease_stat,
                             p.adjust.method = "BH")
} #else{}

kt <- kruskal.test(shannon ~ disease_stat, data =metadata_alpha)



 strip_chart <- metadata_alpha %>% 
   select(disease_stat, sobs, shannon, invsimpson) %>% 
   pivot_longer(-disease_stat, names_to = "metric", values_to = "values") %>% 
   mutate(metric = recode(metric, 
                          sobs = "Observed Richness",
                          shannon = "Shannon Diversity",
                          invsimpson = "Inverse Simpson"),
          metric = factor(metric, 
                          levels = c("Observed Richness","Inverse Simpson", "Shannon Diversity"))) %>% 
   ggplot(aes(x = disease_stat, y = values, fill = disease_stat))+
   stat_summary(fun = median, geom="crossbar" , show.legend = F, width = 0.4)+
   geom_jitter(show.legend = F, width = 0.25, shape =21, color = "black")+
   facet_wrap(~metric,nrow = 3, scales = "free_y", strip.position = "left")+ # Cool strip.position to change the location of facet lables
   
  scale_x_discrete(#guide = guide_axis(angle = 40),
                   breaks = c("NonDiarrhealControl",
                              "DiarrhealControl",
                              "Case"),
                   labels = c(glue("Health<br>(N={healthy_n})"), 
                   glue("Diarrhea and<br>*C.difficile negative*<br>(N={diarrhea_n})"),
                   glue("Diarrhea and<br>*C.difficile positive*<br>\\
                              (N={case_n})")))+
  scale_y_continuous(limits = c(0, NA))+
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
  theme(axis.text.x = element_markdown(),
        strip.placement = "outside",
        strip.background = element_rect(colour = NA)) 
   

line <- tibble(
  metric = factor(c( "Observed Richness",
              "Observed Richness",
              "Shannon Diversity",
              "Shannon Diversity",
              "Inverse Simpson",
              "Inverse Simpson"),
              levels = c("Observed Richness",
                         "Inverse Simpson",
                         "Shannon Diversity"
                         )),
  x = c(2,1,2,1,2,1),
  xend = c(3,2.5,3,2.5,3,2.5),
  y = c(125,190,4,4.5,23,33),
  yend = y
) 


stars <- tibble(
  metric = factor(c( "Observed Richness",
                     "Observed Richness",
                     "Shannon Diversity",
                     "Shannon Diversity",
                     "Inverse Simpson",
                     "Inverse Simpson"),
                  levels = c("Observed Richness",
                             "Inverse Simpson",
                             "Shannon Diversity"
                  )),
  x = c(2.5,1.75,2.5,1.75,2.5,1.75),
  xend = c(3,2.5,3,2.5,3,2.7),
  y = c(131.5,195,4.2,4.6,24.5,33.25),
  yend = y,
  label =c("n.s","*", "n.s","*","n.s","*")
) 

  
strip_chart + 
  geom_segment(data = line, aes(x = x, y = y, 
                                xend =xend, yend = yend
                                ),
               inherit.aes = F)+
  geom_text(data = stars, aes(x=x, y=y, label =label),
            inherit.aes = F)
  # geom_text(data = tibble(x = 1.5, y=25), 
  #           aes(x = x, y =y, label = "n.s"),
  #           size = 6,
  #           inherit.aes = F)+
  # geom_text(data = tibble(x = 2.25, y=33.6), 
  #           aes(x = x, y =y, label = "*"),
  #           size = 10,
  #           inherit.aes = F)
  
  

ggsave("figures/cc095.png", width = 9, height = 7.5)


















