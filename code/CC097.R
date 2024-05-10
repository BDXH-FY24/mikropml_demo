#### CC096
library(tidyverse)
library(readxl)
library(ggtext)
library(glue)



get_sens_spec <- function(x, data){
  predicted <- x > data$invsimpson
  
  tp <- sum(predicted & data$disease_stat)
  tn <- sum(!predicted & !data$disease_stat)
  fp <- sum(predicted & !data$disease_stat)
  fn <- sum(!predicted & data$disease_stat)
  
  specificity <- tn/(tn+fp)
  sensitivity <- tp/(tp+fn)

  tibble("sensitivity" = sensitivity,"specificity" =  specificity)
  
}

set.seed(1234)
metadata <- read_excel(path="raw_data/schubert.metadata.xlsx", na="NA")

alpha_diversity <- read_tsv("raw_data/schubert.groups.ave-std.summary") %>%
  filter(method == "ave") %>%
  select(-label, -method)

disease_invsimpson <- inner_join(metadata,
                                 alpha_diversity,
                                 by = c("sample_id" = "group")) %>% 
  select(disease_stat, invsimpson)

negative <- "NonDiarrhealControl"
positive <- "Case"


get_roc_data <- function(negative, positive){
disease_invsimpson %>% 
  filter(disease_stat == negative | disease_stat == positive ) %>% 
  mutate(disease_stat = recode(disease_stat,
                                "{negative}" := FALSE,
                                "{positive}" := TRUE)) %>% 
  mutate(sens_spec = map_dfr( invsimpson, get_sens_spec,.) )%>% 
  unnest(sens_spec) %>% 
  mutate(comparison = glue("{negative}_{positive}"))
  
}

roc_data <- bind_rows(
get_roc_data("NonDiarrhealControl","DiarrhealControl"),
get_roc_data("NonDiarrhealControl", "Case"),
get_roc_data("DiarrhealControl", "Case")
) %>% 
  arrange(invsimpson)




roc_data %>% 
  separate(comparison, into = c("negative", "positive"), remove = FALSE)%>% 
  mutate(disease_stat = if_else(disease_stat,positive, negative)) %>% 
  
  
  
  ggplot(aes(x = 1-specificity, y = sensitivity, 
             group = comparison,
             color = disease_stat))+
  geom_line(linewidth =1, lineend = "round")+
  geom_abline(slope = 1, intercept = 0, color = "black", alpha = 0.6)+
  theme_minimal()
  


ggsave("figures/cc097.png", width = 7, height = 5)


















