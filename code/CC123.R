source("code/CC123_A.R")



get_sens_spec <- function(threshold, score, actual, direction){
  
  # threshold <- 100
  # score <- test$score
  # actual <- test$srn
  # direction <- "greaterthan"
  
  predicted <- if(direction == "greaterthan"){
    score > threshold
 } else{
      score < threshold 
    }
  
  tp <- sum(predicted & actual)
  tn <- sum(!predicted & !actual)
  fp <- sum(predicted & !actual)
  fn <- sum(!predicted & actual)
  
  specificity <- tn/(tn+fp)
  sensitivity <- tp/(tp+fn)
  tibble("specificity" = specificity, "sensitivity" = sensitivity)
  
}


# Run fucntion  
get_sens_spec(100, test$score,test$srn, "greaterthan" )


test <- composite %>% 
  inner_join(sig_genera, by = "taxonomy") %>%
  select(group, taxonomy, rel_abund, fit_result, srn) %>% 
  pivot_wider(names_from = taxonomy, 
              values_from = rel_abund) %>% 
  pivot_longer(cols = -c(group, srn), 
               names_to = "metric", 
               values_to = "score") %>% 
  filter(metric =="fit_result")


  
  
 