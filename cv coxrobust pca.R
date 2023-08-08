
#### CV for PCA ROBUST COX ####
library(coxrobust)
library(tidyverse)
library(survival)

liv <- read_csv("C:/Users/kange/Desktop/Thesis/data/liver_clinic.csv")  #keep barcode for merging only, not in model
liv <- liv %>% select(-c(1))
cn <- read_csv("C:/Users/kange/Desktop/Uoft/Uoft Courses/practicum/week13/cn_pc_nolast.csv")
dm <- read_csv("C:/Users/kange/Desktop/Uoft/Uoft Courses/practicum/week13/dm_genomic_pcs.csv")
ge <- read_csv("C:/Users/kange/Desktop/Uoft/Uoft Courses/practicum/week13/ge_pc_nolast.csv")

allpcs <- cn %>% inner_join(dm, by="barcode") %>% inner_join(ge, by="barcode")

liv_final <- liv %>% inner_join(allpcs, by="barcode")

liv_final <- liv_final %>% mutate(stage_recode = ifelse(stage == "Stage I" | stage == "Stage II", "Early Stage", 
                                                        "Not Early Stage"))


repcrossval <- function(repeats, folds, weight){
  cv_conc <- rep(0, repeats)
  for (k in 1:repeats){
    splits <- split(liv_final, sample(1:folds, nrow(liv_final), replace = T))  
    
    c <- rep(0, folds)
    for (i in 1:folds){
      d_minusi <- bind_rows(splits[-c(i)])  # dataset without ith fold 
      fit <- coxr(Surv(OS.time, OS) ~age+pc1_cn+pc5_cn+pc1_dm+pc5_dm+pc4_ge+stage_recode, 
                  data = d_minusi, f.weight = weight) # fit the "training model" 
      
      # predict new data's HR/risk using training model
      somedf <- splits[[i]] %>% mutate(stage_recodeNotEarly=ifelse(stage_recode=="Not Early Stage", 1, 0))
      somedf <- somedf %>% select(c(age, pc1_cn, pc5_cn, pc1_dm, pc5_dm, pc4_ge, stage_recodeNotEarly))
      preds <- data.frame(preds = exp(as.matrix(somedf) %*% fit$coefficients))
      
      # bind preds to testing data
      test <- splits[[i]]
      test_preds <- cbind(test, preds)
      
      # useful info
      surv_time <- test$OS.time
      censor_status <- test$OS
      risk <- test_preds$preds
      
      # call concordance
      
      
      pairs <- data.frame(t(combn(c(1:nrow(test)), 2)))
      colnames(pairs) <- c("obs1", "obs2")
      
      obs1 <- pairs$obs1
      obs2 <- pairs$obs2
      time1 <- surv_time[obs1]
      time2 <- surv_time[obs2]
      status1 <- censor_status[obs1]
      status2 <- censor_status[obs2]
      risk1 <- risk[obs1]
      risk2 <- risk[obs2]
      good_df <- data.frame(cbind(obs1, obs2, time1, time2, status1, status2, risk1, risk2))
      # good_df contains pairs that can be compared
      good_df <- good_df %>% filter(!(status1==0 & status2==0)) %>%  # not both censored
        filter(!(time1<=time2 & status1==0 & status2==1)) %>%   # time1<time2, time1 censored, time2 not 
        filter(!(time1>=time2 & status2==0 & status1==1)) %>%  # time2<time1, time2 censored, time1not
        rowid_to_column()
      
      # tied predictors - just see if any duplications in dataframe, all.vars are the variables in the model
      # duplicated(liv_final[, all.vars(fit$call[[2]])[-c(1:2)]])
      
      # tied observed event times (a_y)
      #tied_obs_ind <- good_df %>% filter(time1==time2) %>% select(rowid)
      #ay <- nrow(tied_obs)
      # tied predictor and event time
      
      # conc vs discordant
      conc_disc <- good_df %>% 
        filter(!(time1==time2)) %>% 
        mutate(conc = case_when(
          (time1 < time2 & risk1 > risk2) ~ 1,
          (time2 < time1 & risk2 > risk1) ~ 1,
          TRUE ~ 0
        ))
      concordance <- mean(conc_disc$conc)
      c[i] <- concordance
      
    }
    cross_val_score <- mean(c)
    cv_conc[k] <- cross_val_score
  }
  mean(cv_conc)
}

# deafult is quadratic