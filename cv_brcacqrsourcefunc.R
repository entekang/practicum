
crossval_cqr <- function(repeats, folds, quantile){
  
  cv_conc <- rep(0, repeats)
  for (k in 1:repeats){
    splits <- split(brca1, sample(1:folds, nrow(brca1), replace = T))
    
    c <- rep(0, folds)
    for (i in 1:folds){
      d_minusi <- bind_rows(splits[-c(i)])  # dataset without ith fold 
      fit <- crq(Surv(pfs_time, pfs_out) ~BRCA_Overall2 + LengthofSurgery_hour+ADJcycles+Ascites+DSscore, 
                 data = d_minusi, method = "PengHuang")    # fit the "training model" 
      
      # predict new data's risk/time using training model
      q1coef <- coef(fit, taus = c(quantile))
      #model matrix for data to be predicted, splits[[i]]
      somedf <- splits[[i]] %>% mutate(BRCA_Overall2Positive = ifelse(BRCA_Overall2=="Positive", 1, 0), 
                                       AscitesYes = ifelse(Ascites=="Yes", 1, 0), DSscoreLow = ifelse(DSscore=="Low", 1, 0), 
                                       DSscoreModerate=ifelse(DSscore=="Moderate", 1, 0))
      somedf <- somedf %>% select(c(BRCA_Overall2Positive, LengthofSurgery_hour, ADJcycles, 
                                    AscitesYes, DSscoreLow, DSscoreModerate))
      int <- rep(1, nrow(splits[[i]]))
      modelmatrix <- cbind(int, somedf)
      preds1 <- data.frame(preds=as.matrix(modelmatrix) %*% q1coef)
      
      # bind preds to testing data
      test <- splits[[i]]
      test_preds <- cbind(test, preds1)
      
      # useful info
      surv_time <- test$pfs_time
      censor_status <- test$pfs_out
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
      risk1 <- risk[obs1]  # predicted survival time
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
          (time1 < time2 & risk1 < risk2) ~ 1,  # need to reverse direction of "risk" now, which means predicted surv time
          (time2 < time1 & risk2 < risk1) ~ 1,
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
