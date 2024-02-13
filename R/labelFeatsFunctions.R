
labelFeatsManual <- function(row_data){
  mzbounds <- pmppm(row_data$mzmed, 10)
  rtbounds <- row_data$rtmed+c(-1, 1)
  eic <- msdata$MS1[mz%between%mzbounds][rt%between%rtbounds] %>% 
    arrange(rt)
  dev.new(width=6, height=4)
  plot.new()
  plot.window(xlim=rtbounds, ylim=c(0, max(eic$int)), mar=c(1.1, 0.1, 0.1, 0.1))
  for(j in unique(eic$filename)){
    lines(eic[eic$filename==j, c("rt", "int")])
  }
  axis(1)
  title(paste(row_data$id, round(mzbounds[1], 7)))
  abline(v=rtbounds, col="red")
  keyinput <- getGraphicsEvent(prompt = "", onKeybd = function(x){
    return(x)
  })
  dev.off()
  if(keyinput=="ctrl-["){ # aka Esc button
    break
  }
  feat_class <- switch(
    keyinput,
    "Right" = "Good",
    "Left" = "Bad",
    "Up" = "Ambiguous",
    "Down" = "Stans only"
  )
  # print(feat_class)
}

# feat_class <- rep(NA, nrow(feat_feats))
# names(feat_class) <- feat_feats$feature
# # feat_class["FT038"] <- "Good"
# # feat_class["FT001"] <- "Bad"
# while(TRUE){
#   feat_i <- sample(feat_feats$feature[is.na(feat_class)], 1)
#   feat_data_i <- feat_data[feat_data$id==feat_i,]
#   feat_class[feat_data_i$id] <- getPeakClass(feat_data_i)
#   
#   model_df <- feat_feats %>%
#     bind_cols(feat_class=feat_class) %>%
#     mutate(feat_class=case_when(
#       feat_class=="Good"~1,
#       feat_class=="Bad"~0,
#       feat_class=="Meh"~0.5,
#       feat_class=="Ambiguous"~NA
#     )) %>%
#     drop_na()
#   glmodel <- glm(formula=feat_class~med_cor+med_snr+PC1, data = model_df, family = binomial)
#   model_df$pred <- predict(glmodel, newdata=model_df, type = "response")
#   
#   print(table(model_df$feat_class, model_df$pred>0.5))
#   feat_feats %>%
#     mutate(pred_val=predict(glmodel, newdata=feat_feats, type = "response")) %>%
#     with(hist(pred_val, breaks=100)) %>%
#     plot(xlim=c(0, 1))
# }


labelFeatsLasso <- function(){
  
}