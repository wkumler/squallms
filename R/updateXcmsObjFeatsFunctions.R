
logModelFeatProb <- function(feature_metrics, feature_labels, 
                             log_formula=feat_class~med_cor+med_snr,
                             verbosity=2){
  # Todo: add checks for formula vars existing in feat_metrics
  model_df <- feature_metrics %>%
    left_join(data.frame(feat_class=feature_labels) %>% rownames_to_column("feature"),
              by="feature")
  glmodel <- model_df %>%
    filter(!is.na(feat_class)) %>%
    mutate(feat_class=ifelse(feat_class=="Bad", 0, 1)) %>%
    glm(formula=log_formula, family = binomial)
  if(verbosity>0){
    message("Logistic model regression coefficients:")
    print(glmodel$coefficients)
  }
  
  if(verbosity>1){
    gp <- model_df %>%
      mutate(pred=predict(glmodel, newdata=model_df, type = "response")) %>%
      mutate(feat_class=ifelse(is.na(feat_class), "Unclassified", feat_class)) %>%
      ggplot(aes(x=med_cor, y=med_snr, color=pred, shape=feat_class, label=feature)) +
      geom_point(size=3) +
      scale_shape_manual(breaks=c("Good", "Unclassified", "Bad"),
                         values=c(19, 1, 4)) +
      scale_color_gradientn(colors = c("#a41118", "#f4bb23", "#028e34"), limits=c(0, 1)) +
      theme_bw() +
      labs(x="Peak shape coefficient", y="Signal-to-noise ratio", shape="Lasso class",
           color="Predicted\nclass")
    print(gp)
  }
  
  preds <- predict(glmodel, newdata=model_df, type = "response")
  names(preds) <- model_df$feature
  preds
}

logModelFeatQuality <- function(feature_metrics, feature_labels, 
                                log_formula=feat_class~med_cor+med_snr,
                                likelihood_threshold=0.5, 
                                verbosity=2){
  feature_preds <- logModelFeatProb(feature_metrics, feature_labels, 
                                    log_formula=log_formula,
                                    verbosity=verbosity)
  pred_df <- feature_metrics %>%
    left_join(data.frame(feat_class=feature_labels) %>% rownames_to_column("feature"),
              by="feature") %>%
    left_join(data.frame(pred_prob=feature_preds) %>% rownames_to_column("feature"),
              by="feature") %>%
    mutate(pred_class=ifelse(pred_prob>likelihood_threshold, "Good", "Bad"))
  
  if(verbosity>0){
    print(caret::confusionMatrix(factor(pred_df$pred_class), factor(pred_df$feat_class)))
  }
  
  keep_feats <- pred_df %>%
    filter(pred_class=="Good") %>%
    pull(feature)
}

updateXcmsObjFeats <- function(xcms_obj, feature_metrics, feature_labels,
                               likelihood_threshold=0.5, verbosity=2){
  good_features <- logModelFeatQuality(feature_metrics, feature_labels, 
                                       likelihood_threshold = likelihood_threshold, 
                                       verbosity = 2)
  featureDefinitions(xcms_obj) <- featureDefinitions(xcms_obj)[good_features,]
  
  xcms_obj@.processHistory[[length(xcms_obj@.processHistory)+1]] <- list(
    type="Low-quality feature removal via squallms",
    date=date(),
    info="Nonstandard processing step",
    fileIndex=seq_along(fileNames(xcms_filled)),
    `Parameter class`="none",
    `MS level(s)`=1
  )
  # ph <- xcms:::ProcessHistory(date. = date(),
  #                             type. = "Blah",
  #                             fileIndex = 1:length(fileNames(xcms_filled)))
  # xph <- xcms:::XProcessHistory(param = NULL, date. = date(), 
  #                              type. = "Unknown", 
  #                              fileIndex = 1:length(fileNames(xcms_filled)), 
  #                              msLevel = 1)
  # xcms:::addProcessHistory(xcms_filled, ph) %>%
  #   processHistory()
  xcms_obj
}
