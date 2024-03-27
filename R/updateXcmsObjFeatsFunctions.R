

logModelFeatQuality <- function(peak_data, feature_labels, 
                                log_formula=feat_class~med_cor+med_snr,
                                likelihood_threshold=0.5, 
                                verbosity=2){
  # Todo: add checks for formula vars existing in feat_metrics
  model_df <- peak_data %>%
    group_by(feature) %>%
    summarise(mzmed=unique(feat_mzmed), rtmed=unique(feat_rtmed)) %>%
    left_join(feat_metrics, by="feature") %>%
    mutate(feat_class=class_labels)
  if(verbosity>1){
    print(ggplot(model_df, aes(x=rtmed, y=mzmed, color=feat_class, label=feature)) + geom_point())
    print(ggplot(model_df, aes(x=med_cor, y=med_snr, color=feat_class, label=feature)) + geom_point())
    print(ggplot(model_df, aes(x=PC1, y=PC2, color=feat_class, label=feature)) + geom_point())
  }
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
  
  keep_feats <- model_df %>%
    mutate(pred=predict(glmodel, newdata=model_df, type = "response")) %>%
    filter(pred>likelihood_threshold) %>%
    pull(feature)
}

updateXcmsObjFeats(xcms_obj, peak_data, feature_labels,
                   likelihood_threshold=0.5, verbosity=2){
  good_features <- logModelFeatQuality(peak_data, feature_labels, 
                                       likelihood_threshold = likelihood_threshold, 
                                       verbosity = verbosity)
  featureDefinitions(xcms_obj) <- featureDefinitions(xcms_obj)[good_features,]
  xcms_obj
}