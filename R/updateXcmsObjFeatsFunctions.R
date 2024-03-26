

updateXcmsObjFeats <- function(xcms_obj, peak_data, feature_labels,
                               likelihood_threshold=0.5, verbosity=2){
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
    mutate(feat_onehot=ifelse(feat_class=="Bad", 0, 1)) %>%
    filter(!is.na(feat_class)) %>%
    glm(formula=feat_onehot~med_cor+med_snr, family = binomial)
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
  featureDefinitions(xcms_obj) <- featureDefinitions(xcms_obj)[keep_feats,]
  xcms_obj
}
