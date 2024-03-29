
#' Model feature quality using a logistic regression
#' 
#' After metrics have been extracted (typically with `extractChromMetrics`) and 
#' labeling has occurred (typically with `labelFeatsManual` or `labelFeatsLasso`),
#' a model can be built to robustly classify the peaks that were not labeled
#' based on the metrics. The default regression equation fits the med_cor and
#' med_snr columns in feature_metrics to the feature_labels, but additional
#' metrics can be passed as well (but be careful of overfitting!).
#'
#' @param feature_metrics A data.frame with columns used to construct the 
#' feature quality model.
#' @param feature_labels A character vector named with feature IDs and entries
#' corresponding to the peak quality (either "Good", "Bad", or NA).
#' @param log_formula The formula to use when predicting feature quality from
#' the feature metrics. This formula is passed to `glm` as-is, so make sure that
#' the predictive features exist in the feature_metrics data.frame.
#' @param verbosity Scalar value between zero and two determining how much 
#' diagnostic information is produced. Zero should return nothing, one should
#' return text-based progress markers, and 2 will return diagnostic plots.
#'
#' @return A numeric vector of probabilities returned by the logistic model
#' named by feature ID
#' @export
#'
#' @examples
#' msnexp_filled <- readRDS(system.file("extdata", "intro_xcms_filled.rds", package="squallms"))
#' peak_data <- makeXcmsObjFlat(msnexp_filled)
#' feat_metrics <- extractChromMetrics(peak_data, verbosity = 0)
#' lasso_classes <- readRDS(system.file("extdata", "intro_lasso_labels.rds", package="squallms"))
#' feat_probs <- logModelFeatProb(feat_metrics, lasso_classes)
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

#' Turn 0-1 likelihood values into categorical (good/bad) classifications
#' 
#' This function wraps `logModelFeatProb` for modeling of feature quality
#' to estimate the quality of every feature in the dataset and classify
#' them as good/bad based on whether they exceed the provided 
#' `likelihood_threshold`. Lower likelihood thresholds (0.01, 0.1) will produce
#' more false positives (noise peaks included when they shouldn't be). Higher 
#' likelihood thresholds (0.9, 0.99) will produce more false negatives (good
#' peaks removed when they shouldn't be).
#'
#' @inheritParams logModelFeatProb
#' @param likelihood_threshold A scalar numeric above which features will be
#' kept if their predicted probability exceeds. 
#'
#' @return A character vector of feature quality assessments returned by the 
#' logistic model and named by feature ID 
#' @export
#'
#' @examples
#' msnexp_filled <- readRDS(system.file("extdata", "intro_xcms_filled.rds", package="squallms"))
#' peak_data <- makeXcmsObjFlat(msnexp_filled)
#' feat_metrics <- extractChromMetrics(peak_data, verbosity = 0)
#' lasso_classes <- readRDS(system.file("extdata", "intro_lasso_labels.rds", package="squallms"))
#' feat_classes <- logModelFeatQuality(feat_metrics, lasso_classes)
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
    print(caret::confusionMatrix(factor(pred_df$pred_class), factor(pred_df$feat_class),
                                 positive="Good"))
  }
  
  keep_feats <- pred_df %>%
    filter(pred_class=="Good") %>%
    pull(feature)
}

#' Update features in an XCMS object
#' 
#' After metrics have been extracted (typically with `extractChromMetrics`) and 
#' labeling has occurred (typically with `labelFeatsManual` or `labelFeatsLasso`),
#' low quality features can be removed from the XCMS object to improve downstream
#' processing. This function wraps `logModelFeatQuality` and automatically edits
#' the features within the provided XCMS object.
#'
#' @inheritParams logModelFeatQuality
#' @param xcms_obj The XCMS object from which features were initially extracted,
#' usually with `extractChromMetrics`
#'
#' @return The initial xcms_obj but only containing features that exceed the
#' provided quality threshold (as established in likelihood space)
#' @export
#'
#' @examples
#' msnexp_filled <- readRDS(system.file("extdata", "intro_xcms_filled.rds", package="squallms"))
#' peak_data <- makeXcmsObjFlat(msnexp_filled)
#' feat_metrics <- extractChromMetrics(peak_data, verbosity = 0)
#' lasso_classes <- readRDS(system.file("extdata", "intro_lasso_labels.rds", package="squallms"))
#' msnexp_filled <- updateXcmsObjFeats(msnexp_filled, feat_metrics, lasso_classes)
updateXcmsObjFeats <- function(xcms_obj, feature_metrics, feature_labels,
                               log_formula=feat_class~med_cor+med_snr,
                               likelihood_threshold=0.5, verbosity=2){
  good_features <- logModelFeatQuality(feature_metrics, feature_labels, 
                                       likelihood_threshold = likelihood_threshold, 
                                       verbosity = verbosity)
  featureDefinitions(xcms_obj) <- featureDefinitions(xcms_obj)[good_features,]
  
  # xcms_obj@.processHistory[[length(xcms_obj@.processHistory)+1]] <- list(
  #   type="Low-quality feature removal via squallms",
  #   date=date(),
  #   info="Nonstandard processing step",
  #   fileIndex=seq_along(fileNames(xcms_filled)),
  #   `Parameter class`="none",
  #   `MS level(s)`=1
  # )
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
