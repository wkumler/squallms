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
#' diagnostic information is produced. 0 should return nothing, 1 should
#' return text-based progress markers, and 2 will return diagnostic plots if
#' available.
#'
#' @return A numeric vector of probabilities returned by the logistic model
#' named by feature ID
#' @export
#'
#' @examples
#' library(xcms)
#' library(MSnbase)
#' mzML_files <- system.file("extdata", package = "RaMS") |>
#'     list.files(full.names = TRUE, pattern = "[A-F].mzML")
#' register(BPPARAM = SerialParam())
#' cwp <- CentWaveParam(snthresh = 0, extendLengthMSW = TRUE, integrate = 2)
#' obp <- ObiwarpParam(binSize = 0.1, response = 1, distFun = "cor_opt")
#' pdp <- PeakDensityParam(
#'     sampleGroups = 1:3, bw = 12, minFraction = 0,
#'     binSize = 0.001, minSamples = 0
#' )
#' xcms_filled <- mzML_files |>
#'     readMSData(msLevel. = 1, mode = "onDisk") |>
#'     findChromPeaks(cwp) |>
#'     adjustRtime(obp) |>
#'     groupChromPeaks(pdp) |>
#'     fillChromPeaks(FillChromPeaksParam(ppm = 5))
#' peak_data <- makeXcmsObjFlat(xcms_filled)
#' feat_metrics <- extractChromMetrics(peak_data, verbosity = 0)
#'
#' # Load demo labels previously assigned using the lasso method
#' lasso_classes <- readRDS(system.file("extdata", "intro_lasso_labels.rds", package = "squallms"))
#' feat_probs <- logModelFeatProb(feat_metrics, lasso_classes)
logModelFeatProb <- function(feature_metrics, feature_labels,
                             log_formula = feat_class ~ med_cor + med_snr,
                             verbosity = 2) {
    # Todo: add checks for formula vars existing in feat_metrics
    model_df <- feature_metrics %>%
        left_join(data.frame(feat_class = feature_labels) %>% rownames_to_column("feature"),
            by = "feature"
        )
    if(all(is.na(model_df$feat_class))){
        stop("All feature labels seem to be NA!")
    }
    glmodel <- model_df %>%
        filter(!is.na(feat_class)) %>%
        mutate(feat_class = ifelse(feat_class == "Bad", 0, 1)) %>%
        glm(formula = log_formula, family = binomial)
    if (verbosity > 0) {
        message("Logistic model regression coefficients:")
        show(glmodel$coefficients)
    }

    if (verbosity > 1) {
        gp <- model_df %>%
            mutate(pred = predict(glmodel, newdata = model_df, type = "response")) %>%
            mutate(feat_class = ifelse(is.na(feat_class), "Unclassified", feat_class)) %>%
            ggplot(aes(x = med_cor, y = med_snr, color = pred, shape = feat_class, label = feature)) +
            geom_point(size = 3) +
            scale_shape_manual(
                breaks = c("Good", "Unclassified", "Bad"),
                values = c(19, 1, 4)
            ) +
            scale_color_gradientn(colors = c("#a41118", "#f4bb23", "#028e34"), limits = c(0, 1)) +
            theme_bw() +
            labs(
                x = "Peak shape coefficient", y = "Signal-to-noise ratio", shape = "Lasso class",
                color = "Predicted\nclass"
            )
        show(gp)
    }

    preds <- predict(glmodel, newdata = model_df, type = "response")
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
#' library(xcms)
#' library(dplyr)
#' library(MSnbase)
#' mzML_files <- system.file("extdata", package = "RaMS") %>%
#'     list.files(full.names = TRUE, pattern = "[A-F].mzML")
#' register(BPPARAM = SerialParam())
#' cwp <- CentWaveParam(snthresh = 0, extendLengthMSW = TRUE, integrate = 2)
#' obp <- ObiwarpParam(binSize = 0.1, response = 1, distFun = "cor_opt")
#' pdp <- PeakDensityParam(
#'     sampleGroups = 1:3, bw = 12, minFraction = 0,
#'     binSize = 0.001, minSamples = 0
#' )
#' xcms_filled <- mzML_files %>%
#'     readMSData(msLevel. = 1, mode = "onDisk") %>%
#'     findChromPeaks(cwp) %>%
#'     adjustRtime(obp) %>%
#'     groupChromPeaks(pdp) %>%
#'     fillChromPeaks(FillChromPeaksParam(ppm = 5))
#' peak_data <- makeXcmsObjFlat(xcms_filled)
#' feat_metrics <- extractChromMetrics(peak_data, verbosity = 0)
#'
#' # Load demo labels previously assigned using the lasso method
#' lasso_classes <- readRDS(system.file("extdata", "intro_lasso_labels.rds", package = "squallms"))
#' feat_classes <- logModelFeatQuality(feat_metrics, lasso_classes)
logModelFeatQuality <- function(feature_metrics, feature_labels,
                                log_formula = feat_class ~ med_cor + med_snr,
                                likelihood_threshold = 0.5,
                                verbosity = 2) {
    feature_preds <- logModelFeatProb(feature_metrics, feature_labels,
        log_formula = log_formula,
        verbosity = verbosity
    )
    pred_df <- feature_metrics %>%
        left_join(data.frame(feat_class = feature_labels) %>% rownames_to_column("feature"),
            by = "feature"
        ) %>%
        left_join(data.frame(pred_prob = feature_preds) %>% rownames_to_column("feature"),
            by = "feature"
        ) %>%
        mutate(pred_class = ifelse(pred_prob > likelihood_threshold, "Good", "Bad"))

    if (verbosity > 0) {
        show(caret::confusionMatrix(factor(pred_df$pred_class), factor(pred_df$feat_class),
            positive = "Good"
        ))
    }

    keep_feats <- pred_df %>%
        filter(pred_class == "Good") %>%
        pull(feature)
}
