
test_that("makeXcmsObjFlat works", {
  expect_contains(names(peak_data), "feature")
  expect_contains(names(peak_data), "mz")
  expect_contains(names(peak_data), "mzmin")
  expect_contains(names(peak_data), "mzmax")
  expect_contains(names(peak_data), "rt")
  expect_contains(names(peak_data), "rtmin")
  expect_contains(names(peak_data), "rtmax")
  expect_contains(names(peak_data), "filename")
  expect_contains(names(peak_data), "filepath")
})

test_that("makeXcmsObjFlat demo specifics", {
  expect_equal(nrow(peak_data), 626)
  expect_equal(length(unique(peak_data$feature)), 195)
  expect_equal(peak_data$feature[626], "FT195")
})

test_that("pickyPCA works", {
  pickyoutput <- pickyPCA(peak_data, msdata$MS1)
  
  expect_length(pickyoutput, 2)
  expect_named(pickyoutput, c("interp_df", "pcamat"))
  expect_s3_class(pickyoutput$interp_df, "data.frame")
  expect_contains(class(pickyoutput$pcamat), "matrix")
  
  # expect_equal(colnames(pickyoutput$pcamat), unique(peak_data$feature))
  expect_contains(unique(peak_data$feature), colnames(pickyoutput$pcamat))
  
  expect_named(pickyoutput$interp_df, c("feature", "filename", "approx_int", "approx_rt"))
  expect_contains(unique(peak_data$feature), unique(pickyoutput$interp_df$feature))
})

test_that("extractChromMetrics works", {
  init_met <- extractChromMetrics(peak_data)
  recalc_met <- extractChromMetrics(peak_data, recalc_betas = TRUE)
  msdata_met <- extractChromMetrics(peak_data, ms1_data = msdata$MS1)
  
  expect_identical(init_met, recalc_met)
  expect_identical(init_met, msdata_met)
  
  expect_identical(colnames(init_met), c("feature", "med_mz", "med_rt",
                                         "med_cor", "med_snr"))
  expect_identical(unique(peak_data$feature), unique(init_met$feature))
  
  expect_lt(max(init_met$med_cor, na.rm = TRUE), 1)
  expect_gt(min(init_met$med_cor, na.rm = TRUE), -1)
  expect_gt(min(init_met$med_snr, na.rm = TRUE), 0)
})

test_that("logModelFeatProb works", {
  # lasso_classes <- labelFeatsLasso(peak_data, msdata$MS1)
  lasso_classes <- readRDS(system.file("extdata", "intro_lasso_labels.rds", 
                                       package="squallms"))
  feat_metrics <- extractChromMetrics(peak_data, verbosity = 0)
  
  expect_identical(length(lasso_classes), nrow(feat_metrics))
  
  lasso_classes <- lasso_classes[!is.na(feat_metrics$med_cor)]
  feat_metrics <- feat_metrics[!is.na(feat_metrics$med_cor),]
  
  model_probs <- logModelFeatProb(feature_metrics = feat_metrics, 
                                  feature_labels = lasso_classes)
  
  expect_lt(max(model_probs, na.rm = TRUE), 1)
  expect_gt(min(model_probs, na.rm = TRUE), 0)
  
  expect_gt(mean(model_probs[lasso_classes=="Good"], na.rm=TRUE), 0.5)
  expect_lt(mean(model_probs[lasso_classes=="Bad"], na.rm=TRUE), 0.5)
})
