scale_zero_one <- function(x) (x - min(x)) / (max(x) - min(x))

# makeXcmsObjFlat ----
#' Make an XCMS object flat
#'
#' XCMSnExp objects are complicated S4 objects that make
#' it difficult to access information about the features and their associated
#' peaks. This function turns the output into a "flat" file format (a
#' data.frame) so it can interact with tidyverse functions more easily.
#'
#' @param xcms_obj A XCMSnExp object produced by the typical XCMS workflow which
#' should include retention time correction and peak correspondence and filling
#' @param revert_rts Scalar boolean controlling whether the adjusted retention
#' times found in the XCMS object are propagated or returned as-is
#' @param verbosity Scalar numberic. Will report XCMS info if greater than zero
#' or run silently if not (the default).
#'
#' @return A data.frame with columns for feature information (from
#' featureDefinitions: feature, feat_mzmed, feat_rtmed, feat_npeaks, and
#' peakidx) and peak information (from chromPeaks: mz, mzmin, mzmax, rt, rtmin,
#' rtmax, into, intb, maxo, sn, sample) as well as the full path to the
#' associated file and the file name alone (filepath and filename).
#'
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
#' makeXcmsObjFlat(xcms_filled)
makeXcmsObjFlat <- function(xcms_obj, revert_rts = TRUE, verbosity = 0) {
    peak_data_long <- xcms_obj %>%
        chromPeaks() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        mutate(peakidx = row_number()) %>%
        mutate(filepath = fileNames(xcms_obj)[sample]) %>%
        mutate(filename = basename(filepath))
    feat_data <- xcms_obj %>%
        featureDefinitions() %>%
        as.data.frame() %>%
        select(mzmed, rtmed, npeaks, peakidx) %>%
        rownames_to_column("id") %>%
        mutate(rtmed = rtmed / 60)
    peak_data <- feat_data %>%
        unnest_longer(peakidx) %>%
        rename_with(~ paste0("feat_", .x), .cols = -peakidx) %>%
        dplyr::rename(feature = "feat_id") %>%
        left_join(peak_data_long, by = join_by(peakidx))

    if (revert_rts) {
        if (hasAdjustedRtime(xcms_obj)) {
            peak_data <- backToRawRTs(peak_data, xcms_obj, verbosity)
        } else {
            stop("Unable to find adjusted RTs in xcms_obj")
        }
    } else {
        if (verbosity > 0) {
            warn_msg <- paste("revert_rts has been set to FALSE",
                "Please confirm that the xcms_obj RTs match the raw data",
                sep = "\n"
            )
            warning(warn_msg)
        }
    }
    peak_data %>%
        mutate(across(starts_with("rt"), function(x) x / 60))
}

backToRawRTs <- function(peak_data, xcms_obj, verbosity = 0) {
    if (verbosity > 0) {
        dropped_rts <- rtime(dropAdjustedRtime(xcms_obj))
    } else {
        dropped_rts <- rtime(suppressMessages(dropAdjustedRtime(xcms_obj)))
    }
    rt_cors <- data.frame(
        sample = fromFile(xcms_obj),
        adj_rt = rtime(xcms_obj),
        raw_rt = dropped_rts
    )
    # peak_data %>%
    #   select(sample, rtmin) %>%
    #   left_join(rt_cors, by=c("sample", rtmin="adj_rt")) %>%
    #   mutate(rtmin=raw_rt) %>%
    #   select(-raw_rt)
    # Non-equi join below handles inexact matching due to gap filled medians
    peak_data %>%
        left_join(rt_cors, by = join_by("sample", closest(rtmin <= "adj_rt"))) %>%
        mutate(rtmin = raw_rt) %>%
        select(-adj_rt, -raw_rt) %>%
        left_join(rt_cors, by = join_by("sample", closest(rtmax >= "adj_rt"))) %>%
        mutate(rtmax = raw_rt) %>%
        select(-adj_rt, -raw_rt) %>%
        left_join(rt_cors, by = join_by("sample", closest(rt >= "adj_rt"))) %>%
        mutate(rt = raw_rt) %>%
        select(-adj_rt, -raw_rt)

    # Unit test
    # xcms_obj <- readRDS('demodata/falkor/msnexp_filled.rds')
    # msdata <- grabMSdata(fileNames(xcms_obj))
    # msdata$MS1[mz%between%pmppm(118.0865, 10)] %>% qplotMS1data() +
    #   geom_vline(xintercept = unadjusted_rts$rtmin/60, color="green") +
    #   geom_vline(xintercept = unadjusted_rts$rtmax/60, color="red")
    # adjusted_rts <- peak_data %>%
    #   filter(feat_mzmed%between%pmppm(118.0865)) %>%
    #   select(sample, mz, rt, rtmin, rtmax)
    # unadjusted_rts <- peak_data %>%
    #   backToRawRTs(xcms_obj) %>%
    #   filter(feat_mzmed%between%pmppm(118.0865)) %>%
    #   select(sample, mz, rt, rtmin, rtmax)
}

# pickyPCA ----
#' Perform a PCA on multi-file chromatographic data
#'
#' Internal function, mostly.
#'
#' @inheritParams extractChromMetrics
#' @param rt_window_width The width of the retention time window that should
#' be used for PCA construction, in minutes.
#' @param ppm_window_width The width of the m/z window that should be used for
#' PCA construction, in parts per million.
#'
#' @return A list with two named components, interp_df and pcamat. interp_df
#' is the MS1 data associated with each peak interpolated to a retention time
#' spacing shared across the feature. pcamat is the result of cleaning this
#' object up, pivoting it wider, and converting it to a matrix so that a PCA
#' can be performed.
#'
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
#' msdata <- RaMS::grabMSdata(unique(peak_data$filepath), grab_what = "MS1", verbosity = 0)
#' pixel_pca <- pickyPCA(peak_data, msdata$MS1)
pickyPCA <- function(peak_data, ms1_data, rt_window_width = NULL,
                     ppm_window_width = NULL, verbosity = 1) {
    if (is.null(rt_window_width)) {
        rt_window_width <- peak_data %>%
            mutate(rtemp = (rtmax - rtmin) * 2) %>%
            pull(rtemp) %>%
            median(na.rm = TRUE)
    }
    if (is.null(ppm_window_width)) {
        ppm_window_width <- peak_data %>%
            mutate(mtemp = (mzmax - mzmin) * 1e6 / mzmax * 2) %>%
            pull(mtemp) %>%
            median(na.rm = TRUE)
    }

    join_args <- join_by("filename", between(y$rt, x$rtmin, x$rtmax), between(y$mz, x$mzmin, x$mzmax))

    if (verbosity > 0) {
        message("Constructing interpolated data frame")
    }
    interp_df <- peak_data %>%
        select(feature, filename, rt, mz) %>%
        mutate(rtmin = rt - rt_window_width / 2, rtmax = rt + rt_window_width / 2) %>%
        mutate(mzmin = mz - mz * ppm_window_width / 1e6, mzmax = mz + mz * ppm_window_width / 1e6) %>%
        select(-mz, -rt) %>%
        left_join(ms1_data, join_args) %>%
        select(feature, filename, mz, rt, int) %>%
        drop_na() %>%
        group_by(feature, filename) %>%
        filter(n() > 2) %>%
        summarise(approx_int = list(
            approx(rt, int, xout = seq(min(rt), max(rt), length.out = 50), ties = max, rule = 2)$y
        ), .groups = "drop") %>%
        unnest(approx_int) %>%
        group_by(feature, filename) %>%
        mutate(approx_rt = row_number()) %>%
        mutate(approx_int = scale_zero_one(approx_int)) %>%
        ungroup()

    if (verbosity > 0) {
        message("Performing PCA")
    }
    pcamat <- interp_df %>%
        complete(feature, filename, approx_rt) %>%
        group_by(feature, filename) %>%
        mutate(approx_rt = row_number()) %>%
        group_by(feature, approx_rt) %>%
        mutate(approx_int = ifelse(is.na(approx_int), mean(approx_int, na.rm = TRUE), approx_int)) %>%
        ungroup() %>%
        pivot_wider(names_from = feature, values_from = approx_int) %>%
        arrange(filename, approx_rt) %>%
        select(-approx_rt, -filename) %>%
        data.matrix()

    return(list(interp_df = interp_df, pcamat = pcamat))
}

# updateXcmsObjFeats ----
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
#' lasso_classes <- readRDS(system.file("extdata", "intro_lasso_labels.rds", package = "squallms"))
#' xcms_filled <- updateXcmsObjFeats(xcms_filled, feat_metrics, lasso_classes)
updateXcmsObjFeats <- function(xcms_obj, feature_metrics, feature_labels,
                               log_formula = feat_class ~ med_cor + med_snr,
                               likelihood_threshold = 0.5, verbosity = 2) {
    good_features <- logModelFeatQuality(feature_metrics, feature_labels,
        likelihood_threshold = likelihood_threshold,
        verbosity = verbosity
    )
    featureDefinitions(xcms_obj) <- featureDefinitions(xcms_obj)[good_features, ]

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
    #                             fileIndex = seq_along(fileNames(xcms_filled)))
    # xph <- xcms:::XProcessHistory(param = NULL, date. = date(),
    #                              type. = "Unknown",
    #                              fileIndex = seq_along(fileNames(xcms_filled)),
    #                              msLevel = 1)
    # xcms:::addProcessHistory(xcms_filled, ph) %>%
    #   processHistory()
    xcms_obj
}
