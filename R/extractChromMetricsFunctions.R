

backToRawRTs <- function(peak_data, xcms_obj){
  rt_cors <- data.frame(
    sample=fromFile(xcms_obj),
    adj_rt=rtime(xcms_obj),
    raw_rt=rtime(suppressMessages(dropAdjustedRtime(xcms_obj)))
  )
  # peak_data %>%
  #   select(sample, rtmin) %>%
  #   left_join(rt_cors, by=c("sample", rtmin="adj_rt")) %>%
  #   mutate(rtmin=raw_rt) %>%
  #   select(-raw_rt)
  # Non-equi join below handles inexact matching due to gap filled medians
  peak_data %>%
    left_join(rt_cors, by=join_by("sample", closest(rtmin<="adj_rt"))) %>%
    mutate(rtmin=raw_rt) %>%
    select(-adj_rt, -raw_rt) %>%
    left_join(rt_cors, by=join_by("sample", closest(rtmax>="adj_rt"))) %>%
    mutate(rtmax=raw_rt) %>%
    select(-adj_rt, -raw_rt) %>%
    left_join(rt_cors, by=join_by("sample", closest(rt>="adj_rt"))) %>%
    mutate(rt=raw_rt) %>%
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
calcBetaCoefs <- function(peak_data, ms1_data, verbosity=1){
  join_args <- join_by("filename", between(y$rt, x$rtmin, x$rtmax), between(y$mz, x$mzmin, x$mzmax))
  joined_df <- peak_data %>%
    select(feature, filename, rtmin, rtmax, mzmin, mzmax) %>%
    left_join(ms1_data, join_args) %>%
    distinct()
  beta_val_df <- joined_df %>%
    group_by(feature, filename) %>%
    summarise(peak_mz=weighted.mean(mz, int), peak_rt=median(rt), 
              beta_vals=list(qscoreCalculator(rt, int)), .groups = "drop_last") %>%
    unnest_wider(beta_vals, simplify = FALSE) %>%
    summarise(med_mz=median(peak_mz, na.rm=TRUE),
              med_rt=median(peak_rt, na.rm=TRUE),
              med_cor=median(beta_cor, na.rm=TRUE), 
              med_snr=median(beta_snr, na.rm=TRUE))
}
qscoreCalculator <- function(rt, int, na.rm=TRUE){
  rt <- rt[!is.na(rt)]
  int <- int[!is.na(int)]
  #Check for bogus EICs
  if(length(rt)<5){
    return(c(beta_snr=NA, beta_cor=NA))
  }
  #Calculate where each rt would fall on a beta dist (accounts for missed scans)
  scaled_rts <- (rt-min(rt))/(max(rt)-min(rt))
  
  # Create a couple different skews and test fit
  maybe_skews <- c(2.5,3,4,5) #Add 7 to catch more multipeaks and more noise
  #Add 2 to catch very slopey peaks and more noise
  best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
  peak_cor <- cor(perf_peak, int)
  
  
  #Calculate the normalized residuals
  residuals <- int/max(int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(int)-min(int))/sd(norm_residuals*max(int))
  #Return the quality score
  return(c(beta_snr=SNR, beta_cor=peak_cor))
}
scale_zero_one <- function(x)(x-min(x))/(max(x)-min(x))
pickPCAPixels <- function(peak_data, ms1_data, verbosity=1){
  pickyPCAoutput <- pickyPCA(peak_data, ms1_data, rt_window_width, ppm_window_width)
  interp_df <- pickyPCAoutput$interp_df
  pcamat <- pickyPCAoutput$pcamat
  pcamat_nounitvar <- pcamat[,!apply(pcamat, 2, var)==0]
  pcafeats <- prcomp(pcamat_nounitvar, center = TRUE, scale. = TRUE)
  if(verbosity>1){
    plot(pcafeats)
    pcafeats$x[,1] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC1")
    pcafeats$x[,2] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC2")
    pcafeats$x[,3] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC3")
    pcafeats$x[,4] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC4")
  }
  pcafeats %>%
    .data$rotation %>%
    as.data.frame() %>%
    select(1:5) %>%
    rownames_to_column("feature")
}
#' Extract metrics of chromatographic peak quality
#' 
#' This function takes flat-form XC-MS data (i.e. the format produced by
#' `makeXcmsObjectFlat`) and calculates metrics of peak shape and
#' similarity for downstream processing (e.g. by \code{\link{updateXcmsObjFeats}}).
#' The core metrics are those described in \url{https://doi.org/10.1186/s12859-023-05533-4},
#' corresponding to peak shape (correlation coefficient between the raw data
#' and an idealized bell curve, beta_cor) and within-peak signal-to-noise (
#' maximum intensity divided by the standard deviation of the residuals after
#' the best-fitting bell is subtracted from the raw data). Additionally,
#' the function interpolates the raw data to fixed retention time intervals
#' and performs a principal components analysis (PCA) across a file-by-RT matrix
#' which typically extracts good-looking peaks in the first component (PC1).
#' These functions are calculated on the raw MS data as obtained via 
#' \pkg{RaMS}, so either the filepaths must be included in the peak_data object
#' or supplied as ms1_data.
#'
#' @param peak_data Flat-form XC-MS data with columns for the bounding box of
#' a chromatographic peak (mzmin, mzmax, rtmin, rtmax) as grouped by a feature 
#' ID. Must be provided WITHOUT retention time correction for proper matching
#' to the values in the raw data.
#' @param recalc_betas Scalar boolean controlling whether existing beta values
#' (as calculated in the XCMS object when peakpicked with \code{CentWaveParam(verboseBetaColumns=TRUE))}
#' should be used as-is (the default) or recalculated. See 
#' \url{https://github.com/sneumann/xcms/pull/685} for differences in these 
#' implementations.
#' @param ms1_data Optional data.table object produced by RaMS containing MS1
#' data with columns for filename, rt, mz, and int. If not provided, the files
#' are detected from the filepath column in peak_data.
#' @param verbosity Scalar value between zero and two determining how much 
#' diagnostic information is produced. Zero should return nothing, one should
#' return text-based progress markers, and 2 will return diagnostic plots.
#'
#' @return A data.frame containing one row for each feature in peak_data, with
#' columns containing the median peak shape value (med_cor), and the median SNR
#' value (med_snr).
#' @export
#'
#' @examples
#' \dontrun{
#' library(xcms)
#' library(dplyr)
#' mzML_files <- system.file("extdata", package = "RaMS") %>%
#'   list.files(full.names=TRUE, pattern="[A-F].mzML")
#' register(BPPARAM = SerialParam())
#' cwp <- CentWaveParam(snthresh = 0, extendLengthMSW = TRUE, integrate = 2)
#' obp <- ObiwarpParam(binSize = 0.1, response = 1, distFun = "cor_opt")
#' pdp <- PeakDensityParam(sampleGroups = 1:3, bw = 12, minFraction = 0, 
#'                         binSize = 0.001, minSamples = 0)
#' xcms_filled <- mzML_files %>%
#'   readMSData(msLevel. = 1, mode = "onDisk") %>%
#'   findChromPeaks(cwp) %>%
#'   adjustRtime(obp) %>%
#'   groupChromPeaks(pdp) %>%
#'   fillChromPeaks(FillChromPeaksParam(ppm = 5))
#' peak_data <- makeXcmsObjFlat(msnexp_filled)
#' feat_metrics <- extractChromMetrics(peak_data, recalc_betas = TRUE)
#' }
extractChromMetrics <- function(peak_data, recalc_betas=FALSE, ms1_data=NULL, 
                                verbosity=0){
  if(!is(peak_data, "data.frame")){
    stop("peak_data must be a data.frame object")
  }
  necessary_colnames <- c("feature", "mz", "mzmin", "mzmax", "rt", "rtmin", 
                          "rtmax", "filepath", "filename")
  missing_colnames <- setdiff(necessary_colnames, colnames(peak_data))
  if(length(missing_colnames)!=0){
    stop(paste("peak_data is missing necessary column(s):", 
               paste(missing_colnames, collapse = ", ")))
  }
  
  if(is.null(ms1_data)){
    if(verbosity>0){
      message("Grabbing raw MS1 data")
    }
    ms1_data <- grabMSdata(unique(peak_data$filepath), grab_what = "MS1", verbosity=verbosity)$MS1
  }
  
  if(!"beta_cor"%in%names(peak_data) | recalc_betas){
    if(verbosity>0){
      message("Recalculating beta coefficients")
    }
    beta_df <- calcBetaCoefs(peak_data, ms1_data = ms1_data, verbosity=verbosity)
  } else {
    beta_df <- peak_data %>% 
      group_by(feature) %>%
      summarise(med_mz=median(feat_mzmed),
                med_rt=median(feat_rtmed),
                med_cor=median(beta_cor, na.rm=TRUE), 
                med_snr=median(beta_snr, na.rm=TRUE))
  }
  
  # if(verbosity>0){
  #   message("Constructing pixel matrix and performing PCA")
  # }
  # pca_df <- pickPCAPixels(peak_data, ms1_data = ms1_data, 
  #                         rt_window_width = rt_window_width, 
  #                         ppm_window_width = ppm_window_width,
  #                         verbosity=verbosity)
  
  # full_join(beta_df, pca_df, by="feature")
  beta_df
}
