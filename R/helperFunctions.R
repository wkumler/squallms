
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
#'
#' @return A data.frame with columns for feature information (from 
#' featureDefinitions: feature, feat_mzmed, feat_rtmed, feat_npeaks, and
#' peakidx) and peak information (from chromPeaks: mz, mzmin, mzmax, rt, rtmin,
#' rtmax, into, intb, maxo, sn, sample) as well as the full path to the 
#' associated file and the file name alone (filepath and filename).
#' @export
#'
#' @examples
#' msnexp_filled <- readRDS("inst/extdata/msnexp_filled.rds")
#' makeXcmsObjFlat(msnexp_filled)
makeXcmsObjFlat <- function(xcms_obj, revert_rts=TRUE){
  peak_data_long <- xcms_obj %>%
    chromPeaks() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate(peakidx=row_number()) %>%
    mutate(filepath=fileNames(xcms_obj)[sample]) %>%
    mutate(filename=basename(filepath))
  feat_data <- xcms_obj %>%
    featureDefinitions() %>%
    as.data.frame() %>%
    select(mzmed, rtmed, npeaks, peakidx) %>%
    rownames_to_column("id") %>%
    mutate(rtmed=rtmed/60)
  peak_data <- feat_data %>%
    unnest_longer(peakidx) %>%
    rename_with(~paste0("feat_", .x), .cols = -peakidx) %>%
    dplyr::rename(feature="feat_id") %>%
    left_join(peak_data_long, by = join_by(peakidx))
  
  if(revert_rts){
    if(hasAdjustedRtime(xcms_obj)){
      peak_data <- backToRawRTs(peak_data, xcms_obj)
    } else {
      stop("Unable to find adjusted RTs in xcms_obj")
    }
  } else {
    warn_msg <- paste("revert_rts has been set to FALSE", 
                      "Please confirm that the xcms_obj RTs match the raw data",
                      sep = "\n")
    warning(warn_msg)
  }
  peak_data %>%
    mutate(across(starts_with("rt"), function(x)x/60))
}

# pickyPCA ----
#' Perform a PCA on multi-file chromatographic data
#' 
#' Internal function, mostly.
#'
#' @inheritParams extractChromMetrics
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
pickyPCA <- function(peak_data, ms1_data, rt_window_width=NULL, 
                     ppm_window_width=NULL, verbosity=1){
  if(is.null(rt_window_width)){
    rt_window_width <- peak_data %>% 
      mutate(rtemp=(rtmax-rtmin)*2) %>% 
      pull(rtemp) %>% 
      median(na.rm = TRUE)
  }
  if(is.null(ppm_window_width)){
    ppm_window_width <- peak_data %>% 
      mutate(mtemp=(mzmax-mzmin)*1e6/mzmax*2) %>% 
      pull(mtemp) %>% 
      median(na.rm = TRUE)
  }
  
  join_args <- join_by("filename", between(y$rt, x$rtmin, x$rtmax), between(y$mz, x$mzmin, x$mzmax))
  
  if(verbosity>0){
    message("Constructing interpolated data frame")
  }
  interp_df <- peak_data %>%
    select(feature, filename, rt, mz) %>%
    mutate(rtmin=rt-rt_window_width/2, rtmax=rt+rt_window_width/2) %>%
    mutate(mzmin=mz-mz*ppm_window_width/1e6, mzmax=mz+mz*ppm_window_width/1e6) %>%
    select(-mz, -rt) %>%
    left_join(ms1_data, join_args) %>%
    select(feature, filename, mz, rt, int) %>%
    group_by(feature, filename) %>%
    filter(n()>2) %>%
    summarise(approx_int=list(
      approx(rt, int, xout=seq(min(rt), max(rt), length.out=50), ties = max, rule = 2)$y
    ), .groups = "drop") %>%
    unnest(approx_int) %>%
    group_by(feature, filename) %>%
    mutate(approx_rt=row_number()) %>%
    mutate(approx_int=scale_zero_one(approx_int)) %>%
    ungroup()
  
  if(verbosity>0){
    message("Performing PCA")
  }
  pcamat <- interp_df %>%
    complete(feature, filename, approx_rt) %>%
    group_by(feature, filename) %>%
    mutate(approx_rt=row_number()) %>%
    group_by(feature, approx_rt) %>%
    mutate(approx_int=ifelse(is.na(approx_int), mean(approx_int, na.rm=TRUE), approx_int)) %>%
    ungroup() %>%
    pivot_wider(names_from=feature, values_from = approx_int) %>%
    arrange(filename, approx_rt) %>%
    select(-approx_rt, -filename) %>%
    data.matrix()
  
  return(list(interp_df=interp_df, pcamat=pcamat))
}

# Import area ----

#' @import xcms
#' @import RaMS
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import shiny
#' @importFrom shinyjs useShinyjs, extendShinyjs
#' @importFrom stats approx
NULL
