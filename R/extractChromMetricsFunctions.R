
makeXcmsObjFlat <- function(xcms_obj, revert_rts=TRUE){
  peak_data_long <- xcms_obj %>%
    chromPeaks() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate(peakidx=row_number()) %>%
    mutate(filename=basename(fileNames(xcms_obj))[sample])
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
  beta_val_df <- peak_data %>%
    select(feature, filename, rtmin, rtmax, mzmin, mzmax) %>%
    left_join(ms1_data, join_args) %>%
    group_by(filename, feature) %>%
    summarise(peak_mz=weighted.mean(mz, int), peak_rt=median(rt), 
              beta_vals=list(qscoreCalculator(rt, int)), .groups = "drop") %>%
    unnest_wider(beta_vals) %>%
    group_by(feature) %>%
    summarise(med_mz=median(peak_mz),
              med_rt=median(peak_rt),
              med_cor=median(beta_cor, na.rm=TRUE), 
              med_snr=median(beta_snr, na.rm=TRUE)
    )
}
qscoreCalculator <- function(rt, int){
  #Check for bogus EICs
  if(length(rt)<5){
    return(list(SNR=NA, peak_cor=NA))
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
pickPCAPixels <- function(peak_data, ms1_data, rt_window_width=1, ppm_window_width=5, verbosity=1){
  join_args <- join_by("filename", between(y$rt, x$rtmin, x$rtmax), between(y$mz, x$mzmin, x$mzmax))
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
  
  pcamat_nounitvar <- pcamat[,!apply(pcamat, 2, var)==0]
  pcafeats <- prcomp(pcamat_nounitvar, center = TRUE, scale. = TRUE)
  if(verbosity>1){
    plot(pcafeats)
    graphics::layout(matrix(c(1,2,3,4), ncol = 2, byrow = TRUE))
    pcafeats$x[,1] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC1")
    pcafeats$x[,2] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC2")
    pcafeats$x[,3] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC3")
    pcafeats$x[,4] %>% matrix(ncol = length(unique(interp_df$filename))) %>% 
      matplot(type="l", col="black", main="PC4")
    graphics::layout(1) # Fix this later
  }
  pcafeats %>%
    pluck("rotation") %>%
    as.data.frame() %>%
    select(1:5) %>%
    rownames_to_column("feature")
}
extractChromMetrics <- function(xcms_obj, recalc_betas=FALSE, verbosity=0, ms1_data=NULL){
  peak_data <- makeXcmsObjFlat(xcms_obj, revert_rts=TRUE)
  
  if(verbosity>0){
    message("Grabbing raw MS1 data")
  }
  if(is.null(ms1_data)){
    ms1_data <- grabMSdata(fileNames(xcms_obj), grab_what = "MS1", verbosity=verbosity)$MS1
  }
  
  if(!"beta_cor"%in%names(peak_data) | recalc_betas){
    if(verbosity>0){
      message("Recalculating beta coefficients")
    }
    beta_df <- calcBetaCoefs(peak_data, ms1_data = msdata$MS1, verbosity=verbosity)
  } else {
    beta_df <- peak_data %>% 
      group_by(feature) %>%
      summarise(med_mz=median(feat_mzmed),
                med_rt=median(feat_rtmed),
                med_cor=median(beta_cor, na.rm=TRUE), 
                med_snr=median(beta_snr, na.rm=TRUE))
  }
  
  if(verbosity>0){
    message("Constructing pixel matrix and performing PCA")
  }
  pca_df <- pickPCAPixels(peak_data, ms1_data = ms1_data, verbosity=verbosity)
  
  full_join(beta_df, pca_df, by="feature")
}
