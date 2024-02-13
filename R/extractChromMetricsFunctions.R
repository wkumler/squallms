
extractChromMetrics <- function(xcms_obj=NULL, peakbound_df=NULL, mzML_files=NULL, 
                                pixelPCA=FALSE, recalc_betas=FALSE, verbosity=1){
  if(!is.null(xcms_obj)){
    if(!is.null(peakbound_df)){
      warning(paste0("Both `xcms_obj` and `peakbound_df` were supplied\n",
                     "Using `xcms_obj` and ignoring `peakbound_df`"))
    }
    if(FALSE){ # Fix this later
      stop(paste0("`xcms_obj` does not seem to contain beta_cor or beta_snr\n",
                  "Consider setting recalc_betas=TRUE or rerunning XCMS with ",
                  "verboseBetaColumns=TRUE in the peakpicking step."))
    }
    if(FALSE){ # Fix this later
      stop("`xcms_obj` does not seem to contain groups\n")
    }
    if(pixelPCA & is.null(mzML_files)){
      stop(paste0("pixelPCA requires access to the raw data\n",
                  "Please supply filepaths to the mzML_files argument"))
    }
    flat_msnexp_data <- makeMSnExpFlat(xcms_obj)
    beta_val_df <- flat_msnexp_data %>%
      group_by(feature) %>%
      summarise(med_mz=unique(feat_mzmed),
                med_rt=unique(feat_rtmed),
                med_cor=median(beta_cor, na.rm=TRUE), 
                med_snr=median(beta_snr, na.rm=TRUE)
      )
    if(recalc_betas || pixelPCA){
      message("Extracting raw MS1 data from files")
      ms1_data <- grabMSdata(mzML_files, grab_what = "MS1", verbosity = 1)$MS1
      if(hasAdjustedRtime(xcms_obj)){
        ms1_data <- correctRawRTs(ms1_data, xcms_obj)
      }
    }
    if(recalc_betas){
      peakbound_df <- flat_msnexp_data %>%
        mutate(across(starts_with("rt"), ~.x/60)) %>%
        mutate(filename=basename(mzML_files)[sample])
      beta_val_df <- calcBetaCoefs(peakbound_df, ms1_data)
    }
    if(pixelPCA){
      peak_center_df <- flat_msnexp_data %>% 
        distinct(feature, mzmed=feat_mzmed, rtmed=feat_rtmed)
      message("Extracting pixel matrix from files and performing PCA")
      pca_pixelvals <- pickPCAPixels(peak_center_df, ms1_data, verbosity = verbosity)
      beta_val_df <- beta_val_df %>% left_join(pca_pixelvals, by="feature")
    }
  } else if(!is.null(peakbound_df)){
    if(is.null(mzML_files))stop("If supplying `peakbound_df`, mzML_files must not be supplied")
    if(is.null(peakbound_df$filename))stop("`peakbound_df` must have a `filename` column")
    if(is.null(peakbound_df$feature))stop("`peakbound_df` must have a `feature` column")
    if(is.null(peakbound_df$rtmin))stop("`peakbound_df` must have an `rtmin` column")
    if(is.null(peakbound_df$rtmax))stop("`peakbound_df` must have an `rtmax` column")
    if(is.null(peakbound_df$mzmin))stop("`peakbound_df` must have an `mzmin` column")
    if(is.null(peakbound_df$mzmax))stop("`peakbound_df` must have an `mzmax` column")
    # Check that mzML_files matches unique(peakbound_df$filename) exactly
    message("Extracting raw MS1 data from files")
    ms1_data <- grabMSdata(mzML_files, grab_what = "MS1", verbosity = 1)$MS1
    beta_val_df <- calcBetaCoefs(peakbound_df, ms1_data)
    if(pixelPCA){
      peak_center_df <- peakbound_df %>%
        group_by(feature) %>%
        summarise(mzmed=mean(c(mzmin, mzmax)), rtmed=mean(c(rtmin, rtmax)))
      message("Extracting pixel matrix from files and performing PCA")
      pca_pixelvals <- pickPCAPixels(peak_center_df, ms1_data, verbosity = 1)
      beta_val_df <- beta_val_df %>% left_join(pca_pixelvals, by="feature")
    }
  } else {
    stop("One of `peakbound_df` or `xcms_obj` must not be NULL")
  }
  beta_val_df
}

makeMSnExpFlat <- function(msnexp){
  peak_data_long <- msnexp %>%
    chromPeaks() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate(peakidx=row_number())
  feat_data <- msnexp %>%
    featureDefinitions() %>%
    as.data.frame() %>%
    select(mzmed, rtmed, npeaks, peakidx) %>%
    rownames_to_column("id") %>%
    mutate(rtmed=rtmed/60)
  feat_data %>%
    unnest_longer(peakidx) %>%
    rename_with(~paste0("feat_", .x), .cols = -peakidx) %>%
    dplyr::rename(feature="feat_id") %>%
    left_join(peak_data_long, by = join_by(peakidx))
}


correctRawRTs <- function(ms1_data, msnexp){
  rt_corrections <- data.frame(new_rt=rtime(msnexp)) %>%
    mutate(file_idx=fromFile(msnexp)) %>%
    group_by(file_idx) %>%
    mutate(scan_num=row_number()) %>%
    ungroup() %>%
    mutate(filename=basename(mzML_files)[file_idx]) %>%
    select(filename, scan_num, new_rt)
  ms1_rt_cors <- ms1_data %>%
    distinct(rt, filename) %>%
    group_by(filename) %>%
    mutate(scan_num=row_number()) %>%
    ungroup() %>%
    left_join(rt_corrections, by = join_by(filename, scan_num)) %>%
    mutate(new_rt=new_rt/60)
  ggplot(ms1_rt_cors) + geom_line(aes(x=rt, y=new_rt-rt,group=filename))
  ms1_data <- ms1_data %>%
    left_join(ms1_rt_cors, by = join_by(rt, filename)) %>%
    mutate(rt=new_rt) %>%
    select(-scan_num, -new_rt)
  ms1_data
}

calcBetaCoefs <- function(peakbound_df, ms1_data){
  message("Calculating beta coefficients from EICs")
  within_peaks_eics <- peakboundsToEIC(peakbounds = peakbound_df, ms1_data = ms1_data)
  beta_val_df <- within_peaks_eics %>%
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
peakboundsToEIC <- function(peakbounds, ms1_data, progress=TRUE){
  map2(split(ms1_data, ms1_data$filename), split(peakbounds, peakbounds$filename), 
       function(singlefile_ms1, singlefile_peakbounds){
         pmap(singlefile_peakbounds, getEICfromBounds, ms1_vals=singlefile_ms1)
       }, .progress = progress) %>%
    bind_rows()
}
getEICfromBounds <- function(ms1_vals, mzmin, mzmax, rtmin, rtmax, feature, ...){
  ms1_vals[mz%between%c(mzmin, mzmax)][rt%between%c(rtmin, rtmax)][,feature:=feature]
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

pickPCAPixels <- function(peak_center_df, ms1_data, verbosity=1){
  interp_dt <- pmap(peak_center_df, function(mzmed, rtmed, feature){
    outer_range <- rtmed+c(-0.6, 0.6)
    interp_range <- rtmed+c(-0.5, 0.5)
    interp_points <- seq(interp_range[1], interp_range[2], length.out=50)
    ms1_data[mz%between%pmppm(mzmed, 10)][rt%between%outer_range] %>%
      split(.$filename) %>%
      lapply(function(eic_file){
        if(nrow(eic_file)>2){
          setNames(approx(eic_file$rt, eic_file$int, xout=interp_points, ties = max, rule = 2), c("rt", "int"))
        } else {
          data.frame(rt=numeric(), int=numeric())
        }
      }) %>%
      bind_rows(.id="filename") %>%
      mutate(feature)
  }, .progress = verbosity>0) %>%
    bind_rows() %>%
    group_by(feature, filename) %>%
    mutate(int=int/max(int, na.rm = TRUE))
  pcamat <- interp_dt %>%
    group_by(feature, filename) %>%
    mutate(rt=rank(rt)) %>%
    ungroup() %>%
    complete(feature, filename, rt) %>%
    group_by(feature, rt) %>%
    mutate(int=ifelse(is.na(int), mean(int, na.rm=TRUE), int)) %>%
    ungroup() %>%
    pivot_wider(names_from=feature, values_from = int) %>%
    arrange(filename, rt) %>%
    select(-rt, -filename) %>%
    data.matrix()
  pcamat_nounitvar <- pcamat[,!apply(pcamat, 2, var)==0]
  pcafeats <- prcomp(pcamat_nounitvar, center = TRUE, scale. = TRUE)
  if(verbosity>1){
    plot(pcafeats)
    layout(matrix(c(1,2,3,4), ncol = 2, byrow = TRUE))
    pcafeats$x[,1] %>% matrix(ncol = length(unique(interp_dt$filename))) %>% 
      matplot(type="l", col="black", main="PC1")
    pcafeats$x[,2] %>% matrix(ncol = length(unique(interp_dt$filename))) %>% 
      matplot(type="l", col="black", main="PC2")
    pcafeats$x[,3] %>% matrix(ncol = length(unique(interp_dt$filename))) %>% 
      matplot(type="l", col="black", main="PC3")
    pcafeats$x[,4] %>% matrix(ncol = length(unique(interp_dt$filename))) %>% 
      matplot(type="l", col="black", main="PC4")
    layout(1) # Fix this later
  }
  pcafeats %>%
    pluck("rotation") %>%
    as.data.frame() %>%
    select(1:5) %>%
    rownames_to_column("feature")
}

scale_zero_one <- function(x)(x-min(x))/(max(x)-min(x))
