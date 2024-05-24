#' @md
#' @details
#' squallms is an R package that allows for the identification and removal of
#' low-quality chromatographic features from mass-spectrometry data. This 
#' process involves three steps: first, the calculation of peak quality
#' metrics via [extractChromMetrics()]; second, labeling of high and low quality 
#' features either one at a time with [labelFeatsManual()] or in bulk with
#' [labelFeatsLasso()]; and finally, probabilistic model construction using
#' [logModelFeatQuality()].
#' 
#' The package interfaces neatly with XCMS, providing two additional functions
#' that turn an XCMS object into a flat file ([makeXcmsObjFlat()]) and edit
#' the XCMS object's feature list to contain only high quality features
#' ([updateXcmsObjFeats()]). Data can be provided without using XCMS if it
#' is properly formatted (see the help pages). 
#' 
#' See the package intro on GitHub at https://github.com/wkumler/squallms and
#' explore the vignettes with \code{vignette("intro_to_squallms", package = "squallms")}
#' 
#' @keywords internal
"_PACKAGE"

# Weird NOTE dodges ----
utils::globalVariables("closest") # https://dplyr.tidyverse.org/articles/in-packages.html#join-helpers
utils::globalVariables(c(
    "PC1", "PC2", "adj_rt", "agg_int_avg", "agg_int_iqr",
    "approx_int", "approx_rt", "beta_cor", "beta_snr",
    "cluster", "feat_class", "feat_mzmed", "beta_vals",
    "feat_rtmed", "feature", "filename", "filepath", "int",
    "med_cor", "med_snr", "mtemp", "mzmax", "mzmed", "mzmin",
    "npeaks", "peak_mz", "peak_rt", "peakidx", "pred",
    "pred_class", "pred_prob", "raw_rt", "rtemp", "rtmax",
    "rtmed", "rtmin", "x", "y", ".data", "ppm_window_width",
    "rt_window_width"
))

# Import area ----

# Need to specify a few of these manually to avoid conflicts
#' @rawNamespace import(xcms, except = c(span, groups, collect))
#' @importFrom MSnbase fileNames fromFile readMSData
#' @importFrom RaMS grabMSdata pmppm qplotMS1data
#' @rawNamespace import(dplyr, except = c(between, first, last))
#' @import tidyr
#' @importFrom tibble rownames_to_column
#' @import ggplot2
#' @import shiny
#' @rawNamespace import(plotly, except = c(rename, groups, last_plot, filter))
#' @import data.table
#' @importFrom caret confusionMatrix
#' @rawNamespace import(stats, except = c(lag, filter, smooth, sigma))
#' @rawNamespace import(graphics, except = c(layout))
#' @importFrom grDevices dev.new dev.off getGraphicsEvent
#' @importFrom utils head globalVariables
#' @importFrom MsExperiment readMsExperiment
NULL
