library(xcms)
library(dplyr)
library(RaMS)
mzML_files <- system.file("extdata", package = "RaMS") %>%
    list.files(full.names = TRUE, pattern = "[A-F].mzML")
register(BPPARAM = SerialParam())
cwp <- CentWaveParam(snthresh = 0, extendLengthMSW = TRUE, integrate = 2)
obp <- ObiwarpParam(binSize = 0.1, response = 1, distFun = "cor_opt")
pdp <- PeakDensityParam(
    sampleGroups = 1:3, bw = 12, minFraction = 0,
    binSize = 0.001, minSamples = 0
)
xcms_filled <- mzML_files %>%
    readMSData(msLevel. = 1, mode = "onDisk") %>%
    findChromPeaks(cwp) %>%
    adjustRtime(obp) %>%
    groupChromPeaks(pdp) %>%
    fillChromPeaks(FillChromPeaksParam(ppm = 5))

msdata <- grabMSdata(mzML_files, grab_what = "MS1", verbosity = 0)

peak_data <- makeXcmsObjFlat(xcms_filled)




# mzML_files <- "G:/My Drive/AllMeso/mzMLs/pos/" %>%
#   list.files(full.names=TRUE, pattern="190715_(Smp|Poo)")
# mzML_files <- "G:/My Drive/AllMeso/mzMLs/pos/" %>%
#   list.files(full.names=TRUE, pattern="190715_Poo")
# cwp <- CentWaveParam(snthresh = 0, extendLengthMSW = TRUE, integrate = 2,
#                      prefilter = c(5, 1e6))
# pdp <- PeakDensityParam(sampleGroups = seq_along(mzML_files), bw = 12, minFraction = 0,
#                         binSize = 0.001, minSamples = 0)
# obp <- ObiwarpParam(binSize = 0.1, response = 1, distFun = "cor_opt")
# xcms_filled <- mzML_files %>%
#   readMSData(msLevel. = 1, mode = "onDisk") %>%
#   findChromPeaks(cwp) %>%
#   adjustRtime(obp) %>%
#   groupChromPeaks(pdp) %>%
#   fillChromPeaks(FillChromPeaksParam(ppm = 5))
# msdata <- grabMSdata(mzML_files, grab_what = "MS1", verbosity=0)
# peak_data <- makeXcmsObjFlat(xcms_filled)
