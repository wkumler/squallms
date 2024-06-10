
library(tidyverse)
library(xcms)
library(MsExperiment)
library(RaMS)
library(squallms)

mzML_files <- system.file("extdata", package = "RaMS") %>%
    list.files(full.names = TRUE, pattern = "[A-F].mzML")
mzML_files <- list.files("mzMLs/", full.names = TRUE)

register(BPPARAM = SerialParam())
cwp <- CentWaveParam(snthresh = 0, extendLengthMSW = TRUE, integrate = 2, prefilter = c(5, 1e6))
obp <- ObiwarpParam(binSize = 0.1, response = 1, distFun = "cor_opt")
pdp <- PeakDensityParam(
    sampleGroups = c(1,1,1,2,2,2), bw = 12, minFraction = 0,
    binSize = 0.001, minSamples = 0
)
fpp <- ChromPeakAreaParam()

msexp <- readMsExperiment(mzML_files, msLevel. = 1, mode = "onDisk")

xcms_filled <- msexp %>%
    findChromPeaks(cwp) %>%
    adjustRtime(obp) %>%
    groupChromPeaks(pdp) %>%
    fillChromPeaks(fpp)

feat_peak_info <- makeXcmsObjFlat(xcms_filled) %>%
    select(feature, starts_with("mz"), starts_with("rt"), filename, filepath)

msdata <- grabMSdata(unique(feat_peak_info$filepath), verbosity = 0)
# shape_metrics <- extractChromMetrics(feat_peak_info, ms1_data = msdata$MS1)


pcaoutput <- feat_peak_info %>%
    filter(feature %in% sprintf("FT%03d", 130:147)) %>%
    pickyPCA(ms1_data = msdata$MS1, rt_window_width = 1, ppm_window_width = 5,
             verbosity = 0)
pcaoutput$interp_df <- pcaoutput$interp_df %>%
    mutate(filename=factor(filename, levels=unique(filename), labels=paste("Sample", 1:6)))

library(cowplot)
singlefeat_chrom_gp <- pcaoutput$interp_df %>%
    filter(feature=="FT137") %>%
    ggplot() +
    geom_line(aes(x=approx_rt, y=approx_int, color=filename), lwd=1) +
    scale_x_continuous(breaks = c(1, 25, 50), labels = c("0", "0.5", "1"), 
                       expand = expansion()) +
    scale_color_manual(breaks=paste("Sample", 1:6), 
                       values=c("#a41118", "#e0670b", "#f4bb23", "#028e34", "#0b505c", "#580770")) +
    labs(y="Intensity normalized\nto maximum", color="Filename") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
singlefeat_heat_gp <- pcaoutput$interp_df %>%
    filter(feature=="FT137") %>%
    ggplot() +
    geom_tile(aes(x = approx_rt, y = filename, fill = approx_int)) +
    scale_x_continuous(breaks = c(1, 25, 50), labels = c("0", "0.5", "1"), 
                       expand = expansion()) +
    scale_y_discrete(expand = expansion()) +
    scale_fill_viridis_c(breaks=c(0, 0.5, 1), begin = 0, end = 0.75) +
    labs(x="Scaled retention time", y="Filename", fill="Intensity\nnormalized\nto maximum") +
    theme_bw()
fig_1a <- plot_grid(singlefeat_chrom_gp, singlefeat_heat_gp, ncol = 1)



fig_1b <- pcaoutput$interp_df %>%
    bind_rows(
        pcaoutput$pcamat %>%
            prcomp() %>%
            .$x %>%
            .[, 1:2] %>%
            as.data.frame() %>%
            mutate(filename=rep(paste("Sample", 1:6), each=50)) %>%
            mutate(approx_rt=rep(1:50, 6)) %>%
            pivot_longer(c("PC1", "PC2"), names_to = c("feature"), values_to = "approx_int") %>%
            mutate(approx_int=(approx_int-min(approx_int))/(max(approx_int)-min(approx_int)))
    ) %>%
    ggplot() +
    geom_tile(aes(x = approx_rt, y = filename, fill = approx_int)) +
    facet_wrap(~feature, nrow = 4) +
    scale_x_continuous(breaks = c(1, 25, 50), labels = c("0", "0.5", "1"), 
                       expand = expansion()) +
    scale_y_discrete(expand = expansion()) +
    scale_fill_viridis_c(breaks=c(0, 0.5, 1), begin = 0, end = 0.75) +
    labs(x = "Scaled retention time", y = "Individual files", 
         fill="Normalized\nintensity") +
    theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(),
          strip.text = element_text(size = 10, margin = margin(t = 2, b = 2)))

library(ggrepel)
set.seed(123)
fig_1c <- pcaoutput$pcamat %>%
    prcomp() %>%
    .$rotation %>%
    .[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    ggplot(aes(x = PC1, y = PC2, label = feature, key = feature)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_point() +
    geom_label_repel(min.segment.length = 0, max.overlaps = Inf, fill="#FFFFFFAA", force = 5) +
    scale_x_continuous(expand = expansion(0.2)) +
    labs(x = "Principal component 1", y = "Principal component 2") +
    theme_bw()

# fig_1 <- plot_grid(
#     plot_grid(singlefeat_chrom_gp, singlefeat_heat_gp, ncol = 1),
#     plot_grid(fig_1b, fig_1c, nrow = 1, rel_widths = c(0.6, 0.4)), 
#     ncol = 1, rel_heights = c(0.4, 0.6)
# )
# ggsave("joss_fig1.png", plot = fig_1, device = "png", width = 6.5, height = 6.5,
#        units = "in", dpi = 300)







fig_1 <- egg::ggarrange(singlefeat_chrom_gp, singlefeat_heat_gp, ncol = 1)
ggsave("joss_fig1.png", plot = fig_1, device = "png", width = 6.5, height = 4,
       units = "in", dpi = 300)
fig_2 <- plot_grid(fig_1b, fig_1c, ncol = 1, rel_heights = c(0.45, 0.35))
ggsave("joss_fig2.png", plot = fig_2, device = "png", width = 6.5, height = 5,
       units = "in", dpi = 300)










labelFeatsLasso(feat_peak_info, ms1_data = msdata$MS1)
