labelSingleFeat <- function(row_data, ms1_data) {
    mzbounds <- pmppm(row_data$mzmed, 10)
    rtbounds <- row_data$rtmed + c(-1, 1)
    eic <- ms1_data[mz %between% mzbounds][rt %between% rtbounds] %>%
        arrange(rt)
    if (nrow(eic) < 5) {
        return("Bad")
    }
    dev.new(width = 6, height = 4)
    plot.new()
    plot.window(xlim = rtbounds, ylim = c(0, max(eic$int)), mar = c(1.1, 0.1, 0.1, 0.1))
    for (j in unique(eic$filename)) {
        lines(eic[eic$filename == j, c("rt", "int")])
    }
    abline(v=row_data$rtmed, col="red")
    axis(1)
    title(paste(row_data$feature, round(mzbounds[1], 7)))
    keyinput <- getGraphicsEvent(prompt = "", onKeybd = function(x) {
        return(x)
    })
    dev.off()
    if (length(keyinput)==0){
        return("Quit")
    }
    if (keyinput == "ctrl-[") { # aka Esc button ctrl-[
        return("Quit")
    }
    if( keyinput == "ctrl-H") { # aka Backspace button
        return("Backspace")
    }
    if (!keyinput %in% c("Right", "Left", "Up")) {
        return("Other")
    }
    feat_class <- switch(keyinput,
        "Right" = "Good",
        "Left" = "Bad",
        "Up" = "Revisit"
    )
    feat_class
}
#' Label chromatographic features manually
#'
#' This function provides an interactive interface to a labeling function that
#' plots one chromatographic feature at a time, showing the data from all the
#' files in ms1_data within 1 minute of retention time and 10 ppm m/z space.
#' The plot appears in a new window and is set up to detect key presses which
#' have been bound to specific classifications. Currently, the left arrow key
#' categorizes the feature as "Bad" while the right arrow key categorizes it
#' as "Good". Up and down arrows can be bound to additional classifications if
#' desired by editing the code, though a future update may provide more control.
#' Classifications are returned as a character vector named with the feature
#' IDs even if the full dataset is not classified. Exit the classifier with
#' the Escape key.
#'
#' @param peak_data Flat-form XC-MS data with columns for the bounding box of
#' a chromatographic peak (mzmin, mzmax, rtmin, rtmax) as grouped by a feature
#' ID. Must be provided WITHOUT retention time correction for proper matching
#' to the values in the raw data.
#' @param ms1_data Optional data.table object produced by RaMS containing MS1
#' data with columns for filename, rt, mz, and int. If not provided, the files
#' are detected from the filepath column in peak_data.
#' @param existing_labels A character vector of equal length to the number of
#' features in peak_data named with feature IDs. Can be used to provide a
#' previously-existing partially-labeled dataset or double-check existing
#' classifications.
#' @param selection Either the string "Labeled" or "Unlabeled". If "Unlabeled",
#' the classifier will target entries in the dataset that have not yet received
#' a classification (i.e. those that have NA values in existing_labels). If
#' "Labeled", the classifier will target (and overwrite) existing labels.
#' Otherwise, the classifier will randomly target any feature whether previously
#' classified or not.
#' @param verbosity Scalar value between zero and two determining how much
#' diagnostic information is produced. Zero should return nothing, one should
#' return text-based progress markers, and 2 will return diagnostic plots.
#'
#' @return A character vector named with feature IDs containing the classifications
#' of each peak that was viewed during the interactive phase. NA values indicate
#' those features that were not classified.
#' @export
#'
#' @examples
#' library(xcms)
#' library(dplyr)
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
#' if(interactive()){
#'     manual_labels <- labelFeatsManual(peak_data)
#' }
labelFeatsManual <- function(peak_data, ms1_data = NULL, existing_labels = NULL,
                             selection = "Unlabeled", verbosity = 1) {
    feat_data <- peak_data %>%
        group_by(feature) %>%
        summarise(mzmed = median(mz), rtmed = median(rt))

    if (is.null(existing_labels)) {
        if (selection == "Labeled") {
            message("No existing labels provided, defaulting to unlabeled")
            selection <- "Unlabeled"
        }
    }

    on.exit(return(feat_class_vec))
    plot(1)

    if (is.null(ms1_data)) {
        if (verbosity > 0) {
            message("Grabbing raw MS1 data")
        }
        ms1_data <- grabMSdata(unique(peak_data$filepath), grab_what = "MS1", verbosity = verbosity)$MS1
    }

    feat_class_vec <- rep(NA, nrow(feat_data))

    prev_feat_idx <- numeric()
    backspace_triggered <- FALSE
    while (TRUE) {
        if (selection == "Unlabeled") {
            feature_subset <- which(is.na(feat_class_vec))
        } else if (selection == "Labeled") {
            feature_subset <- which(!is.na(existing_labels))
        } else {
            feature_subset <- seq_along(feat_class_vec)
        }
        if (length(feature_subset) == 0) {
            message("Nothing to label!")
            break
        }
        if (backspace_triggered) {
            chosen_feat_idx <- prev_feat_idx
        } else {
            chosen_feat_idx <- sample(feature_subset, 1)
        }
        if (interactive()) {
            feat_label <- labelSingleFeat(feat_data[chosen_feat_idx, ], ms1_data)
            if (feat_label == "Quit") {
                break
            } else if (feat_label == "Backspace") {
                backspace_triggered <- TRUE
            } else if (feat_label == "Other") {
                init_warn <- getOption("warn")
                options(warn = 1)
                warning("Unrecognized input, skipping")
                options(warn = init_warn)
                backspace_triggered <- FALSE
            } else {
                feat_class_vec[chosen_feat_idx] <- feat_label
                prev_feat_idx <- chosen_feat_idx
                backspace_triggered <- FALSE
            }
        } else {
            warning("Not running in interactive mode, returning NA")
            break
        }
    }
}




plotpeak <- function(feat_ids, interp_df) {
    plot_dt <- interp_df %>%
        filter(feature %in% feat_ids) %>%
        group_by(approx_rt) %>%
        summarise(agg_int_avg = mean(approx_int, na.rm = TRUE), agg_int_iqr = IQR(approx_int, na.rm = TRUE))
    plot_title <- ifelse(length(feat_ids) == 1, feat_ids, "Aggregate")
    par(mar = c(0.1, 0.1, 1.1, 0.1))
    with(plot_dt, plot(approx_rt, agg_int_avg,
        type = "l", lwd = 2, ylim = c(0, 1),
        xlab = "", ylab = "", main = plot_title
    ))
    with(plot_dt, lines(approx_rt, agg_int_avg + agg_int_iqr))
    with(plot_dt, lines(approx_rt, agg_int_avg - agg_int_iqr))
}
classyfeatUI <- function() {
    fluidPage(
        tags$head(tags$script(HTML("Shiny.addCustomMessageHandler('closeWindow', function(m) {window.close();});"))),
        sidebarLayout(
            sidebarPanel(
                h3("Group peakpicker"),
                h4("Settings"),
                numericInput("n_kmeans_groups",
                    label = "Number of k-means groups",
                    value = 4, min = 1, max = 10, step = 1
                ),
                numericInput("n_kmeans_dims",
                    label = "Number of PCs to use for k-means",
                    value = 3, min = 1, max = 10, step = 1
                ),
                actionButton("kmeans_click", label = "Rerun k-means"),
                p(" "),
                plotOutput(outputId = "pcaprops", height = "200px"),
                p(" "),
                actionButton("chosen_good", label = "Flag selection as Good", width = "100%"),
                actionButton("chosen_bad", label = "Flag selection as Bad", width = "100%"),
                actionButton("endsession", label = "Return to R", width = "100%"),
                width = 3
            ),
            mainPanel(
                fluidRow(
                    column(width = 8, plotlyOutput(outputId = "plotlypca")),
                    column(width = 4, plotOutput(outputId = "kmeans_avgpeak"))
                ),
                fluidRow(
                    column(width = 6, plotOutput(outputId = "live_peak", height = "200px")),
                    column(width = 6, plotOutput(outputId = "avg_selected_peak", height = "200px"))
                ),
                width = 9
            )
        )
    )
}
classyfeatServer <- function(input, output, session, pcaoutput, interp_df,
                             feat_vec, verbosity = 0) {
    init_par <- par(no.readonly = TRUE)
    on.exit(par(init_par))

    feat_class_vec <- reactiveVal(setNames(rep(NA, length(feat_vec)), feat_vec))

    output$pcaprops <- renderPlot({
        par(mar = c(2.1, 4.1, 0.1, 0.1))
        perc_exp <- pcaoutput$sdev^2 / sum(pcaoutput$sdev^2)
        barplot(head(perc_exp * 100, 10),
            ylab = "% variance explained", names.arg = paste0("PC", 1:10)
        )
        exp_thresholds <- c(0.2, 0.5, 0.8)
        PCs_to_explain_perc <- vapply(exp_thresholds, function(exp_threshold) {
            which(cumsum(perc_exp) > exp_threshold)[1]
        }, numeric(1))
        legend_text <- paste(paste0(exp_thresholds * 100, "%: "), PCs_to_explain_perc, "PCs")
        legend("topright", legend = legend_text, bty = "n", bg = "transparent")
    })
    kmeaned_df <- reactive({
        input$kmeans_click
        sel_kmeans <- as.data.frame(pcaoutput$rotation[, seq_len(input$n_kmeans_dims)])
        sel_kmeans %>%
            mutate(cluster = factor(kmeans(sel_kmeans, centers = input$n_kmeans_groups)$cluster)) %>%
            rownames_to_column("feature")
    })
    output$plotlypca <- renderPlotly({
        req(input$n_kmeans_groups)
        req(input$n_kmeans_dims)
        gp <- kmeaned_df() %>%
            ggplot(aes(x = PC1, y = PC2, label = feature, color = cluster, key = feature)) +
            geom_text()
        ggplotly(gp, source = "plotlypca") %>% layout(dragmode = "lasso")
    })
    output$kmeans_avgpeak <- renderPlot({
        req(input$n_kmeans_groups)
        req(input$n_kmeans_dims)
        clustergroups <- split(kmeaned_df(), kmeaned_df()$cluster) %>% lapply(`[[`, "feature")
        plot_dt <- interp_df %>%
            left_join(kmeaned_df(), by = "feature") %>%
            group_by(approx_rt, cluster) %>%
            summarise(
                agg_int_avg = mean(approx_int, na.rm = TRUE),
                agg_int_iqr = IQR(approx_int, na.rm = TRUE),
                .groups = "drop"
            )
        gp <- plot_dt %>%
            ggplot(aes(x = approx_rt, color = cluster)) +
            geom_line(aes(y = agg_int_avg), linewidth = 1) +
            geom_line(aes(y = agg_int_avg + agg_int_iqr)) +
            geom_line(aes(y = agg_int_avg - agg_int_iqr)) +
            facet_wrap(~cluster) +
            coord_cartesian(ylim = c(0, 1), clip = "on") +
            theme_bw() +
            theme(legend.position = "none") +
            labs(x = "Normalized retention time", y = "Normalized intensity")
        show(gp)
    })
    output$live_peak <- renderPlot({
        ed_hover <- event_data(source = "plotlypca", event = c("plotly_hover"))
        req(ed_hover)
        plotpeak(ed_hover$key, interp_df)
    })
    output$avg_selected_peak <- renderPlot({
        ed_selected <- event_data(source = "plotlypca", event = c("plotly_selected"))
        req(ed_selected)
        plotpeak(ed_selected$key, interp_df)
    })
    observeEvent(input$chosen_good, {
        good_feats <- event_data(source = "plotlypca", event = "plotly_selected")$key
        init_feat_classes <- feat_class_vec()
        init_feat_classes[good_feats] <- "Good"
        feat_class_vec(init_feat_classes)
        if (verbosity > 0) {
            good_msg <- paste(length(good_feats), "features assigned 'Good' classification")
            message(good_msg)
        }
    })
    observeEvent(input$chosen_bad, {
        bad_feats <- event_data(source = "plotlypca", event = "plotly_selected")$key
        init_feat_classes <- feat_class_vec()
        init_feat_classes[bad_feats] <- "Bad"
        feat_class_vec(init_feat_classes)
        if (verbosity > 0) {
            bad_msg <- paste(length(bad_feats), "features assigned 'Bad' classification")
            message(bad_msg)
        }
    })
    observeEvent(input$endsession, {
        session$sendCustomMessage(type = "closeWindow", message = "message")
        stopApp(feat_class_vec())
    })
    session$onSessionEnded(function() {
        session$sendCustomMessage(type = "closeWindow", message = "message")
        stopApp(isolate(feat_class_vec()))
    })
}
#' Label similar chromatographic features in bulk via interactive selection
#'
#' This function interpolates multi-file chromatograms to a shared set of retention
#' time points then performs a PCA to place similar chromatograms near each other
#' in a reduced dimensionality space. Features can then be labeled in groups
#' instead of one at a time, massively reducing the burden of creating a
#' high-quality dataset. This implementation relies on R's \pkg{shiny} package
#' to provide interactive support in a browser environment and the
#' \pkg{plotly} package for selection tools. Classified features are then
#' returned as a simple R object for downstream use.
#'
#' @inheritParams pickyPCA
#'
#'
#' @return A character vector named with feature IDs containing the classifications
#' of each peak that was viewed during the interactive phase. NA values indicate
#' those features that were not classified.
#' @export
#'
#' @examples
#' library(xcms)
#' library(dplyr)
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
#' if(interactive()){
#'     lasso_labels <- labelFeatsLasso(peak_data)
#' }
labelFeatsLasso <- function(peak_data, ms1_data = NULL, rt_window_width = 1,
                            ppm_window_width = 5, verbosity = 1) {
    if (verbosity > 0) {
        message("Loading libraries")
    }
    if (is.null(ms1_data)) {
        if (verbosity > 0) {
            message("Grabbing raw MS1 data")
        }
        ms1_data <- grabMSdata(unique(peak_data$filepath), grab_what = "MS1", verbosity = verbosity)$MS1
    }

    pickyPCAoutput <- pickyPCA(peak_data, ms1_data, rt_window_width, ppm_window_width)
    interp_df <- pickyPCAoutput$interp_df
    pcaoutput <- prcomp(pickyPCAoutput$pcamat)

    if (verbosity > 0) {
        message("Launching Shiny app")
    }
    shinydef <- shinyApp(
        ui = classyfeatUI,
        server = function(input, output, session) {
            classyfeatServer(input, output, session, pcaoutput, interp_df, feat_vec, verbosity)
        },
        options = c(launch.browser = TRUE)
    )
    feat_vec <- unique(peak_data$feature)
    if (interactive()) {
        feat_class_vec <- runApp(shinydef)
    } else {
        warning("Not running in interactive mode, returning NA")
        feat_class_vec <- rep(NA, length(feat_vec))
        names(feat_class_vec) <- feat_vec
    }

    if (verbosity > 0) {
        message("Returning classification data")
    }
    return(feat_class_vec)
}
