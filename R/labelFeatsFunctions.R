
manualFeatUI <- function(){fluidPage(
    tags$head(tags$script(HTML("Shiny.addCustomMessageHandler('closeWindow', function(m) {window.close();});"))),
    useKeys(),
    keysInput("keys", c(
        letters, 0:9, "left", "right", "up", "down", "backspace", "esc", "delete"
    )),
    sidebarLayout(
        sidebarPanel(
            h4("Manual feature labeling tool"),
            numericInput("rt_window", label = "Retention time window (minutes)",
                         value = 1, min = 0, step = 1),
            numericInput("ppm", label = "PPM mass error", value = 2.5, min = 0),
            checkboxInput("showint", label = "Show intensity axis?"),
            radioButtons("sequence_method", "What labeling method?",
                         choices = c("Random", "Sequential"), selected = "Random"),
            textInput("spec_feature", label="Jump to a specific feature?", width="100%"),
            checkboxInput("show_keys", label = "Show keybindings?", value = FALSE),
            uiOutput("keybindings"),
            hr(style = "border-top: 1px solid #000000;"),
            actionButton("endsession", label = "Return to R", width = "100%"),
            style="padding-top: 5px; margin-top: 15px;"
        ),
        mainPanel(
            plotOutput("featureplot")
        )
    )
)}
manualFeatServer <- function(input, output, session, feat_data, ms1_data) {
    feature_labels <- reactiveVal(
        data.frame(feature=feat_data$feature, label=NA)
    )
    current_feature_id <- reactiveVal(feat_data[1,]$feature)
    
    
    prev_feature_id <- reactiveVal()
    
    observeEvent(input$spec_feature, current_feature_id(input$spec_feature), ignoreInit = TRUE)
    
    output$featureplot <- renderPlot({
        req(input$ppm)
        req(input$rt_window)
        req(current_feature_id()%in%feat_data$feature)
        
        feature_data_i <- feat_data[feat_data$feature==current_feature_id(),]
        mzbounds <- pmppm(feature_data_i$mzmed, input$ppm)
        rtbounds <- feature_data_i$rtmed + c(-1, 1) * input$rt_window / 2
        eic <- ms1_data[mz %between% mzbounds][rt %between% rtbounds][order(rt)]
        if(nrow(eic)<5){
            plot(1, type="n")
            text(x=1, y=1, labels="EIC does not contain enough data points")
            title(paste(feature_data_i$feature, round(mzbounds[1], 7)))
        } else {
            plot.new()
            plot.window(xlim = rtbounds, ylim = c(0, max(eic$int)), mar = c(1.1, 0.1, 0.1, 0.1))
            for (j in unique(eic$filename)) {
                lines(eic[eic$filename == j, c("rt", "int")])
            }
            abline(v = feature_data_i$rtmed, col = "red")
            axis(1)
            if(input$showint)axis(2)
            title(paste(feature_data_i$feature, round(mzbounds[1], 7)))
        }
    })
    observeEvent(input$keys, {
        if(!remap()){
            key_action <- current_keys()[current_keys()$Key==input$keys,"Action"]
            
            req(key_action)
            if(key_action%in%c("Good", "Bad", "Revisit", "Ignore")){
                prev_feature_id(current_feature_id())
                feature_df <- feature_labels()
                feature_df$label[feature_df$feature==current_feature_id()] <- key_action
                feature_labels(feature_df)
                next_feat_options <- feature_labels()$feature[is.na(feature_labels()$label)]
                if(length(next_feat_options)==0){
                    session$sendCustomMessage(type = "closeWindow", message = "message")
                    output_vec <- feature_labels()$label
                    names(output_vec) <- feature_labels()$feature
                    message("All features labeled!")
                    stopApp(output_vec)
                }
                
                if(input$sequence_method=="Random"){
                    current_feature_id(sample(next_feat_options, size = 1))
                } else {
                    current_feature_id(next_feat_options[1])
                }
            }
            if(key_action=="Undo once"){
                if(!is.null(prev_feature_id())){
                    current_feature_id(prev_feature_id())
                }
            }
            if(key_action=="Return to R"){
                session$sendCustomMessage(type = "closeWindow", message = "message")
                output_vec <- feature_labels()$label
                names(output_vec) <- feature_labels()$feature
                stopApp(output_vec)
            }
        } else {
            if(remap()){
                keymap <- current_keys()
                if(!input$keys%in%keymap$Key){
                    keymap[remap_idx(),"Key"] <- input$keys
                    current_keys(keymap)
                    if(remap_idx()<6){
                        remap_idx(remap_idx()+1)
                    } else {
                        remap(FALSE)
                        remap_idx(1)
                        updateCheckboxInput(inputId = "show_keys", value = FALSE)
                    }
                }
            }
        }
    })
    
    observeEvent(input$endsession, {
        session$sendCustomMessage(type = "closeWindow", message = "message")
        output_vec <- feature_labels()$label
        names(output_vec) <- feature_labels()$feature
        stopApp(output_vec)
    })
    
    # Remap things
    output$keybindings <- renderUI({
        if(input$show_keys){
            tagList(
                h4("Current keys are:"),
                tableOutput("keymap"),
                uiOutput("remapbutton")
            )
        }
    })
    
    
    current_keys <- reactiveVal(data.frame(
        Action=c("Bad", "Good", "Revisit", "Ignore", "Undo once", "Return to R"),
        Key=c("left", "right", "up", "down", "backspace", "esc")
    ))
    output$keymap <- renderTable(current_keys())
    remap <- reactiveVal(FALSE)
    remap_idx <- reactiveVal(1)
    observeEvent(input$rekey, {
        remap(TRUE)
        current_keys(data.frame(
            Key=c("", "", "", "", "", ""),
            Action=c("Bad", "Good", "Revisit", "Ignore", "Undo once", "Return to R")
        ))
    })
    output$remapbutton <- renderUI({
        if(remap()){
            NULL
        } else {
            actionButton("rekey", "Reassign?", width = "100%")
        }
    })
}

#' Label chromatographic features manually one at a time via interactive interface
#'
#' This function allows the user to view and label individual chromatographic
#' features using the keyboard. Running labelFeatsManual launches a browser
#' which shows a single feature extracted from multiple files. The feature can
#' then be classified by pressing a key (defaults are left=bad, right=good) on 
#' the keyboard which is recorded and a new feature is automatically shown.This 
#' implementation relies on R's \pkg{shiny} package
#' to provide interactive support in a browser environment and the
#' \pkg{keys} package for key binding within a browser. Classified features are then
#' returned as a simple R object for downstream use.
#'
#' @param peak_data Flat-form XC-MS data with columns for the bounding box of
#' a chromatographic peak (mzmin, mzmax, rtmin, rtmax) as grouped by a feature
#' ID. Must be provided WITHOUT retention time correction for proper matching
#' to the values in the raw data.
#' @param ms1_data Optional data.table object produced by RaMS containing MS1
#' data with columns for filename, rt, mz, and int. If not provided, the files
#' are detected from the filepath column in peak_data.
#' @param verbosity Scalar value between zero and two determining how much
#' diagnostic information is produced. 0 should return nothing while 1 will
#' report diagnostic messages.
#' 
#' @return A character vector named with feature IDs containing the classifications
#' of each peak that was viewed during the interactive phase. NA values indicate
#' those features that were not classified.
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
#' if (interactive()) {
#'     manual_labels <- labelFeatsManual(peak_data)
#' }
labelFeatsManual <- function(peak_data, ms1_data = NULL, verbosity = 1) {
    feat_data <- peak_data %>%
        group_by(feature) %>%
        summarise(mzmed = median(mz), rtmed = median(rt))
    
    if (is.null(ms1_data)) {
        if (verbosity > 0) {
            message("Grabbing raw MS1 data")
        }
        ms1_data <- grabMSdata(unique(peak_data$filepath), grab_what = "MS1", 
                               verbosity = verbosity)$MS1
    }
    
    if (verbosity > 0) {
        message("Launching Shiny app")
    }

    manualdef <- shinyApp(
        ui = manualFeatUI,
        server = function(input, output, session) {
            manualFeatServer(input, output, session, feat_data, ms1_data)
        },
        options = c(launch.browser = TRUE)
    )
    if (interactive()) {
        feat_class_vec <- runApp(manualdef)
    } else {
        warning("Not running in interactive mode, returning NA")
        feat_class_vec <- rep(NA, nrow(feat_data))
        names(feat_class_vec) <- feat_class_vec
    }
    
    if (verbosity > 0) {
        message("Returning classification data")
        print(table(feat_class_vec))
    }
    return(feat_class_vec)
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
                    label = "Number of k-means groups (max=10)",
                    value = 4, min = 1, max = 10, step = 1
                ),
                numericInput("n_kmeans_dims",
                    label = "Number of PCs to use for k-means (max=10)",
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
        sel_kmeans <- as.data.frame(pcaoutput$rotation[, seq_len(min(input$n_kmeans_dims, 10))])
        sel_kmeans %>%
            mutate(cluster = factor(kmeans(sel_kmeans, centers = min(input$n_kmeans_groups, 10))$cluster)) %>%
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
#' if (interactive()) {
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
        ms1_data <- grabMSdata(unique(peak_data$filepath), grab_what = "MS1", 
                               verbosity = verbosity)$MS1
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
