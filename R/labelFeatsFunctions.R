
labelSingleFeat <- function(row_data, ms1_data){
  mzbounds <- pmppm(row_data$mzmed, 10)
  rtbounds <- row_data$rtmed+c(-1, 1)
  eic <- ms1_data[mz%between%mzbounds][rt%between%rtbounds] %>% 
    arrange(rt)
  if(nrow(eic)<5)return("Bad")
  dev.new(width=6, height=4)
  plot.new()
  plot.window(xlim=rtbounds, ylim=c(0, max(eic$int)), mar=c(1.1, 0.1, 0.1, 0.1))
  for(j in unique(eic$filename)){
    lines(eic[eic$filename==j, c("rt", "int")])
  }
  axis(1)
  title(paste(row_data$feature, round(mzbounds[1], 7)))
  abline(v=rtbounds, col="red")
  keyinput <- getGraphicsEvent(prompt = "", onKeybd = function(x){
    return(x)
  })
  dev.off()
  if(keyinput=="ctrl-["){ # aka Esc button
    break
  }
  feat_class <- switch(
    keyinput,
    "Right" = "Good",
    "Left" = "Bad",
    "Up" = "Ambiguous",
    "Down" = "Stans only"
  )
  feat_class
}
labelFeatsManual <- function(peak_data, ms1_data, existing_labels=NULL, selection="Unlabeled"){
  feat_data <- peak_data %>%
    group_by(feature) %>%
    summarise(mzmed=unique(feat_mzmed), rtmed=unique(feat_rtmed))
  
  if(is.null(existing_labels)){
    if(selection=="Labeled"){
      message("No existing labels provided, defaulting to unlabeled")
      selection <- "Unlabeled"
    }
    feat_class_vec <- rep(NA, nrow(feat_data))
  } else {
    feat_class_vec <- existing_labels
  }

  on.exit(return(feat_class_vec))
  plot(1)
  
  while(TRUE){
    if(selection=="Unlabeled"){
      feature_subset <- which(is.na(feat_class_vec))
    } else if(selection=="Labeled"){
      feature_subset <- which(!is.na(feat_class_vec))
    } else {
      feature_subset <- seq_along(feat_class_vec)
    }
    if(length(feature_subset)==0){
      message("Nothing to label!")
      break
    }
    chosen_feat_idx <- sample(feature_subset, 1)
    feat_class_vec[chosen_feat_idx] <- labelSingleFeat(feat_data[chosen_feat_idx,], ms1_data)
  }
}




plotpeak <- function(feat_ids, interp_df){
  plot_dt <- interp_df %>%
    filter(feature%in%feat_ids) %>%
    group_by(approx_rt) %>%
    summarise(agg_int_avg=mean(approx_int, na.rm=TRUE), agg_int_iqr=IQR(approx_int, na.rm=TRUE))
  plot_title <- ifelse(length(feat_ids)==1, feat_ids, "Aggregate")
  par(mar=c(0.1, 0.1, 1.1, 0.1))
  with(plot_dt, plot(approx_rt, agg_int_avg, type="l", lwd=2, ylim=c(0, 1), 
                     xlab="", ylab="", main=plot_title))
  with(plot_dt, lines(approx_rt, agg_int_avg+agg_int_iqr))
  with(plot_dt, lines(approx_rt, agg_int_avg-agg_int_iqr))
}
classyfeatUI <- function(){
  fluidPage(
    useShinyjs(),
    extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }", 
                  functions = c("closeWindow")),
    sidebarLayout(
      sidebarPanel(
        h3("Group peakpicker"),
        h4("Settings"),
        numericInput("n_kmeans_groups", label="Number of k-means groups",
                     value = 4, min = 1, max = 10, step = 1),
        numericInput("n_kmeans_dims", label="Number of PCs to use for k-means",
                     value = 3, min = 1, max = 10, step = 1),
        actionButton("kmeans_click", label = "Rerun k-means"),
        p(" "),
        plotOutput(outputId = "pcaprops", height = "200px"), 
        p(" "),
        actionButton("chosen_good", label = "Flag selection as Good"),
        actionButton("chosen_bad", label = "Flag selection as Bad"),
        actionButton("endsession", label = "Return to R"),
        width = 3
      ),
      mainPanel(
        fluidRow(
          column(width=8, plotlyOutput(outputId = "plotlypca")),
          column(width=4, plotOutput(outputId = "kmeans_avgpeak"))
        ),
        fluidRow(
          column(width = 6, plotOutput(outputId = "live_peak", height = "200px")),
          column(width = 6, plotOutput(outputId = "avg_selected_peak", height = "200px"))
        ),
        width=9
      )
    )
  )
}
classyfeatServer <- function(input, output, session, pcaoutput, interp_df, feat_vec, verbosity=0){
  init_par <- par(no.readonly = TRUE)
  on.exit(par(init_par))
  
  feat_class_vec <- reactiveVal(setNames(rep(NA, length(feat_vec)), feat_vec))
  
  output$pcaprops <- renderPlot({
    par(mar=c(2.1, 4.1, 0.1, 0.1))
    perc_exp <- pcaoutput$sdev^2/sum(pcaoutput$sdev^2)
    barplot(head(perc_exp*100, 10), 
            ylab = "% variance explained", names.arg = paste0("PC", 1:10))
    exp_thresholds <- c(0.2, 0.5, 0.8)
    PCs_to_explain_perc <- sapply(exp_thresholds, function(exp_threshold){
      which(cumsum(perc_exp)>exp_threshold)[1]
    })
    legend_text <- paste(paste0(exp_thresholds*100, "%: "), PCs_to_explain_perc, "PCs")
    legend("topright", legend = legend_text, bty='n', bg="transparent")
  })
  kmeaned_df <- reactive({
    input$kmeans_click
    pcaoutput$rotation[,1:input$n_kmeans_dims] %>%
      as.data.frame() %>%
      mutate(cluster=factor(kmeans(., centers=input$n_kmeans_groups)$cluster)) %>%
      rownames_to_column("feature")
  })
  output$plotlypca <- renderPlotly({
    req(input$n_kmeans_groups)
    req(input$n_kmeans_dims)
    gp <- kmeaned_df() %>%
      ggplot(aes(x=PC1, y=PC2, label=feature, color=cluster, key=feature)) +
      geom_text()
    ggplotly(gp, source = "plotlypca") %>% layout(dragmode="lasso")
  })
  output$kmeans_avgpeak <- renderPlot({
    req(input$n_kmeans_groups)
    req(input$n_kmeans_dims)
    clustergroups <- kmeaned_df() %>% 
      split(.$cluster) %>%
      lapply(`[[`, "feature")
    plot_dt <- interp_df %>%
      left_join(kmeaned_df(),by="feature") %>%
      group_by(approx_rt, cluster) %>%
      summarise(agg_int_avg=mean(approx_int, na.rm=TRUE), 
                agg_int_iqr=IQR(approx_int, na.rm=TRUE),
                .groups = "drop")
    gp <- plot_dt %>%
      ggplot(aes(x=approx_rt, color=cluster)) +
      geom_line(aes(y=agg_int_avg), linewidth=1) +
      geom_line(aes(y=agg_int_avg+agg_int_iqr)) +
      geom_line(aes(y=agg_int_avg-agg_int_iqr)) +
      facet_wrap(~cluster) +
      coord_cartesian(ylim=c(0, 1), clip="on") +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x="Normalized retention time", y="Normalized intensity")
    print(gp)
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
    if(verbosity>0){
      good_msg <- paste(length(good_feats), "features assigned 'Good' classification")
      message(good_msg)
    }
  })
  observeEvent(input$chosen_bad, {
    bad_feats <- event_data(source = "plotlypca", event = "plotly_selected")$key
    init_feat_classes <- feat_class_vec()
    init_feat_classes[bad_feats] <- "Bad"
    feat_class_vec(init_feat_classes)
    if(verbosity>0){
      bad_msg <- paste(length(bad_feats), "features assigned 'Bad' classification")
      message(bad_msg)
    }
  })
  observeEvent(input$endsession, {
    js$closeWindow()
    stopApp(feat_class_vec())
  })
  session$onSessionEnded(function() {
    stopApp(feat_class_vec())
  })
}
labelFeatsLasso <- function(peak_data, ms1_data, rt_window_width=1, 
                            ppm_window_width=5, verbosity=1){
  if(verbosity>0){
    message("Loading libraries")
  }
  library(shiny)
  library(shinyjs)
  library(plotly)

  if(verbosity>0){
    message("Constructing interpolated data frame")
  }
  join_args <- join_by("filename", between(y$rt, x$rtmin, x$rtmax), between(y$mz, x$mzmin, x$mzmax))
  feat_vec <- unique(peak_data$feature)
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
  pcaoutput <- interp_df %>%
    complete(feature, filename, approx_rt) %>%
    group_by(feature, filename) %>%
    mutate(approx_rt=row_number()) %>%
    group_by(feature, approx_rt) %>%
    mutate(approx_int=ifelse(is.na(approx_int), mean(approx_int, na.rm=TRUE), approx_int)) %>%
    ungroup() %>%
    pivot_wider(names_from=feature, values_from = approx_int) %>%
    arrange(filename, approx_rt) %>%
    select(-approx_rt, -filename) %>%
    data.matrix() %>%
    prcomp(center = TRUE, scale. = TRUE)
  
  if(verbosity>0){
    message("Launching Shiny app")
  }
  shinydef <- shinyApp(
    ui=classyfeatUI, 
    server = function(input, output, session){
      classyfeatServer(input, output, session, pcaoutput, interp_df, feat_vec, verbosity)
    }, 
    options = c(launch.browser=TRUE))
  feat_class_vec <- runApp(shinydef)
  
  if(verbosity>0){
    message("Saving classification data")
  }
  return(feat_class_vec)
}
