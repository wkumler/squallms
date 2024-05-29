
# library(xcms)
# library(dplyr)
# library(MSnbase)
# mzML_files <- system.file("extdata", package = "RaMS") %>%
#   list.files(full.names = TRUE, pattern = "[A-F].mzML")
# register(BPPARAM = SerialParam())
# cwp <- CentWaveParam(snthresh = 0, extendLengthMSW = TRUE, integrate = 2)
# obp <- ObiwarpParam(binSize = 0.1, response = 1, distFun = "cor_opt")
# pdp <- PeakDensityParam(
#   sampleGroups = 1:3, bw = 12, minFraction = 0,
#   binSize = 0.001, minSamples = 0
# )
# xcms_filled <- mzML_files %>%
#   readMSData(msLevel. = 1, mode = "onDisk") %>%
#   findChromPeaks(cwp) %>%
#   adjustRtime(obp) %>%
#   groupChromPeaks(pdp) %>%
#   fillChromPeaks(FillChromPeaksParam(ppm = 5))
# peak_data <- makeXcmsObjFlat(xcms_filled)
# 
# feat_data <- peak_data %>%
#   group_by(feature) %>%
#   summarise(mzmed = median(mz), rtmed = median(rt))
# selection <- "Unlabeled"





library(keys)
library(shiny)

mapped_keys <- c(
  letters, 0:9, "left", "right", "up", "down", "backspace", "esc", "delete"
)

manualFeatUI <- fluidPage(
  tags$head(tags$script(HTML("Shiny.addCustomMessageHandler('closeWindow', function(m) {window.close();});"))),
  useKeys(),
  keysInput("keys", mapped_keys),
  sidebarLayout(
    sidebarPanel(
      h4("Manual feature labeling tool"),
      numericInput("rt_window", label = "Retention time window (minutes)",
                   value = 1, min = 0, step = 1),
      numericInput("ppm", label = "PPM mass error", value = 2.5, min = 0),
      checkboxInput("showint", label = "Show intensity axis?"),
      radioButtons("sequence_method", "What labeling method?",
                   choices = c("Random", "Sequential"), selected = "Sequential"),
      checkboxGroupInput("grab_from", label="Grab a new feature from?",
                         choices = c("Unlabeled", "Good", "Bad"), 
                         selected = "Unlabeled"),
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
)


manualFeatServer <- function(input, output, session, feat_data, input_labeled_vec=NULL) {
  if(!is.null(input_labeled_vec)){
    print("Nonnull input vec")
    if(length(input_labeled_vec)!=nrow(feat_data)){
      stop("input_labeled_vec should have as many entries as feat_data has rows")
    }
    if(!all(names(input_labeled_vec)%in%feat_data$feature)){
      stop("input_labeled_vec should be named using features in feat_data")
    }
    feature_labels <- reactiveVal(
      data.frame(feature=names(input_labeled_vec), label=input_labeled_vec)
    )
    current_feature_id <- reactiveVal(feat_data[which(is.na(input_labeled_vec))[1],]$feature)
    # print(isolate(current_feature_id()))
  } else {
    feature_labels <- reactiveVal(
      data.frame(feature=feat_data$feature, label=NA)
    )
    current_feature_id <- reactiveVal(feat_data[1,]$feature)
  }

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
      
      prev_feature_id(current_feature_id())
      
      print(key_action)
      
      if(key_action%in%c("Good", "Bad", "Revisit", "Ignore")){
        if(length(input$grab_from)==3){
          next_feat_options <- feature_labels()$feature
        } else if(length(input$grab_from)==1){
          if(input$grab_from=="Unlabeled"){
            next_feat_options <- feature_labels()$feature[is.na(feature_labels()$label)]
          } else if(input$grab_from=="Good"){
            next_feat_options <- feature_labels()$feature[feature_labels()$label=="Good"]
          } else if(input$grab_from=="Bad"){
            next_feat_options <- feature_labels()$feature[feature_labels()$label=="Bad"]
          }
        } else if(length(input$grab_from==2)){
          if(input$grab_from==c("Good", "Bad")){
            next_feat_options <- feature_labels()$feature[!is.na(feature_labels()$label)]
          } else {
            next_feat_options <- "Fakefeat"
          }
        } else {
          stop("Unsupported feature selection combination")
        }
        
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
      
      feature_df <- feature_labels()
      feature_df$label[feature_df$feature==current_feature_id()] <- key_action
      feature_labels(feature_df)
    } else {
      if(remap()){
        keymap <- current_keys()
        keymap[remap_idx(),"Key"] <- input$keys
        current_keys(keymap)
        if(remap_idx()<6){
          remap_idx(remap_idx()+1)
        } else {
          remap(FALSE)
          remap_idx(1)
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

manualdef <- shinyApp(
  ui = manualFeatUI,
  server = function(input, output, session) {
    # manualFeatServer(input, output, session, feat_data, input_labeled_vec = feat_class_vec)
    manualFeatServer(input, output, session, feat_data)
  },
  options = c(launch.browser = TRUE)
)

feat_class_vec <- runApp(manualdef)
