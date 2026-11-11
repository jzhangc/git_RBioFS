#' @title rbioFS_app
#'
#' @description The web app version of \code{\link{rbioFS}}.
#' @import ggplot2
#' @import shiny
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom randomForest randomForest importance
#' @importFrom parallel detectCores makeCluster stopCluster parApply parLapply
#' @importFrom rpart rpart prune
#' @examples
#' \dontrun{
#' rbioFS_app() # launch the app version by running the function
#' }
#' @export
rbioFS_app <- function(){
  app <- shinyApp(

    ui = fluidPage(
      ## App title
      titlePanel(h1("Function: rbioFS()")),

      ## Sidebar layout with input and output definitions
      sidebarLayout(

        ## Sidebar panel for inputs
        sidebarPanel(
          # adjust the size and scroll
          tags$head(
            tags$style(type = "text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
            tags$style(type = "text/css", "select { max-width: 200px; }"),
            tags$style(type = "text/css", "textarea { max-width: 185px; }"),
            tags$style(type = "text/css", ".jslider { max-width: 200px; }"),
            tags$style(type = "text/css", ".well { max-width: 310px; }"), # size
            tags$style(type = "text/css", ".well { min-width: 310px; }"), # size
            tags$style(type = "text/css", ".span4 { max-width: 310px; }"),
            tags$style(type = "text/css", "form.well { max-height: 95vh; overflow-y: auto; }") # scroll
          ),

          # Input: Select a file
          fileInput("file1", h2("Input CSV File"), # first quotation has the name of the input argument: input$file1. Same as below
                    multiple = TRUE,
                    accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv")),

          # clear screen
          actionButton("clear_ui", "Clear", icon = icon("refresh")),

          # exit
          actionButton("close", "Close App", icon = icon("exclamation"),
                       onclick = "setTimeout(function(){window.close();}, 100);"),

          # Horizontal line
          tags$hr(),

          ## Input block
          h2("Input file settings"),

          # Input: Select separator
          radioButtons("sep",
                       "Separator",
                       choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                       selected = ","), # selected = "," term sets the default value

          # Input: Select number of rows to display
          radioButtons("disp", "Display", choices = c(Head = "head", All = "all"),
                       selected = "head"),

          # column number for group varible
          numericInput(inputId = "targetVar", label = "Column number for group variable",
                       value = 1, step = 1, min = 1),

          # column range for annotatin varibles
          sliderInput(inputId = "annoVar", label = "Column range for annotation variables",
                      min = 1, max = 100, value = c(1, 2)),

          # Horizontal line
          tags$hr(),
          h2("Data imputation and normalization"),
          checkboxInput("impute", "Data imputation", FALSE),
          checkboxInput("quantile", "Quantile normalization", FALSE),
          numericInput(inputId = "impt_anno", label = "Column number for sample variable",
                       value = 1, step = 1, min = 1),   # column number for smaple varible
          radioButtons("imputeMethod",
                       "Data imputation method (if applicable)",
                       choices = c(`Random Forest` = "rf", Mean = "mean", `Random` = "random"),
                       selected = "rf"),
          numericInput(inputId = "imputeIter", label = "Imputation interation (for Random Forest method)",
                       value = 50, step = 5),
          numericInput(inputId = "imputeNtree", label = "Imputation ntree (for Random Forest method)",
                       value = 501, step = 10),

          div(style = "display:inline-block", downloadButton("dlImp", "Save processed data")),

          # Horizontal line
          tags$hr(),

          ## FS
          h2("Overall FS settings"),
          checkboxInput("multicore", "Parallel computing", FALSE),
          numericInput(inputId = "nTree", label = "ntree", value = 1001, step = 100),
          numericInput(inputId = "nTimes", label = "RF iteration",
                       value = 50, step = 10),

          # Horizontal line
          tags$hr(),
          h2("Inital selection"),
          div(style = "display:inline-block", actionButton("run_initial_FS", "Run Initial FS", icon = icon("space-shuttle"))),
          div(style = "display:inline-block", downloadButton("dlInitial", "Save results")),

          h3("Inital selection plot"),
          # column number for group varible
          numericInput(inputId = "initalFS_n", label = "Number of features to plot",
                       value = 1, step = 1, min = 1),
          # Button
          div(style = "display:inline-block",  actionButton("run_initial_FS_plot", "Generate plot", icon = icon("bar-chart"))),
          div(style = "display:inline-block", downloadButton("initialFS_dlPlot", "Save plot")),

          # title
          textInput("initalFS_Title", "Plot title", value = NULL, width = NULL, placeholder = NULL),

          # error bar
          h4("Error bar"),
          radioButtons("initialFS_errorbar", "Type", choices = c(SEM = "sem", SD = "sd"),
                       selected = "sem"),
          numericInput(inputId = "initialFS_errorbarWidth", label = "Error bar width",
                       value = 0.2, step = 0.05),

          # axis
          h4("Axis"),
          textInput("initialFS_xLabel", "x-axis label", value = NULL, width = NULL, placeholder = NULL),
          numericInput(inputId = "initialFS_xTxtSize", label = "x-axis font size",
                       value = 10),
          textInput("initialFS_yLabel", "y-axis label", value = NULL, width = NULL, placeholder = NULL),
          numericInput(inputId = "initialFS_yTxtSize", label = "y-axis font size",
                       value = 10),

          # Space
          tags$br(),

          # Plot: size
          numericInput(inputId = "plotWidth", label = "Plot width",
                       value = 800, step = 10),
          numericInput(inputId = "plotHeight", label = "Plot height",
                       value = 600, step = 10),

          # Horizontal line
          tags$hr(),
          h2("SFS"),
          radioButtons("SFS_mTry",
                       "mtry method",
                       choices = c(Recursion = "recur_default", RF = "rf_default"),
                       selected = "recur_default"),
          div(style = "display:inline-block",  actionButton("run_SFS", "Run SFS", icon = icon("space-shuttle"))),
          div(style = "display:inline-block", downloadButton("dlSFS", "Save results")),

          h3("SFS plot"),

          # Button
          div(style = "display:inline-block", actionButton("run_SFS_plot", "Generate plot", icon = icon("line-chart"))),
          div(style = "display:inline-block", downloadButton("SFS_dlPlot", "Save plot")),

          # title
          textInput("SFS_Title", "Plot title", value = NULL, width = NULL, placeholder = NULL),
          numericInput(inputId = "SFS_symbolSize", label = "Symbol size",
                       value = 2, step = 1),

          # error bar
          h4("Error bar"),
          radioButtons("SFS_errorbar", "Type", choices = c(SEM = "sem", SD = "sd"),
                       selected = "sem"),
          numericInput(inputId = "SFS_errorbarWidth", label = "Error bar width",
                       value = 0.2, step = 0.05),

          # Plot
          h4("Axis"),
          textInput("SFS_xLabel", "x-axis label", value = NULL, width = NULL, placeholder = NULL),
          numericInput(inputId = "SFS_xTxtSize", label = "x-axis font size",
                       value = 10),

          textInput("SFS_yLabel", "y-axis label", value = NULL, width = NULL, placeholder = NULL),
          numericInput(inputId = "SFS_yTxtSize", label = "y-aix font size",
                       value = 10),

          # Space
          tags$br(),
          # Plot: size
          numericInput(inputId = "SFS_plotWidth", label = "Plot width",
                       value = 800, step = 10),
          numericInput(inputId = "SFS_plotHeight", label = "Plot height",
                       value = 600, step = 10)

        ),

        ## Main panel for displaying outputs
        mainPanel(
          # set up tabs
          tabsetPanel(type = "tabs",
                      tabPanel("Raw data", tableOutput("contents")), # "contents" means go to output to find the variable output$contents
                      tabPanel("Processed data (if applicable)", tableOutput("impt")),
                      tabPanel("Initial FS results", verbatimTextOutput("initalFSsum")),
                      tabPanel("Initial FS plot", plotOutput("initalFSplot", height = 480, width = 550)),
                      tabPanel("SFS results", verbatimTextOutput("SFSsum")),
                      tabPanel("SFS plot", plotOutput("SFSplot", height = 480, width = 550)))
        ), fluid = FALSE
      )
    ),

    server = function(input, output, session){
      ## input data check
      # input$file1 will be NULL initially.
      data <- reactive({
        req(input$file1)
        df <- read.table(file = input$file1$datapath, header = TRUE, sep = input$sep,
                         na.strings = "NA", stringsAsFactors = FALSE,
                         check.names = FALSE)
        df[[1]] <- factor(df[[1]], levels = c(unique(df[[1]]))) # avoid R's automatic re-ordering the factors automatically - it will keep the "typed-in" order

        tgt <- factor(as.character(df[, input$targetVar]), levels = unique(df[, input$targetVar]))
        if (input$impute){ # imputation
          if (TRUE %in% apply(df, 2, function(x) any(is.na(x)))){
            impt <- RBioFS::rbioIMP(dfm = df[, -c(input$annoVar[1]:input$annoVar[2])], method = input$imputeMethod,
                                    iter = input$imputeIter, ntree = input$imputeNtree,
                                    fct = tgt, annot = df[, input$impt_anno], transpo = FALSE)
            impt <- data.frame(df[, c(input$annoVar[1]:input$annoVar[2])], impt, check.names = FALSE)
          } else {
            impt <- df
          }
        } else {
          impt <- df
        }

        if (input$quantile){ # normalization
          impt[, -c(input$annoVar[1]:input$annoVar[2])] <- t(rbioNorm(RawData = t(impt[, -c(input$annoVar[1]:input$annoVar[2])]),
                                                                      correctBG = FALSE))
        }

        out <- list(raw = df, imputed = impt)
        return(out)
      })

      ## to display raw data
      # After the user selects and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      output$contents <- renderTable({
        if(input$disp == "head"){
          if (dim(data()$raw)[2] > 5){
            return(head(data()$raw[, 1:5]))
          } else {
            return(head(data()$raw))
          }
        }
        else {
          return(data()$raw)
        }
      }, digits = 5)

      ## to display imputed data
      # After the user selects and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      output$impt <- renderTable({
        if(input$disp == "head"){
          if (dim(data()$imputed)[2] > 5){
            return(head(data()$imputed[, 1:5]))
          } else {
            return(head(data()$imputed))
          }
        }
        else {
          return(data()$imputed)
        }
      }, digits = 5)

      output$dlImp <- downloadHandler(
        filename = function(){paste(substr(noquote(input$file1), 1, nchar(input$file1) - 4),".processed.csv", sep = "")},
        content = function(file){
          write.csv(data()$imputed, file, quote = FALSE, na = "NA", row.names = FALSE)
        }
      )

      ## to set the range
      observe({ # update the selection range for group variable
        updateNumericInput(session = session, inputId = "targetVar", min = 1, max = dim(data()$raw)[2])
      })
      observe({ # update the selection range for group variable
        updateNumericInput(session = session, inputId = "impt_anno", min = 1, max = dim(data()$raw)[2])
      })
      observe({ # update the selection upper range for annotatin variables
        updateSliderInput(session = session, inputId = "annoVar", min = 1, max = dim(data()$raw)[2] - 1,
                          value = c(1, 2))
      })


      ## initial FS
      initialFS_data <- eventReactive(input$run_initial_FS, {
        ## progress bar
        #Create a Progress object
        progress <- Progress$new() # from Shiny package
        progress$set(message = "Computing data", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        # Create a callback function to update progress.
        # Each time this is called:
        # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
        #   distance. If non-NULL, it will set the progress to that value.
        # - It also accepts optional detail text.
        updateProgress <- function(value = NULL, detail = NULL) {
          if (is.null(value)) {
            value <- progress$getValue()
            value <- value + (progress$getMax() - value) / 5
          }
          progress$set(value = value, detail = detail)
        }

        ## pepare the target variable
        tgt <- factor(as.character(data()$imputed[, input$targetVar]), levels = unique(data()$imputed[, input$targetVar]))

        # subset
        x <- data()$imputed[, -c(input$annoVar[1]:input$annoVar[2])] # slider has [1] for min and [2] for max

        ## validate
        # check numbers of features
        validate(need(ncol(x) > 1, "Error: \n
                      Only one feature found in input data. No need to select.\n")) # feature count check

        # check NAs
        validate(need(!TRUE %in% apply(data()$imputed, 2, function(x) any(is.na(x))), "Error: \n
                      NA found in the input data. Try imputation.\n")) # feature count check

        ## FS
        # prepare the data
        training <- as.matrix(x)

        # draw size
        nlvl <- length(levels(tgt))
        size <- min(as.vector(table(tgt))) # down-sampling
        drawSize <- rep(size, nlvl)

        # pre-set an empty matrix with the number of columns same as the number of RF runs
        # note that nrow is the number of features, hence the ncol of the training set
        vimtx <- matrix(nrow = ncol(training), ncol = input$nTimes)
        errmtx <- matrix(nrow = 1, ncol = input$nTimes)

        # RF
        if (!input$multicore){
          tmpFunc <- function(n, m, tmptimes, tmpvimtx, tmperrmtx, tmpTraining, tmpTgt,
                              tmpTree, tmpTry, tmpSize,
                              updateProgress = NULL){ # temp function for recursive RF
            tmploclEnv <- environment() # save the environment local to tmpFunc

            if (n == 0){
              rownames(tmpvimtx) <- colnames(tmpTraining)
              colnames(tmpvimtx) <- c(paste("vi", seq(m - 1), sep = "_"))
              rownames(tmperrmtx) <- "OOB_error_rate"
              colnames(tmperrmtx) <- c(paste("OOB_error_tree", seq(m - 1), sep = "_"))

              tmplst <- list(raw_vi = tmpvimtx, raw_OOB_error = tmperrmtx)
              return(tmplst)

            } else {
              rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, mtry = tmpTry, importance = TRUE,
                                 proximity = TRUE, drawSize = tmpSize)
              impt <- importance(rf, type = 1)
              tmpvimtx[, m] <- impt[, 1] # fill the vi matrix
              tmperrmtx[, m] <- rf$err.rate[tmptimes, 1] # fill the OOB error rate

              # update progress bar
              # If we were passed a progress update function, call it
              if (is.function(updateProgress)){
                text <- paste("Processing RF iteration: ", m, sep = "")
                updateProgress(detail = text)
              }

              tmpFunc(n - 1, m + 1, tmptimes, tmpvimtx, tmperrmtx, tmpTraining, tmpTgt,
                      tmpTree, tmpTry, tmpSize, updateProgress = updateProgress)
            }
          }
          lst <- tmpFunc(n = input$nTimes, m = 1, tmptimes = input$nTree, tmpvimtx = vimtx, tmperrmtx = errmtx, tmpTraining = training, tmpTgt = tgt,
                         tmpTree = input$nTree, tmpTry = max(floor(ncol(x)/3), 2), tmpSize = drawSize,
                         updateProgress = updateProgress)
          phase0mtx_vi <- lst$raw_vi
          phase0mtx_OOB_err <- lst$raw_OOB_error

        } else { # parallel computing
          # recursive RF using par-apply functions
          tmpfunc2 <- function(i){
            rf <- randomForest::randomForest(x = training, y = tgt, ntree = input$nTree, mtry = max(floor(ncol(x)/3), 2), importance = TRUE,
                                             proximity = TRUE, drawSize = drawSize)

            impt <- randomForest::importance(rf, type = 1)
            tmpvimtx <- impt[, 1] # fill the vi matrix
            tmperrmtx <- rf$err.rate[input$nTree, 1] # fill the OOB error rate
            lst <- list(tmpvimtx = tmpvimtx, tmperrmtx = tmperrmtx)
          }

          ## parallel computing
          # set up cpu cluster
          n_cores <- detectCores() - 1
          cl <- makeCluster(n_cores)
          registerDoParallel(cl)

          # foreach parallel
          tmp <- foreach(i = 1:input$nTimes, .export = c("input", "isolate")) %dopar% isolate(tmpfunc2(i))
          vimtx <- foreach(i = 1:input$nTimes, .combine = cbind) %dopar% tmp[[i]]$tmpvimtx
          errmtx <- foreach(i = 1:input$nTimes, .combine = cbind) %dopar% tmp[[i]]$tmperrmtx

          stopCluster(cl) # close connect when done

          rownames(vimtx) <- colnames(training)
          colnames(vimtx) <- c(paste("vi", seq(input$nTimes), sep = "_"))

          rownames(errmtx) <- "OOB_error_rate"
          colnames(errmtx) <- c(paste("OOB_error_tree", seq(input$nTimes), sep = "_"))

          phase0mtx_vi <- vimtx
          phase0mtx_OOB_err <- errmtx
        }

        ## prepare the vi dataframe
        fName_vi <- rownames(phase0mtx_vi)
        fMean_vi <- rowMeans(phase0mtx_vi)
        fSD_vi <- apply(phase0mtx_vi, 1, sd)
        fSEM_vi <- sapply(fSD_vi, function(x)x/sqrt(ncol(phase0mtx_vi)))
        tmpdfm_vi <- data.frame(Target = fName_vi, Mean = fMean_vi, SD = fSD_vi, SEM = fSEM_vi, stringsAsFactors = FALSE)
        tmpdfm_vi <- tmpdfm_vi[order(tmpdfm_vi$Mean), ]
        tmpdfm_vi$Target <- factor(tmpdfm_vi$Target, levels = unique(tmpdfm_vi$Target))

        fMean_OOB_err <- rowMeans(phase0mtx_OOB_err)
        fSD_OOB_err <- apply(phase0mtx_OOB_err, 1, sd)
        fSEM_OOB_err <- fSD_OOB_err/sqrt(ncol(phase0mtx_OOB_err))

        # ranked vi dataframe
        outdfm_vi <- data.frame(tmpdfm_vi[order(tmpdfm_vi$Mean, decreasing = TRUE), ],
                                Rank = c(1:nrow(tmpdfm_vi))) # make sure to resort the dataframe in descenting order.

        # OOB dtaframe
        outdfm_OOB_err <- data.frame(Mean = fMean_OOB_err, SD = fSD_OOB_err, SEM = fSEM_OOB_err, stringsAsFactors = FALSE)
        rownames(outdfm_OOB_err) <- paste(input$nTimes, "trees_OOB_err", sep = "_")

        ## initial feature elimination
        cartTree <- rpart(SD ~ Rank, data = outdfm_vi, cp = 0, minsplit = 2) # CART modelling: classify Rank by SD. Using ANOVA (regression) method.
        mincp <- cartTree$cptable[which(cartTree$cptable[, 4] == min(cartTree$cptable[, 4])) ,1] # extract the minimum cp value
        cartprune <- prune(cartTree, cp = mincp) # prune the tree so that SD values that won't impact Rank classfication are discarded
        minpredv <- min(predict(cartprune)) # obatain the minimum prediciton value (predicted SD) as the SD threshold for Mean

        if (length(which(outdfm_vi$Mean < minpredv)) == 0){ # in the case of VI values don't meet the cut.
          thsd <- ncol(training)
        } else {
          thsd <- min(which(outdfm_vi$Mean < minpredv)) - 1 # compare Mean and SD. Discard all the features with a mean < minimum predicted SD.
        }

        feature_initFS <- as.character(outdfm_vi$Target[1:thsd]) # extract selected features
        training_initFS <- training[, feature_initFS, drop = FALSE] # subsetting the input matrix

        ## return the vi ranking and OOB err dataframes for the initial feature elimination
        outlst <- list(matrix_initial_FS = training_initFS,
                       feature_initial_FS = feature_initFS,
                       recur_vi_summary = outdfm_vi,
                       recur_OOB_err_summary = outdfm_OOB_err)

        return(outlst)
      })

      ## display initial FS resutls
      observeEvent(input$run_initial_FS, {
        output$initalFSsum <- renderPrint(
          initialFS_data()
        )
      })

      # download the initial FS summary
      output$dlInitial <- downloadHandler(
        filename = function(){paste(substr(noquote(input$file1), 1, nchar(input$file1) - 4), ".initialFS.txt", sep = "")},
        content = function(file){
          sink(file, append = FALSE)
          print(initialFS_data())
          sink() # end the dump
        }
      )

      ## initial FS plot
      observe({ # update the selection upper range for annotatin variables
        updateNumericInput(session = session, inputId = "initalFS_n", min = 1, max = nrow(initialFS_data()$recur_vi_summary),
                           value = nrow(initialFS_data()$recur_vi_summary))
      })

      ggplotdata_initialFS <- eventReactive(input$run_initial_FS_plot, {
        loclEnv <- environment()

        pltdfm <- initialFS_data()$recur_vi_summary[order(initialFS_data()$recur_vi_summary$Rank, decreasing = TRUE), ]
        pltdfm <- tail(pltdfm, input$initalFS_n)

        # plotting
        baseplt <- ggplot(pltdfm, aes(x = Target, y = Mean), environment = loclEnv) +
          geom_bar(position="dodge", stat="identity", color="black", fill = "gray66")+
          scale_x_discrete(expand = c(0.01, 0)) +
          scale_y_continuous(expand = c(0.01, 0)) +
          ggtitle(input$initalFS_Title) +
          xlab(input$initialFS_yLabel) + # the arguments for x and y labls are switched as the figure will be rotated
          ylab(input$initialFS_xLabel) + # the arguments for x and y labls are switched as the figure will be rotated
          geom_hline(yintercept = 0) +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
                plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.title = element_text(face = "bold"),
                legend.position = "bottom",
                legend.title = element_blank(),
                axis.text.x = element_text(size = input$initialFS_xTxtSize, angle = 0, hjust = 0.5), # x and y not reversed as they are not associated with the roation of the axes.
                axis.text.y = element_text(size = input$initialFS_yTxtSize, hjust = 0.5)) +
          coord_flip()

        if (input$initialFS_errorbar == "sem"){
          plt <- baseplt +
            geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = input$initialFS_errorbarWidth,
                          position = position_dodge(0.9))
        } else if (input$initialFS_errorbar == "sd") {
          plt <- baseplt +
            geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = input$initialFS_errorbarWidth,
                          position = position_dodge(0.9))
        }

        # extract gtable
        pltgtb <- ggplot_gtable(ggplot_build(plt))

        # add the right side y axis
        Aa <- which(pltgtb$layout$name == "axis-l")
        pltgtb_a <- pltgtb$grobs[[Aa]]
        axs <- pltgtb_a$children[[2]]
        axs$widths <- rev(axs$widths)
        axs$grobs <- rev(axs$grobs)
        axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
        Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
        pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
        pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

        return(pltgtb)
      })

      observeEvent(input$run_initial_FS_plot, {
        output$initalFSplot <- renderPlot({
          grid.draw(ggplotdata_initialFS())
        }, width = input$plotWidth, height = input$plotHeight)
      })

      output$initialFS_dlPlot <- downloadHandler(
        filename = function(){paste(substr(noquote(input$file1), 1, nchar(input$file1) - 4),".initialFS.pdf", sep = "")},
        content = function(file) {
          ggsave(file, plot = ggplotdata_initialFS(),
                 width = (input$plotWidth * 25.4) / 72, height = (input$plotHeight * 25.4) / 72, units = "mm", dpi = 600, device = "pdf")
        }
      )

      ## SFS
      SFS_data <- eventReactive(input$run_SFS, {
        ## progress bar
        #Create a Progress object
        progress <- Progress$new() # from Shiny package
        progress$set(message = "Computing data", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        # Create a callback function to update progress.
        # Each time this is called:
        # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
        #   distance. If non-NULL, it will set the progress to that value.
        # - It also accepts optional detail text.
        updateProgress <- function(value = NULL, detail = NULL) {
          if (is.null(value)) {
            value <- progress$getValue()
            value <- value + (progress$getMax() - value) / 5
          }
          progress$set(value = value, detail = detail)
        }

        ## pepare the target variable
        tgt <- factor(as.character(data()$imputed[, input$targetVar]), levels = unique(data()$imputed[, input$targetVar]))

        # subset
        y <- initialFS_data()$matrix_initial_FS # slider has [1] for min and [2] for max

        # prepare the data
        trainingsfs <- as.matrix(y)

        # SFS
        # prepare draw size. this uses down-sampling if the samples are unbalanced
        nlvl <- length(levels(tgt))
        size <- min(as.vector(table(tgt))) # down-sampling
        drawSize <- rep(size, nlvl)

        # prepare blank tree OOB error matrics
        singleerrmtx <- matrix(nrow = 1, ncol = input$nTimes) # for the recursive OOB error rates from a single tree
        ooberrmtx <- matrix(nrow = ncol(trainingsfs), ncol = input$nTimes) # for the recursive OOB error rates from all trees.

        if (!input$multicore){ ## signle core computing: recursive structure
          tmpFunc <- function(n, m, tmperrmtx, tmpTraining, tmpTgt,
                              tmpTree, tmpTry, tmpSize){
            if (n == 0){
              return(tmperrmtx)

            } else {
              if (tmpTry == "recur_default"){

                if (ncol(tmpTraining) < 4){
                  rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, importance = TRUE,
                                     proximity = TRUE, drawSize = tmpSize)
                } else {
                  rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, mtry = max(floor(ncol(tmpTraining) / 3), 2),
                                     importance = TRUE,
                                     proximity = TRUE, drawSize = tmpSize)
                }

              } else if (tmpTry == "rf_default"){
                rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree,
                                   importance = TRUE,
                                   proximity = TRUE, drawSize = tmpSize)
              }

              tmperrmtx[, m] <- tail(rf$err.rate[, 1], n = 1) # fill the OOB error rate
              tmpFunc(n - 1, m + 1, tmperrmtx, tmpTraining, tmpTgt,
                      tmpTree, tmpTry, tmpSize)
            }
          }

          tmpFunc2 <- function(i, j, tmp2mtx, updateProgress = NULL,
                               ...){
            if (i == 0){
              rownames(tmp2mtx) <- seq(j - 1)
              colnames(tmp2mtx) <- c(paste("OOB_error_tree_rep", seq(input$nTimes), sep = "_"))

              return(tmp2mtx)
            } else {
              tmp2mtx[j, ] <- tmpFunc(n = input$nTimes, m = 1, tmperrmtx = singleerrmtx,
                                      tmpTraining = trainingsfs[, 1:j, drop = FALSE], tmpTgt = tgt, tmpTree = input$nTree, tmpTry = input$SFS_mTry,
                                      tmpSize = drawSize)

              if (is.function(updateProgress)){  # update progress bar
                text <- paste("Processing SFS iteration: ", j, sep = "")
                updateProgress(detail = text)
              }

              tmpFunc2(i - 1, j + 1, tmp2mtx, updateProgress = updateProgress, ...)
            }
          }

          mtxforfunc2 <- ooberrmtx
          ooberrmtx <- tmpFunc2(i = ncol(trainingsfs), j = 1, tmp2mtx = mtxforfunc2, updateProgress = updateProgress) # j is the tree index

        } else { # parallel computing. TBC
          # set up cpu cluster
          n_cores <- detectCores() - 1
          cl <- makeCluster(n_cores)
          registerDoParallel(cl)

          # recursive RF using par-apply functions
          ooberrmtx <- foreach(i = 1:ncol(trainingsfs), .combine = rbind, .export = c("isolate", "input")) %dopar% isolate({
            vct <- vector(length = input$nTimes)
            for (o in 1:input$nTimes){
              rf <- randomForest::randomForest(x = trainingsfs[, 1:i, drop = FALSE], y = tgt, ntree = input$nTree, importance = TRUE,
                                               proximity = TRUE, drawSize = drawSize)
              vct[o] <- tail(rf$err.rate[, 1], n = 1) # compute the OOB error rate
            }
            vct
          })

          rownames(ooberrmtx) <- seq(ncol(trainingsfs))
          colnames(ooberrmtx) <- c(paste("OOB_error_tree_rep", seq(50), sep = "_"))

          stopCluster(cl) # close connect when exiting the function
        }

        ## perpare the summary dataframe for OOB error rates
        ooberrnames <- rownames(ooberrmtx)
        ooberrmean <- rowMeans(ooberrmtx)
        ooberrSD <- apply(ooberrmtx, 1, sd)
        ooberrSEM <- sapply(ooberrSD, function(x)x / sqrt(ncol(ooberrmtx)))
        ooberrsummary <- data.frame(Features = ooberrnames, Mean = ooberrmean, SD = ooberrSD,
                                    SEM = ooberrSEM, stringsAsFactors = FALSE)
        ooberrsummary$Features <- factor(ooberrsummary$Features, levels = unique(ooberrsummary$Features))

        ## output
        mean_min_idx <- which.min(ooberrsummary$Mean)  # index for the minimum mean oob feature group
        sd_min <- ooberrsummary$SD[mean_min_idx]  # oob SD for the feature group above
        minerrsd <- with(ooberrsummary, which(Mean <= (Mean[mean_min_idx] + sd_min)))  # 1sd minimum selection

        minfeatures <- colnames(trainingsfs)[1:min(minerrsd)]
        sfsmatrix <- trainingsfs[, 1:min(minerrsd), drop = FALSE]

        outlst <- list(selected_features = minfeatures,
                       feature_subsets_with_min_OOBerror_plus_1SD = minerrsd,
                       OOB_error_rate_summary = ooberrsummary,
                       SFS_matrix = sfsmatrix)
        return(outlst)
      })

      ## display SFS resutls
      observeEvent(input$run_SFS, {
        output$SFSsum <- renderPrint(
          SFS_data()
        )
      })

      # download the SFS summary
      output$dlSFS <- downloadHandler(
        filename = function(){paste(substr(noquote(input$file1), 1, nchar(input$file1) - 4), ".SFS.txt", sep = "")},
        content = function(file){
          sink(file, append = FALSE)
          print(SFS_data())
          sink() # end the dump
        }
      )

      ## SFS Plot
      ggplotdata_SFS <- eventReactive(input$run_SFS_plot, {
        # validate
        validate(need(nrow(SFS_data()$OOB_error_rate_summary) > 1, "Error: \n
                      Only one feature found in input data. No need to plot.\n")) # feature count check

        # plot
        loclEnv <- environment()

        pltdfm <- SFS_data()$OOB_error_rate_summar

        sfsbaseplt <- ggplot(pltdfm, aes(x = Features, y = Mean, group = 1), environment = loclEnv) +
          geom_line() +
          geom_point(size = input$SFS_symbolSize) +
          scale_x_discrete(expand = c(0.01, 0)) +
          ggtitle(input$SFS_Title) +
          xlab(input$SFS_xLabel) + # the arguments for x and y labls are switched as the figure is rotated
          ylab(input$SFS_yLabel) + # the arguments for x and y labls are switched as the figure is rotated
          geom_vline(xintercept = min(SFS_data()$feature_subsets_with_min_OOBerror_plus_1SD), linetype = "dashed") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
                plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.title = element_text(face = "bold"),
                legend.position = "bottom",
                legend.title = element_blank(),
                axis.text.x = element_text(size = input$SFS_xTxtSize),
                axis.text.y = element_text(size = input$SFS_yTxtSize, hjust = 0.5))

        if (input$SFS_errorbar == "sem"){
          plt <- sfsbaseplt +
            geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = input$SFS_errorbarWidth, position = position_dodge(0.9)) +
            scale_y_continuous(expand = c(0, 0),
                               limits = c(with(pltdfm, min(Mean - SEM) * 0.6),
                                          with(pltdfm, max(Mean + SEM) * 1.2)))
        } else if (input$SFS_errorbar == "sd") {
          plt <- sfsbaseplt +
            geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = input$SFS_errorbarWidth, position = position_dodge(0.9)) +
            scale_y_continuous(expand = c(0, 0),
                               limits = c(with(pltdfm, min(Mean - SD) * 0.6),
                                          with(pltdfm, max(Mean + SD) * 1.2)))
        }

        # extract gtable
        pltgtb <- ggplot_gtable(ggplot_build(plt))

        # add the right side y axis
        Aa <- which(pltgtb$layout$name == "axis-l")
        pltgtb_a <- pltgtb$grobs[[Aa]]
        axs <- pltgtb_a$children[[2]]
        axs$widths <- rev(axs$widths)
        axs$grobs <- rev(axs$grobs)
        axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
        Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
        pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
        pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

        return(pltgtb)
      })

      observeEvent(input$run_SFS_plot, {
        output$SFSplot <- renderPlot({
          grid.draw(ggplotdata_SFS())
        }, width = input$SFS_plotWidth, height = input$SFS_plotHeight)
      })

      output$SFS_dlPlot <- downloadHandler(
        filename = function(){paste(substr(noquote(input$file1), 1, nchar(input$file1) - 4),".SFS.pdf", sep = "")},
        content = function(file) {
          ggsave(file, plot = ggplotdata_SFS(),
                 width = (input$SFS_plotWidth * 25.4) / 72, height = (input$SFS_plotHeight * 25.4) / 72, units = "mm", dpi = 600, device = "pdf")
        }
      )

      # clear button
      observeEvent(input$clear_ui, {
        output$initalFSplot <- renderPlot({})
      })
      observeEvent(input$clear_ui, {
        output$initalFSsum <- renderPrint({cat("")})
      })

      observeEvent(input$clear_ui, {
        output$SFSplot <- renderPlot({})
      })
      observeEvent(input$clear_ui, {
        output$SFSsum <- renderPrint({cat("")})
      })

      observeEvent(input$clear_ui, {
        updateRadioButtons(session = session, inputId = "sep", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                           selected = ",")
      })
      observeEvent(input$clear_ui, {
        updateRadioButtons(session = session, inputId = "disp", choices = c(Head = "head", All = "all"),
                           selected = "head")
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "targetVar", min = 1, value = 1)
      })
      observeEvent(input$clear_ui, {
        updateSliderInput(session = session, inputId = "annoVar", min = 1, max = 100, value = c(1, 2))
      })
      observeEvent(input$clear_ui, {
        updateCheckboxInput(session = session, inputId = "impute", value = FALSE)
      })
      observeEvent(input$clear_ui, {
        updateCheckboxInput(session = session, inputId = "quantile", value = FALSE)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "impt_anno", min = 1, value = 1)
      })
      observeEvent(input$clear_ui, {
        updateRadioButtons(session = session, inputId = "imputeMethod", choices = c(`Random Forest` = "rf", Mean = "mean", `Random` = "random"),
                           selected = "rf")
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "imputeIter", value = 50, step = 5)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "imputeNtree", value = 501, step = 10)
      })
      observeEvent(input$clear_ui, {
        updateCheckboxInput(session = session, inputId = "multicore", value = FALSE)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "nTree", value = 1001, step = 100)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "nTimes", value = 50, step = 10)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "initalFS_n", value = 1, step = 1, min = 1)
      })
      observeEvent(input$clear_ui, {
        updateTextInput(session = session, inputId = "initalFS_Title", value = "", placeholder = NULL)
      })
      observeEvent(input$clear_ui, {
        updateRadioButtons(session = session, inputId = "initialFS_errorbar", choices = c(SEM = "sem", SD = "sd"),
                           selected = "sem")
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "initialFS_errorbarWidth", value = 0.2, step = 0.05)
      })
      observeEvent(input$clear_ui, {
        updateTextInput(session = session, inputId = "initialFS_xLabel", value = "", placeholder = NULL)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "initialFS_xTxtSize", value = 10)
      })
      observeEvent(input$clear_ui, {
        updateTextInput(session = session, inputId = "initialFS_yLabel", value = "", placeholder = NULL)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "initialFS_yTxtSize", value = 10)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "plotWidth", value = 800, step = 10)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "plotWidth", value = 600, step = 10)
      })
      observeEvent(input$clear_ui, {
        updateRadioButtons(session = session, inputId = "SFS_mTry", choices = c(Recursion = "recur_default", RF = "rf_default"),
                           selected = "recur_default")
      })
      observeEvent(input$clear_ui, {
        updateTextInput(session = session, inputId = "SFS_Title", value = "", placeholder = NULL)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "SFS_symbolSize", value = 2, step = 1)
      })
      observeEvent(input$clear_ui, {
        updateRadioButtons(session = session, inputId = "SFS_errorbar", choices = c(SEM = "sem", SD = "sd"),
                           selected = "sem")
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "SFS_errorbarWidth", value = 0.2, step = 0.05)
      })
      observeEvent(input$clear_ui, {
        updateTextInput(session = session, inputId = "SFS_xLabel", value = "", placeholder = NULL)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "SFS_xTxtSize", value = 10)
      })
      observeEvent(input$clear_ui, {
        updateTextInput(session = session, inputId = "SFS_yLabel", value = "", placeholder = NULL)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "SFS_yTxtSize", value = 10)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "SFS_plotWidth", value = 800, step = 10)
      })
      observeEvent(input$clear_ui, {
        updateNumericInput(session = session, inputId = "SFS_plotWidth", value = 600, step = 10)
      })

      # stop and close window
      observe({
        if (input$close > 0) stopApp()  # stop shiny
      })
    }
  )
  runApp(app, launch.browser = TRUE)
}
