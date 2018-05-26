#' @title rbioFS_PCA_app
#'
#' @description The web app version of \code{\link{rbioFS_PCA}}.
#' @import ggplot2
#' @import shiny
#' @examples
#' \dontrun{
#' rbioFS_PCA_app() # launch the app version by running the function
#' }
#' @export
rbioFS_PCA_app <- function(){
  app <- shinyApp(
    ui = fluidPage(
      ## App title
      titlePanel(h1("Function: rbioFS_PCA()")),

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
          numericInput(inputId = "groupIDVar", label = "Column number for group variable",
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

          ## PCA
          h2("Overall PCA settings"),
          checkboxInput("scaleData", "Scale data", TRUE),
          div(style = "display:inline-block", actionButton("run_PCA", "Run and plot PCA", icon = icon("space-shuttle"))),

          # Horizontal line
          tags$hr(),
          numericInput(inputId = "xTickLblSize", label = "X-axis tick label size",
                       value = 10, step = 0.5, min = 1),
          numericInput(inputId = "yTickLblSize", label = "Y-axis tick label size",
                       value = 10, step = 0.5, min = 1),
          textInput("fontType", "Font type", value = "sans", width = NULL, placeholder = NULL),
          actionButton(inputId = "fontTable", "Font table", icon = icon("th"), onclick = "window.open('http://kenstoreylab.com/wp-content/uploads/2015/08/R-font-table.png', '_blank')"),

          ## PC boxplot
          h3("PC boxplot"),

          # Button
          div(style = "display:inline-block", downloadButton("boxplot_dlPlot", "Save plot")),

          # other settings
          h4("Detailed settings"),
          textInput("boxplot.Title", "Plot title", value = NULL, width = NULL, placeholder = NULL),
          numericInput(inputId = "boxplot.Width", label = "Plot width", value = 800, step = 10),
          numericInput(inputId = "boxplot.Height", label = "Plot height", value = 600, step = 10),

          ## biplot
          h3("Biplot"),
          checkboxInput("loadingplot", "Superimpose loading plot", TRUE),
          numericInput(inputId = "PC_1st", label = "First PC to plot", value = 1, step = 1, min = 1),
          numericInput(inputId = "PC_2nd", label = "Second PC to plot", value = 2, step = 1, min = 1),
          checkboxInput("ellipse", "ellipase", FALSE),
          numericInput(inputId = "ellipse_conf", label = "ellipse confidence",
                       value = 0.6, step = 0.05, max = 1, min = 0),

          # Button
          div(style = "display:inline-block", downloadButton("biplot_dlPlot", "Save plot")),

          # other settings
          h4("Detailed settings"),
          textInput("biplot.Title", "Plot title", value = NULL, width = NULL, placeholder = NULL),
          numericInput(inputId = "loadingSize", label = "loading plot label size",
                       value = 3, step = 1, max = 1, min = 0),
          numericInput(inputId = "biplot.SymbolSize", label = "biplot symbol size",
                       value = 2, step = 0.5, min = 1),

          # Space
          tags$br(),
          # Plot: size
          numericInput(inputId = "biplot.Width", label = "Plot width",
                       value = 800, step = 10),
          numericInput(inputId = "biplot.Height", label = "Plot height",
                       value = 600, step = 10)

        ),

        ## Main panel for displaying outputs
        mainPanel(
          # set up tabs
          tabsetPanel(type = "tabs",
                      tabPanel("Raw data", tableOutput("contents")), # "contents" means go to output to find the variable output$contents
                      tabPanel("Processed data (if applicable)", tableOutput("impt")),
                      tabPanel("Boxplot", plotOutput("PCAboxplot", height = 480, width = 550)),
                      tabPanel("Biplot", plotOutput("PCAbiplot", height = 480, width = 550)))
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

        tgt <- factor(as.character(df[, input$groupIDVar]), levels = unique(df[, input$groupIDVar]))
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
        updateNumericInput(session = session, inputId = "groupIDVar", min = 1, max = dim(data()$raw)[2])
      })
      observe({ # update the selection range for group variable
        updateNumericInput(session = session, inputId = "impt_anno", min = 1, max = dim(data()$raw)[2])
      })
      observe({ # update the selection upper range for annotatin variables
        updateSliderInput(session = session, inputId = "annoVar", min = 1, max = dim(data()$raw)[2] - 1,
                          value = c(1, 2))
      })
      observe({ # update the selection upper range for annotatin variables
        updateNumericInput(session = session, inputId = "PC_1st",
                           min = 1, max = dim(data()$imputed[, -c(input$annoVar[1]:input$annoVar[2])])[2],
                           value = 1)
      })
      observe({ # update the selection upper range for annotatin variables
        updateNumericInput(session = session, inputId = "PC_2nd",
                           min = 1, max = dim(data()$imputed[, -c(input$annoVar[1]:input$annoVar[2])])[2],
                           value = 2)
      })

      ## PCA data
      PCA_data <- eventReactive(input$run_PCA, {
        # subset and setup PCA
        x <- data()$imputed[, -c(input$annoVar[1]:input$annoVar[2])] # slider has [1] for min and [2] for max
        PCA <- prcomp(x, scale. = input$scaleData)

        return(PCA)
      })

      ## boxplot
      ggplotdata_boxplot <- eventReactive(input$run_PCA, {
        varpp <- 100 * summary(PCA_data())$importance[2, ] # extract and calcualte the proportion of variance
        boxdfm <- data.frame(PC = as.numeric(gsub("PC", "", names(varpp))), varpp = varpp)

        boxplt <- ggplot(data = boxdfm, aes(x = PC, y = varpp, group = 1)) +
          geom_bar(position = "dodge", stat = "identity", color = "black", fill = "gray66") +
          scale_x_continuous() +
          scale_y_continuous(expand = c(0, 0),
                             limits = c(0, with(boxdfm, ceiling(max(varpp))) *1.1)) +
          ggtitle(input$boxplot.Title) +
          xlab("PC") +
          ylab("Proportion of variance (%)") +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
                plot.title = element_text(face = "bold", family = input$fontType, hjust = 0.5),
                axis.title = element_text(face = "bold", family = input$fontType),
                legend.position = "bottom",legend.title = element_blank(),legend.key = element_blank(),
                axis.text.x = element_text(size = input$xTickLblSize, family = input$fontType),
                axis.text.y = element_text(size = input$yTickLblSize, family = input$fontType, hjust = 0.5))

        return(boxplt)
      })

      observe({
        output$PCAboxplot <- renderPlot({
          grid.draw(ggplotdata_boxplot())
        }, width = input$boxplot.Width, height = input$boxplot.Height)
      })

      output$boxplot_dlPlot <- downloadHandler(
        filename = function(){paste(substr(noquote(input$file1), 1, nchar(input$file1) - 4),".boxplot.pdf", sep = "")},
        content = function(file) {
          ggsave(file, plot = ggplotdata_boxplot(),
                 width = (input$boxplot.Width * 25.4) / 72, height = (input$boxplot.Height * 25.4) / 72,
                 units = "mm", dpi = 600, device = "pdf")
        }
      )

      ## biplot
      ggplotdata_biplot <- eventReactive(input$run_PCA, {
        varpp <- 100 * summary(PCA_data())$importance[2, ] # extract and calcualte the proportion of variance

        # prepare for scatter plot values (i.e. sample score)
        sampleScore <- data.frame(PCA_data()$x[, c(paste0("PC", as.character(input$PC_1st)), paste0("PC", as.character(input$PC_2nd)))],
                                  check.names = FALSE) # extract rotated sample scores
        sampleScore$Group <- factor(as.character(data()$imputed[, input$groupIDVar]), levels = unique(data()$imputed[, input$groupIDVar]))
        names(sampleScore)[1:2] <- c("axis1", "axis2")

        # prepare for loading plot values (i.e. loading value for variables)
        loadingValue <- data.frame(PCA_data()$rotation[, c(paste0("PC", as.character(input$PC_1st)), paste0("PC", as.character(input$PC_2nd)))],
                                   check.names = FALSE) # extract loading/rotation/eigenvectors for variables
        names(loadingValue) <- c("axis1", "axis2") # give a generic variable name for the ratation dataframe
        loadingScale <- max(max(abs(sampleScore$axis1)) / max(abs(loadingValue$axis1)),
                            max(abs(sampleScore$axis2)) / max(abs(loadingValue$axis2))) * 0.85 # determine scaling for loading values
        loadingValuePlot <- loadingValue * loadingScale
        loadingValuePlot$lbl <- rownames(loadingValuePlot)

        # perpare for the axis labels
        varpp_biplot <- varpp[c(paste0("PC", as.character(input$PC_1st)), paste0("PC", as.character(input$PC_2nd)))] # extract the proportion of variance for the selected PCs
        pc_axis_lbl <- paste(c(paste0("PC", as.character(input$PC_1st)), paste0("PC", as.character(input$PC_2nd))), " (", round(varpp_biplot, digits = 2), "%)", sep = "")

        biplt <- ggplot(sampleScore, aes(x = axis1, y = axis2)) +
          geom_point(aes(shape = Group, colour = Group), size = input$biplot.SymbolSize) + # plot the sample score scatter plot
          ggtitle(input$biplot.Title) +
          xlab(pc_axis_lbl[1]) +
          ylab(pc_axis_lbl[2]) +
          theme_bw() +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
                plot.title = element_text(face = "bold", family = input$fontType, hjust = 0.5),
                axis.title = element_text(face = "bold", family = input$fontType),
                legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
                axis.text.x = element_text(size = input$xTickLblSize, family = input$fontType),
                axis.text.y = element_text(size = input$yTickLblSize, family = input$fontType, hjust = 0.5))

        if (input$ellipse){ # circles
          validate(need(min(summary(sampleScore$Group)) > 3, "Error: \n
                        A minimum of end = 4 is needed for ellipse.\n")) # end count check

          biplt <- biplt +
            stat_ellipse(aes(colour = Group, group = Group), type = "norm", level = input$ellipse_conf)
        }

        if (input$loadingplot){ # superimpose loading plot
          biplt <- biplt +
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_text(data = loadingValuePlot, aes(x = axis1, y = axis2), label = loadingValuePlot$lbl,
                      colour = "gray30", size = input$loadingSize)
        }

        return(biplt)
      })

      observe({
        output$PCAbiplot <- renderPlot({
          grid.draw(ggplotdata_biplot())
        }, width = input$biplot.Width, height = input$biplot.Height)
      })

      output$biplot_dlPlot <- downloadHandler(
        filename = function(){paste(substr(noquote(input$file1), 1, nchar(input$file1) - 4),".biplot.pdf", sep = "")},
        content = function(file) {
          ggsave(file, plot = ggplotdata_biplot(),
                 width = (input$biplot.Width * 25.4) / 72, height = (input$biplot.Height * 25.4) / 72,
                 units = "mm", dpi = 600, device = "pdf")
        }
      )

      # clear button
      observeEvent(input$clear_ui, {
        output$PCAboxplot <- renderPlot({})
      })

      observeEvent(input$clear_ui, {
        output$PCAbiplot <- renderPlot({})
      })

      # stop and close window
      observe({
        if (input$close > 0) stopApp()  # stop shiny
      })
    }
  )
  runApp(app, launch.browser = TRUE)
}
