#' @title rbioReg_plsr
#'
#' @description PLS regression modelling
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is a continuous vairable..
#' @param method PLSR modelling method. Four PLSR algorithms are available: the kernel algorithm ("kernelpls"), the wide kernel algorithm ("widekernelpls"),
#'               SIMPLS ("simpls") and the classical orthogonal scores algorithm (also known as NIPALS) ("oscorespls"). Default is the popular \code{"simpls"}.
#' @param scale Logical, whether to scale the data or not. Default is \code{TRUE}.
#' @param validation Cross validation methods. Options are "none", "CV" (fold), "LOO" (leave-one-out). Default is \code{"CV"}.
#' @param segments Set only when \code{validation = "CV"}, the number of segement to be set. Default is \code{10}.
#' @param segments.type Method to set up the segments. Options are \code{"random", "consecutive", "interleaved"}. Default is \code{"random"}.
#' @param jackknife If to use jack-knife procedure. Default is \code{TRUE}.
#' @param ... Additional arguments for \code{mvr} function from \code{pls} pacakge.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details The center and scale processes are handled by \code{\link{center_scale}} function.
#'          Therefore, the \code{center} and \code{scale} arguments \code{pls::mvr} are set to \code{FALSE}.
#'
#'          The function always centers the data.
#'
#' @return Returns a PLS-DA model object, with classes "mvr" and "rbiomvr".
#'
#'        Additional items for \code{rbiomvr} object to \code{mvr} object from pls package:
#'
#'        \code{centerX}: centered X data if applicable, with scale option.
#'
#'        \code{dummy_y}: NULL
#'
#'        \code{inputX}: raw input predictor data.
#'
#'        \code{inputY}: input group labels.
#'
#'        \code{validation_method}
#'
#'        \code{validation_segments_type}
#'
#'        \code{model.type}
#'
#' @details The data is always centered prior to modelling, i.e. \code{x - col.mean}. Thus no "center = TRUE/FALSE" option is provided.
#'
#'
#'          For sequencing data, the input x needs to be either tranformed with the function like \code{clr_ilr_transfo()} from \code{RBioArray} package,
#'          or normalized using methods like "TMM" or "RLE" implemented in \code{edgeR} pacakge.
#'
#' @importFrom pls plsr
#' @examples
#' \dontrun{
#' rbioClass_plsda(x, y, ncomp = 20)
#' }
#' @export
rbioReg_plsr <- function(x, y, method = "simpls",
                         scale = TRUE,
                         validation = c("none", "CV", "LOO"),
                         segments = 10, segments.type = "random",
                         jackknife = TRUE, ...,
                         verbose = TRUE){
  ## check arguments
  if (!class(x) %in% c("data.frame", "matrix") & !is.null(dim(x))) stop("x needs to be a matrix, data.frame or vector.")
  if (class(x) == "data.frame" | is.vector(x)){
    if (verbose) cat("x converted to a matrix object.\n")
    x <- as.matrix(sapply(x, as.numeric))
  }

  ## data processing
  # X
  if (verbose) cat(paste0("data centered with the scale option ", ifelse(scale, "\"ON\" ", "\"OFF\" "), "prior to modelling..."))
  centered_X <- center_scale(x, scale = scale)  # center data with the option of scaling
  if (verbose) cat("DONE!\n")
  X <- centered_X$centerX
  # Y
  Y <- y
  # constract dataframe for modelling
  model_dfm <- data.frame(n = paste("row", 1:length(Y), sep = ""))
  model_dfm$Y <- Y  # y is the y matrix as a whole, not the content of y
  model_dfm$X <- X  # x is the x matrix as a whole, not the content of x

  ## modelling
  if (validation == "LOO") {
    segments.type <- NA
    out_model <- plsr(Y ~ X, data = model_dfm,
                      method = method, scale = FALSE, center = FALSE,
                      validation = validation, jackknife = TRUE, ...)
  } else {
    out_model <- plsr(Y ~ X, data = model_dfm,
                      method = method, scale = FALSE, center = FALSE,
                      validation = validation, segments = segments, segments.type = segments.type, jackknife = TRUE, ...)
  }
  out_model$centerX <- centered_X
  out_model$dummy_y <- NULL
  out_model$inputX <- x
  out_model$inputY <- y
  out_model$validation_method <- validation
  out_model$validation_segments_type <- segments.type
  out_model$model.type <- "regression"

  class(out_model) <- c("rbiomvr", "mvr")
  return(out_model)
}


#' @title rbioReg_plsr_ncomp_select()
#'
#' @description Optimal number of components selection for PLSR model, with RMSEP plot funcitonality. Selection methods are modified based on \code{selectNcomp()} from \code{pls} pacakge.
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param ... Additional argument for \code{RMSEP} function from \code{pls} package.
#' @param ncomp.selection.method Optimal numbers of components selection method. Options are \code{"min"}, \code{"1err"}, and \code{"randomization"}. Default is \code{"1sd"}.
#' @param randomization.nperm Set only when \code{ncomp.selection.method = "randomization"}, number of permutations. Default is \code{999}.
#' @param randomization.alpha Set only when \code{ncomp.selection.method = "randomization"}, alpha for the p values used during "randomization" selection. Default is \code{0.05}.
#' @param rmsepplot If to generate a RMSEP plot. Default is \code{TRUE}.
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}. Note: doesn't seem to be necessasry as PLS-DA always has at least two y classes.
#' @param plot.optm.ncomp.line If to display the vertical line indicting the optimal number of components. Default is \code{TRUE}.
#' @param multi_plot.ncol Number of columns on one figure page. Default is the number of responding classes, i.e. y.
#' @param multi_plot.nrow Number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"Components"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"RMSEP"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Plot width. Default is \code{170}.
#' @param plot.Height Plot height. Default is \code{150}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Prints the selected number of components for each y class. Returns RMSEP values for each y class to the environment as a \code{rbiomvr_ncomp_select} object,
#'         as well as a pdf file for the RMSEP plot if \code{rmsepplot = TRUE}.
#' @details The RMSEP figure shows both CV estimates and adjusted CV estimates, which is CV estimiates corrected for bias.
#'          Three methods are used for components number selection: \code{"min"} simply chooses the number of components to
#'          reach te minimum RMSEP; \code{"1err"} chooses the number of components when its RMSEP first reaches minimum as well as within one standard error;
#'          For "randomization", see the help file for \code{selectNcomp()} function from  \code{pls} pacakge.
#' @import ggplot2
#' @import foreach
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @importFrom pls RMSEP
#' @examples
#' \dontrun{
#' rbioClass_plsda_ncomp_select(new_model,  multi_plot.ncol = 2, multi_plot.nrow = 2, plot.optm.ncomp.line = T,
#'          ncomp.selection.method = "randomization", randomization.nperm = 999, randomization.alpha = 0.05)
#' }
#' @export
rbioReg_plsr_ncomp_select <- function(object, ...,
                                      ncomp.selection.method = c("1err", "min", "randomization"), randomization.nperm = 999, randomization.alpha = 0.05,
                                      rmsepplot = TRUE,
                                      plot.rightsideY = TRUE,
                                      plot.optm.ncomp.line = TRUE,
                                      multi_plot.ncol = length(dimnames(object$coefficients)[[2]]), multi_plot.nrow = 1,
                                      multi_plot.legend.pos = c("bottom", "top", "left", "right"),
                                      plot.display.Title = TRUE,
                                      plot.SymbolSize = 2,
                                      plot.fontType = "sans",
                                      plot.xLabel = "Components", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                      plot.yLabel = "RMSEP", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                      plot.legendSize = 9,
                                      plot.Width = 170, plot.Height = 150,
                                      verbose = TRUE){
  ## check arguments
  if (!any(class(object) %in% c("mvr", "rbiomvr"))) stop("object has to be a mvr or rbiomvr class.")
  if (!"validation" %in% names(object)) stop("PLS-DA model has to include Cross-Validation.")
  if (object$model.type != "regression") stop("The input object needs to have model.type = \"regression\".")
  # if (!tolower(ncomp.selection.method) %in% c("min", "1err", "randomization")) stop("ncomp.selection.method needs to be \"min\", \"1err\", or \"randomization\" exactly.")
  ncomp.selection.method <- match.arg(tolower(ncomp.selection.method), c("1err", "min", "randomization"))
  multi_plot.legend.pos <- match.arg(tolower(multi_plot.legend.pos), c("bottom", "top", "left", "right"))

  ## calcuate RMSEP
  rmsep <- RMSEP(object, ...)

  ## ncomp selection and plot
  # prepare dataframe for ncomp selection and plotting
  rmsep_dfm_list <- vector(mode = "list", length = dim(rmsep$val)[2])
  plt_list <- vector(mode = "list", length = dim(rmsep$val)[2])
  rmsep_dfm_list[] <- foreach(i = 1:dim(rmsep$val)[2]) %do% {
    dfm <- as.data.frame(t(rmsep$val[,i,]))
    dfm <- data.frame(comps = seq(0, dim(rmsep$val)[3] - 1), dfm, row.names = NULL)
    dfm <- melt(dfm, id = "comps")
  }
  names(rmsep_dfm_list) <- dimnames(rmsep$val)[[2]]  # dimnames(tstrmsep$val)[2] is a list object. use [] to get the characters in string class.

  # ncomp selection
  ncompsel_mtx <- foreach(i = 1:length(rmsep_dfm_list), .combine = "rbind") %do% {
    dfm_cv <- rmsep_dfm_list[[i]][rmsep_dfm_list[[i]]$variable == "CV", ]
    dfm_adjcv <- rmsep_dfm_list[[i]][rmsep_dfm_list[[i]]$variable == "adjCV", ]

    if (ncomp.selection.method == "min"){
      cv_x <- dfm_cv$comps[which.min(dfm_cv$value)]
      adjcv_x <- dfm_adjcv$comps[which.min(dfm_adjcv$value)]
    } else {  # modified from pls::selectNcomp()
      origResponse <- object$inputY
      allresids <- cbind(origResponse - (sum(origResponse) - origResponse) / (length(origResponse) - 1),
                         object$validation$pred[,i,] - origResponse)  # CV prediction for all the comps used
      residsds <- apply(allresids, 2, sd) / sqrt(nrow(allresids))

      if (ncomp.selection.method == "1err"){
        dfm_cv$SD <- residsds
        min_rmsep_cv_idx <- which.min(dfm_cv$value)  # index for the minimum mean oob feature group
        sd_min_cv <- dfm_cv$SD[min_rmsep_cv_idx]  # oob SD for the feature group above
        cv_x <- min(with(dfm_cv, which(value <= (value[min_rmsep_cv_idx] + sd_min_cv))))  # 1sd minimum selection

        dfm_adjcv$SD <- residsds
        min_rmsep_adjcv_idx <- which.min(dfm_adjcv$value)  # index for the minimum mean oob feature group
        sd_min_adjcv <- dfm_adjcv$SD[min_rmsep_adjcv_idx]  # oob SD for the feature group above
        adjcv_x <- min(with(dfm_adjcv, which(value <= (value[min_rmsep_adjcv_idx] + sd_min_adjcv))))  # 1sd minimum selection
      } else if (ncomp.selection.method == "randomization"){
        absBest_cv <- which.min(dfm_cv$value)
        pvals_cv <- sapply(seq_len(absBest_cv - 1), function(m) randomiz.test(allresids[, m], allresids[, absBest_cv], nperm = randomization.nperm))
        idx_cv <- which(pvals_cv > randomization.alpha)
        cv_x <- min(c(idx_cv, absBest_cv)) - 1

        absBest_adjcv <- which.min(dfm_adjcv$value)
        pvals_adjcv <- sapply(seq_len(absBest_adjcv - 1), function(m) randomiz.test(allresids[, m], allresids[, absBest_adjcv], nperm = randomization.nperm))
        idx_adjcv <- which(pvals_cv > randomization.alpha)
        adjcv_x <- min(c(idx_adjcv, absBest_adjcv)) - 1
      }
    }
    c(cv_x, adjcv_x)
  }
  ncompsel_mtx <- matrix(ncompsel_mtx, ncol = 2)
  rownames(ncompsel_mtx) <- names(rmsep_dfm_list)
  colnames(ncompsel_mtx) <- c("CV", "adjCV")
  if (verbose) cat(paste0("Optimial ncomp (selection method: ", ncomp.selection.method, ") : \n", sep = ""))
  if (verbose) print(ncompsel_mtx)

  # plotting
  if (rmsepplot){
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.rmsepplot.pdf...", sep = ""))  # initial message
    plt_list[] <- foreach(i = 1:length(rmsep_dfm_list)) %do% {
      plt <- ggplot(data = rmsep_dfm_list[[i]], aes(x = comps, y = value, colour = variable)) +
        geom_line(aes(linetype = variable)) +
        geom_point(aes(shape = variable), size = plot.SymbolSize) +
        ggtitle(ifelse(plot.display.Title, names(rmsep_dfm_list)[i], NULL)) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              legend.position = "bottom", legend.key = element_blank(),
              legend.text = element_text(size = plot.legendSize), legend.title = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
              axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))

      if (plot.optm.ncomp.line){
        plt <- plt +
          geom_vline(xintercept = ncompsel_mtx[i, ], linetype = "dashed", colour = "red")
      }

      if (plot.rightsideY & length(rmsep_dfm_list) == 1){
        plt <- rightside_y(plt)
      }
      plt
    }
    names(plt_list) <- paste0("g_", names(rmsep_dfm_list))

    if (length(rmsep_dfm_list) > 1) {
      if (multi_plot.ncol * multi_plot.nrow < length(rmsep_dfm_list)){
        stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. make sure they match the number of response groups.")
      } else {
        plt <- RBioplot::multi_plot_shared_legend(plt_list, ncol = multi_plot.ncol, nrow = multi_plot.nrow, position = multi_plot.legend.pos)
      }
    }

    ## save
    # grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.rmsepplot.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    if (verbose) cat("Done!\n")
  }

  ## return RMSEP values
  rmsep_dfm_list$ncomp_selected <- ncompsel_mtx
  class(rmsep_dfm_list) <- "rbiomvr_ncomp_select"
  assign(paste(deparse(substitute(object)), "_plsr_ncomp_select", sep = ""), rmsep_dfm_list, envir = .GlobalEnv)
}

