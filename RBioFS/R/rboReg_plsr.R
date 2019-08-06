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
#' rbioReg_plsr(x, y, ncomp = 20)
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
#' rbioReg_plsr_ncomp_select(new_model,  multi_plot.ncol = 2, multi_plot.nrow = 2, plot.optm.ncomp.line = T,
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


#' @title rbioReg_plsr_perm()
#'
#' @description Permutation test for PLSR models.
#' @param object A \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.
#' @param ncomp Model complexity, i.e. number of components to model. Default is \code{object$ncomp}, i.e. maximum complexity.
#' @param adjCV If to use adjusted CV, i.e. CV adjusted for unbalanced data. Default is \code{FALSE}.
#' @param nperm Number of permutations to run. Default is \code{999}.
#' @param perm.plot Wether to produce a plot or not. Default is \code{TRUE}.
#' @param ... Additional argument for \code{\link{rbioUtil_perm_plot}}.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return The function returns \code{CSV} files for all intermediate permutation RMSEP values as well as the p-value resutls.
#'
#' The results are also exported to the environment as a \code{rbiomvr_perm} object with the following items:
#'
#' \code{nperm} The number of permutation runs.
#'
#' \code{perm.stats} The stats metric used for the permutation test.
#'
#' \code{adjCV} If the adjusted cross-validation stats is used.
#'
#' \code{p.value.summary} P values for the permutation test.
#'
#' \code{poerm.results} The intermediate permutation results, i.e. stats for each permutation test run in a data.frame. \code{nperm = 0} is the original stats.
#'
#' A scatter plot is also generaeted when \code{perm.plot = TRUE}.
#'
#' @details The function uses RMSEP as the stats for comparing original model with permutatsions.
#'
#'          Usually, we use the optimized PLS-DA model for \code{object}, which can be obtained from functions \code{\link{rbioClass_plsda}} and \code{\link{rbioClass_plsda_ncomp_select}}.
#'
#'          Data for permutation are object$centerX$centerX, meaning centered X are used if applicable.
#'
#'          Permutation methods are according to:
#'
#'             Ojala M, Garriga GC. 2010. Permutation test for studying classifier performance.
#'             J Mach Learn Res. 11: 1833 - 63.
#'
#'          For this function, labels (i.e. y) are permutated. A non-signifianct model (permutation p value > alpha, i.e. 0.05) in this case means the data is independent from the groups.
#'
#'
#' @import ggplot2
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pls RMSEP
#' @examples
#' \dontrun{
#' rbioReg_plsr_perm(object = new_model_optm, nperm = 999, adjCV = TRUE, parallelComputing = TRUE)
#' }
#' @export
rbioReg_plsr_perm <- function(object, ncomp = object$ncomp, adjCV = FALSE,
                              perm.plot = TRUE, ...,
                              nperm = 999,
                              parallelComputing = TRUE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                              verbose = TRUE){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr"))) stop("object has to be a \"rbiomvr\" class.")
  if (!"validation" %in% names(object) || is.null(object$validation)) stop("PLS-DA model has to include Cross-Validation.")
  if (!all(ncomp %in% seq(object$ncomp))) stop("ncomp contain non-existant comp.")
  if (length(nperm) != 1) stop("nperm can only contain one integer.")
  if (nperm %% 1 != 0) stop("nperm can only be integer. \n")
  if (nperm < 1) stop("nperm can only take interger equal to or greater than 1.")
  if (parallelComputing){
    clusterType <- match.arg(clusterType, c("PSOCK", "FORK"))
  }
  if (object$model.type != "regression") stop("object needs to have model.type = \"regression\"")

  ## calcuate RMSEP and construct original RMSEP data frame
  rmsep <- pls::RMSEP(object)
  rmsep_dfm <- foreach(i = 1:dim(rmsep$val)[2], .combine = "rbind") %do% {
    dfm <- as.data.frame(t(rmsep$val[, i, ncomp]))
    if (adjCV){
      stats <- dfm[, "adjCV"]
    } else {
      stats <- dfm[, "CV"]
    }
    dfm <- data.frame(comparison = dimnames(rmsep$val)[[2]][i], comps = object$ncomp, RMSEP = stats,
                      adjCV = adjCV, row.names = NULL)
  }

  ## permutation test
  # permutation test functions
  by_y_func <- function(i){
    set.seed(i)
    perm_y <- object$inputY[sample(1:length(object$inputY))]  # sample label permutation
    perm_model <- rbioReg_plsr(x = object$centerX$centerX, y = perm_y,
                               ncomp = ncomp, scale = FALSE, validation = object$validation_method,
                               segments = length(object$validation$segments),
                               segments.type = object$out_model$validation_segments_type,
                               verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
    perm_rmsep <- pls::RMSEP(perm_model)  # permutation model RMSEP

    perm_rmsep_dfm <- foreach(j = 1:dim(perm_rmsep$val)[2], .combine = "rbind") %do% {
      dfm <- as.data.frame(t(perm_rmsep$val[, j, ncomp + 1]))
      if (adjCV){
        stats <- dfm[, "adjCV"]
      } else {
        stats <- dfm[, "CV"]
      }
      dfm <- data.frame(nperm = i, comparison = dimnames(perm_rmsep$val)[[2]][j],
                        comps = ncomp, RMSEP = stats, adjCV = adjCV, row.names = NULL)
    }
    return(perm_rmsep_dfm)
  }

  if (verbose) cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
  if (verbose) cat(paste0("Running permutation test using ", perm.method, " method ","with ", nperm, " permutations (speed depending on hardware configurations)..."))
  # consolidated permutation model RMSEP data frame
  if (!parallelComputing){
    perm_dfm <- foreach(i = 1:nperm, .combine = "rbind") %do% by_y_func(i)
  } else {  # parallel computing
    # set up cpu cluster
    n_cores <- n_cores
    cl <- makeCluster(n_cores, type = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # permutation test
    perm_dfm <- foreach(i = 1:nperm, .combine = "rbind", .packages = c("foreach", "pls", "RBioFS")) %dopar% by_y_func(i)
  }

  # calculate total perm RMSEP
  perm_tot_rmsep_dfm <- foreach(i = 1:nperm, .combine = "rbind") %do% {
    subdfm <- perm_dfm[perm_dfm$nperm == i, , drop = FALSE]
    # below: RMSEP is total RMSEP across Y
    tot_rmsep <- data.frame(nperm = i, comparison = "Y", comps = ncomp, RMSEP = sum(subdfm$RMSEP), adjCV = adjCV)
  }

  # permutation model data frame
  perm_results_p_val <- foreach(i = dimnames(rmsep$val)[[2]], .combine = "rbind") %do% {
    tmp <- perm_tot_rmsep_dfm[perm_tot_rmsep_dfm$comparison == i, ]
    orig_rmsep <- data.frame(nperm = 0, rmsep_dfm[rmsep_dfm$comparison == i, ], row.names = NULL)
    tmp <- rbind(orig_rmsep, tmp)
    rownames(tmp) <- NULL
    write.csv(file = paste0(deparse(substitute(object)),"_", i ,"_vs_rest" , "_perm_intermediate.csv"), tmp, row.names = FALSE)  # export full permutation results to working directory
    orig.val <- tmp$RMSEP[which(tmp$nperm == 0)]
    perm.val <- tmp$RMSEP[which(tmp$nperm != 0)]
    p.val <- (length(which(perm.val <= orig.val)) + 1) / (nperm + 1)
    out <- data.frame(comparison = i, original.RMSEP = orig.val, p.value = p.val,
                      row.names = NULL, stringsAsFactors = FALSE)
  }
  if (verbose) cat("Done!\n")

  ## output
  if (verbose) cat("\n")
  if (verbose) cat("Permutation test restuls: \n")
  if (verbose) print(perm_results_p_val)
  if (verbose) cat("\n")

  # env
  perm_stats <- rbind(data.frame(nperm = rep(0, times = nrow(rmsep_dfm)), comparison = rmsep_dfm$comparison, RMSEP = rmsep_dfm$RMSEP, row.names = NULL),
                      data.frame(nperm = perm_tot_rmsep_dfm$nperm, comparison = perm_tot_rmsep_dfm$comparison, RMSEP = perm_tot_rmsep_dfm$RMSEP, row.names = NULL))

  out <- list(perm.method = "by_y", nperm = nperm, perm.stats = "RMSEP", adjCV = adjCV,
              p.value.summary = perm_results_p_val, perm.results =  perm_stats, model.type = object$model.type)

  class(out) <- "rbiomvr_perm"
  assign(paste(deparse(substitute(object)), "_perm", sep = ""), out, envir = .GlobalEnv)

  # directory
  if (verbose) cat(paste0("Permutation RMSEP values stored in csv files with suffix: _perm_intermediate.csv. \n"))
  if (verbose) cat(paste0("Permutation results stored in csv file: ", paste(deparse(substitute(object)), "_perm.csv", sep = ""), ". \n"))
  write.csv(file = paste(deparse(substitute(object)), "_perm.csv", sep = ""), perm_results_p_val, row.names = FALSE)

  ## plot
  if (perm.plot){
    plsda_permutation_test <- out
    rbioUtil_perm_plot(perm_res = plsda_permutation_test, ..., verbose = verbose)
  }
}


#' @title rbioReg_plsr_vip
#'
#' @description VIP, or variable importance in projection, calcualtion and plotting for PLSR models. This is another FS method, and can be used independently.
#' @param object A \code{mvr} or \code{rbiomvr} object. Make sure the model is built uisng \code{"oscorespls"} method.
#' @param vip.alpha Alpha value (threshold) for VIP values. Any VIP above this is considered important. Defaults is \code{1}.
#' @param comps Integer vector. Components to plot. The index of the components are intergers. The vector length should be between 1 and the total number of components, inclusive. Default is \code{c(1, 2)}.
#' @param bootstrap If to use boostrap for VIP calculation, so that standard deviation on VIP can be estimated. Default is \code{TRUE}.
#' @param boot.n Set only when \code{boostrap = TRUE}, nummbers of iterations for boostrap. Default is \code{50}.
#' @param boot.parallelComputing Set only when \code{boostrap = TRUE}, if to use parallel computering for bootstrap process. Default is \code{TRUE}.
#' @param boot.n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param boot.clusterType Set only when \code{boostrap = TRUE} and \code{boot.parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param plot If to generate a plot. Default is \code{TRUE}.
#' @param ... Additional arguments to \code{\link{rbioFS_plsda_vip_plot}}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Outputs a \code{rbiomvr_vip} class object to the environment.
#'
#'         \code{rbiomvr_vip} items are:
#'
#'         \code{vip_summary}
#'
#'         \code{comps}: number of PLSR components calculated fro VIP
#'
#'         \code{vip.alpha}: cutoff VIP value for a feature to be considered as important
#'
#'         \code{features_above_alpha}: selected features according to vip.alpha, i.e. features with VIP >= vip.alpha
#'
#'         \code{boostrap}
#'
#'         \code{boot.n}: boostrap iterations
#'
#'         \code{bootstrap.iteration.results}: raw VIP values for each bootstrap iteration
#'
#'         \code{mode.type}
#'
#' @details Only works when the plsda model is fitted with the orthorgonal score algorithm, or NIPALS. Such model can be built using \code{\link{rbioClass_plsda}} with \code{method = "oscorespls"}.
#'
#' The \code{vip.alpha} of 1 is the most commonly accepted value. However it is also acceptable to set according to the data and objectives of the study.
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' rbioReg_plsr_vip(object = new_model_optm, boostrap = TRUE, boot.m = 50,
#'                  vip.alpha = 0.8,
#'                  plot = TRUE, plot.title = TRUE, plot.titleSize = 10,
#'                  plot.sig.line = TRUE,
#'                  plot.outlineCol = "black", plot.errorbar = "SEM", plot.errorbarWidth = 0.2,
#'                  plot.errorbarLblSize = 6, plot.fontType = "sans", plot.xLabel = "Features",
#'                  plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
#'                  plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 1, plot.xvAligh = 0.2,
#'                  plot.rightsideY = TRUE,
#'                  plot.yLabel = "Coefficients", plot.yLabelSize = 10, plot.yTickLblSize = 10,
#'                  plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
#'                  plot.legendTtl = FALSE, plot.legendTtlSize = 9, plot.Width = 170,
#'                  plot.Height = 150)
#' }
#' @export
rbioReg_plsr_vip <- function(object, vip.alpha = 1, comps = c(1, 2),
                             bootstrap = TRUE,
                             boot.n = 50, boot.parallelComputing = TRUE, boot.n_cores = parallel::detectCores() - 1, boot.clusterType = "PSOCK",
                             plot = TRUE,...,
                             verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c("rbiomvr", "mvr"))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.")
  if (object$model.type != "regression") stop("object needs to have model.type = \"regression\"")
  if (object$method != "oscorespls") stop("Object needs to fit using oscorespls (i.e. NIPALS) algorithm. Please re-fit using rbioClass_plsda with method = \"oscorespls\".") # only oscorespls algorithm is applicable to VIP
  if (length(comps) > object$ncomp)stop("comps length exceeded the maximum comp length.")
  if (!all(comps %in% seq(object$ncomp)))stop("comps contain non-existant comp.")
  if (bootstrap & (boot.n < 2 | boot.n %% 1 != 0)) stop("when boostrap = TRUE, boot.n needs to be an integer greater than 1.")

  ## VIP process
  if (verbose) cat(paste0("Boostrap: ", ifelse(bootstrap, "ON\n", "OFF\n")))
  if (bootstrap){  # bootstrap VIP process
    # message
    if (verbose) cat(paste0("Parallel computing: ", ifelse(boot.parallelComputing, "ON\n", "OFF\n")))
    if (verbose) cat(paste0("Boostrap interation: ", boot.n, "\n"))
    if (verbose) cat("Boostrap VIP calculation (speed depending on hardware configuration)...")

    # set up data
    orig.dat <- data.frame(Y = object$inputY, object$inputX, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    sum.boot.vip_raw_list <- vector(mode = "list", length = boot.n)

    # set up boot function
    boot.func <- function(i){
      # bootstrap resampling
      idx <- seq(nrow(orig.dat))
      boot.idx <- sample(group.idx, replace = TRUE)
      boot.dat <- orig.dat[boot.idx, ]

      # bootstrap modelling
      boot.X <- boot.dat[, -1]
      boot.Y <- boot.dat[, 1]
      boot.m <- rbioReg_plsr(x = boot.X, y = boot.Y,
                             ncomp = object$ncomp, validation = "CV", method = "oscorespls",
                             verbose = FALSE)

      # boostrap VIP
      boot.score <- boot.m$scores
      boot.lodw <- boot.m$loading.weights
      boot.Wnorm2 <- colSums(boot.lodw^2)

      # vip_raw_list contains VIP for all the components
      boot.vip_raw_list <- vector(mode = "list", length = dim(boot.m$Yloadings)[1])
      boot.vip_raw_list[] <- foreach(m = 1:dim(boot.m$Yloadings)[1]) %do% {
        ylod <- boot.m$Yloadings[m, ] # one y variable at a time
        boot.ss <- c(ylod)^2 * colSums(boot.score^2)
        boot.ssw <- sweep(boot.lodw^2, 2, boot.ss / boot.Wnorm2, "*")
        boot.vip <- sqrt(nrow(boot.ssw) * apply(boot.ssw, 1, cumsum) / cumsum(boot.ss))  # cumsum: cucmulative sum
        if (object$ncomp == 1) {
          boot.vip <- matrix(boot.vip, nrow = 1, dimnames = list('Comp 1', names(boot.vip)))
        }
        boot.vip <- boot.vip[comps, , drop = FALSE]
        return(boot.vip)
      }
      names(boot.vip_raw_list) <- dimnames(boot.m$Yloadings)[[1]]
      boot.vip_raw_list
    }

    # bootstrap modelling and VIP calculation
    if (!boot.parallelComputing){  # non-parallel
      sum.boot.vip_raw_list[] <- foreach(i = 1:boot.n) %do% boot.func(i)
    } else {  # boot parallel computing
      # set up cluster
      n_cores <- boot.n_cores
      cl <- makeCluster(n_cores, type = boot.clusterType)
      registerDoParallel(cl)
      on.exit(stopCluster(cl)) # close connect when exiting the function

      # computing
      sum.boot.vip_raw_list[] <- foreach(i = 1:boot.n, .packages = c("foreach", "RBioFS")) %dopar% boot.func(i)
    }
    names(sum.boot.vip_raw_list) <- paste0("boot ", seq(boot.n))

    # sum.boot.vip_raw_list structure: boot$group$: (row)`comp 1`|`comp 2`|`comp 3`...
    # targeted format: group$`comp 1`: (colums)`boot 1`|`boot 2`|`boot 3`...
    group.comp.boot.vip_list <- vector(mode = "list", length = 1)
    group.comp.boot.vip_list[] <- foreach(i = 1) %do% {
      group.comp.list <- vector(mode = "list", length = length(comps))
      group.comp.list[] <- foreach(j = comps) %do% {
        group.comp.boot_mtx <- foreach(m = 1:boot.n, .combine = "cbind") %do% {
          sum.boot.vip_raw_list[[m]][[i]][j, ]  # i: groups, j: boot, 1: comp
        }
        colnames(group.comp.boot_mtx) <- paste0("boot ", seq(boot.n))
        group.comp.boot_mtx
      }
      names(group.comp.list) <- paste0("comp ", comps)
      group.comp.list
    }
    names(group.comp.boot.vip_list) <- levels(object$inputY)

    # plot list
    boot.vip.plt_dat_list <- vector(mode = "list", length = length(levels(object$inputY)))
    boot.vip.plt_dat_list[] <- foreach(i = 1:length(levels(object$inputY))) %do% {
      boot.plt_list.comp <- vector(mode = "list", length = length(comps))
      boot.plt_list.comp[] <- foreach(j = comps) %do% {
        boot.mean <- rowMeans(group.comp.boot.vip_list[[i]][[j]])
        boot.sd <- matrixStats::rowSds(group.comp.boot.vip_list[[i]][[j]])
        boot.sem <- matrixStats::rowSds(group.comp.boot.vip_list[[i]][[j]])/sqrt(boot.n)
        outdfm <- data.frame(Features = rownames(group.comp.boot.vip_list[[i]][[j]]),
                             VIP = boot.mean,
                             sd = boot.sd,
                             sem = boot.sem,
                             row.names = NULL,
                             check.names = FALSE)
        outdfm <- outdfm[order(outdfm$VIP, decreasing = TRUE),]
        outdfm$Features <- factor(outdfm$Features, levels = unique(outdfm$Features))
        rownames(outdfm) <- NULL
        outdfm
      }
      names(boot.plt_list.comp) <- paste0("comp ", comps)
      boot.plt_list.comp
    }
    names(boot.vip.plt_dat_list) <- levels(object$inputY)
    vip_list <- boot.vip.plt_dat_list
  } else {  # non bootstrap VIP process
    # message
    if (verbose) cat("VIP calculation...")

    # VIP calculation
    score <- object$scores
    lodw <- object$loading.weights
    Wnorm2 <- colSums(lodw^2)

    # vip_raw_list contains VIP for all the components
    vip_raw_list <- vector(mode = "list", length = dim(object$Yloadings)[1])
    vip_raw_list[] <- foreach(i = 1:dim(object$Yloadings)[1]) %do% {
      ylod <- object$Yloadings[i, ] # one y variable at a time
      ss <- c(ylod)^2 * colSums(score^2)
      ssw <- sweep(lodw^2, 2, ss / Wnorm2, "*")
      vip <- sqrt(nrow(ssw) * apply(ssw, 1, cumsum) / cumsum(ss))  # cumsum: cucmulative sum
      vip <- vip[comps, , drop = FALSE]
      vip
    }
    names(vip_raw_list) <- dimnames(object$Yloadings)[[1]]

    vip_list <- vector(mode = "list", length = 1)
    vip_list[] <- foreach(i = 1) %do% {
      plt_list.comp <- vector(mode = "list", length = length(comps))
      plt_list.comp[] <- foreach(j = comps) %do% {
        vip <- vip_raw_list[[i]][j, ]
        outdfm <- data.frame(Features = colnames(vip_raw_list[[i]]),
                             VIP = vip,
                             row.names = NULL,
                             check.names = FALSE)
        outdfm <- outdfm[order(outdfm$VIP, decreasing = TRUE),]
        outdfm$Features <- factor(outdfm$Features, levels = unique(outdfm$Features))
        rownames(outdfm) <- NULL
        outdfm
      }
      names(plt_list.comp) <- paste0("comp ", comps)
      plt_list.comp
    }
    names(vip_list) <- "Y"
  }
  if (verbose) cat("Done!\n")

  ## output
  # important features lsit
  final.ipf_list <- vector(mode = "list", length = length(vip_list))
  final.ipf_list <- foreach(m = 1:length(vip_list)) %do% {
    ipf_list <- vector(mode = "list", length = length(vip_list[[m]]))
    ipf_list[] <- foreach(n = 1:length(vip_list[[m]])) %do% {
      as.character(vip_list[[m]][[n]][which(vip_list[[m]][[n]]$VIP > vip.alpha), 1])
    }
    names(ipf_list) <- names(vip_list[[m]])
    ipf_list
  }
  names(final.ipf_list) <- names(vip_list)

  # output
  out <- list(vip_summary = vip_list,
              vip.alpha = vip.alpha,
              features_above_alpha = final.ipf_list,
              comps = comps,
              bootstrap = bootstrap,
              boot.n = if (bootstrap) boot.n else NULL,
              bootstrap.iteration.results = if (bootstrap) group.comp.boot.vip_list else NULL,
              model.type = object$model.type)
  class(out) <- "rbiomvr_vip"
  assign(paste(deparse(substitute(object)), "_plsda_vip", sep = ""), out, envir = .GlobalEnv)

  ## plot
  if (plot){
    # message
    if (verbose) cat("\n")
    RBioFS::rbioReg_plsr_vip_plot(vip_obj = out, ...)
  }
}


