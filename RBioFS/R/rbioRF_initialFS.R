#' @title rbioRF_initialFS
#'
#' @description Recursive random froest variable importance (vi) and OOB error rate computation.
#' @param objTitle The title for the output data frame. Default is \code{"x_vs_tgt"}
#' @param x Input dataframe or matrix. Make sure to arrange the data with features as column names.
#' @param targetVar The target variable for random forest feature selection. This is a factor object.
#' @param nTimes Number of random forest vi computation runs. Default is \code{50} times.
#' @param nTree Number of trees generated for each random forest run. Default is \code{1001} trees.
#' @param mTry Number of random feature pick when building the tree. Default is \code{max(floor(ncol(dfm) / 3), 2)}.
#' @param multicore If to use parallel computing. Default is \code{TRUE}.
#' @param plot If to plot a bargraph to visualize vi and the ranking. Default is \code{TRUE}
#' @param Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param n Number of features to show. Takes integer numbers. Default is \code{"all"} (make sure to include quotation marks).
#' @param errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @return Outputs a \code{list} object with vi values for each feature and OOB error rate. When \code{TRUE}, bargraph for the vi is also generated and exported as a \code{.pdf} file.
#' @details Make sure to arrange data (dfm) with feature (e.g., gene) as variables (i.e., columns), and rownames as sample names.
#' @import ggplot2
#' @import foreach
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom randomForest randomForest importance
#' @importFrom rpart rpart prune
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @examples
#' \dontrun{
#' rbioRF_initialFS(training_HCvTC, tgtVar_HCvTC, n = 40, errorbar = "SEM", plotWidth = 400, plotHeight = 200)
#' }
#' @export
rbioRF_initialFS <- function(objTitle = "x_vs_tgt",
                             x, targetVar, nTimes = 50, nTree = 1001, mTry = max(floor(ncol(x) / 3), 2),
                             multicore = TRUE,
                             plot = TRUE, n = "all",
                             Title = NULL, xLabel = "Mean Decrease in Accuracy", yLabel = NULL,
                             errorbar = "SEM", errorbarWidth = 0.2,
                             xTxtSize = 10, yTxtSize =10,
                             plotWidth = 170, plotHeight = 150){

  #### check the variables
  if (ncol(x) == 1){
    stop("only one feature detected. No need to select.")
  }

  #### recursive RF
  ### load the dataframe/matrix
  training <- as.matrix(x)

  ### pepare the target variable
  tgt <- factor(as.character(targetVar), levels = unique(targetVar))

  ### prepare draw size. this uses down-sampling if the samples are unbalanced
  nlvl <- length(levels(tgt))
  size <- min(as.vector(table(tgt))) # down-sampling
  drawSize <- rep(size, nlvl)

  ### repeating random forest - recursive approach
  # pre-set an empty matrix with the number of columns same as the number of RF runs
  # note that nrow is the number of features, hence the ncol of the training set
  vimtx <- matrix(nrow = ncol(training), ncol = nTimes)
  errmtx <- matrix(nrow = 1, ncol = nTimes)

  if (!multicore){ ## signle core computing: recursive structure
    tmpFunc <- function(n, m, tmptimes, tmpvimtx, tmperrmtx, tmpTraining, tmpTgt,
                        tmpTree, tmpTry, tmpSize){

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
        tmpFunc(n - 1, m + 1, tmptimes, tmpvimtx, tmperrmtx, tmpTraining, tmpTgt,
                tmpTree, tmpTry, tmpSize)
      }
    }

    lst <- tmpFunc(n = nTimes, m = 1, tmptimes = nTree, tmpvimtx = vimtx, tmperrmtx = errmtx, tmpTraining = training, tmpTgt = tgt,
                   tmpTree = nTree, tmpTry = mTry, tmpSize = drawSize)

    phase0mtx_vi <- lst$raw_vi
    phase0mtx_OOB_err <- lst$raw_OOB_error

  } else { ## parallel computing
    # set up cpu cluster
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # recursive RF using par-apply functions
    tmpfunc2 <- function(i){
      rf <- randomForest::randomForest(x = training, y = tgt, ntree = nTree, mtry = mTry, importance = TRUE,
                                       proximity = TRUE, drawSize = drawSize)

      impt <- randomForest::importance(rf, type = 1)
      tmpvimtx <- impt[, 1] # fill the vi matrix
      tmperrmtx <- rf$err.rate[nTree, 1] # fill the OOB error rate
      lst <- list(tmpvimtx = tmpvimtx, tmperrmtx = tmperrmtx)
    }

    # foreach parallel
    tmp <- foreach(i = 1:nTimes) %dopar% tmpfunc2(i)
    vimtx <- foreach(i = 1:nTimes, .combine = cbind) %dopar% tmp[[i]]$tmpvimtx
    errmtx <- foreach(i = 1:nTimes, .combine = cbind) %dopar% tmp[[i]]$tmperrmtx

    rownames(vimtx) <- colnames(training)
    colnames(vimtx) <- c(paste("vi", seq(nTimes), sep = "_"))

    rownames(errmtx) <- "OOB_error_rate"
    colnames(errmtx) <- c(paste("OOB_error_tree", seq(nTimes), sep = "_"))

    phase0mtx_vi <- vimtx
    phase0mtx_OOB_err <- errmtx
  }

  ####prepare output vi and OOB error dataframes
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
  rownames(outdfm_OOB_err) <- paste(nTimes, "trees_OOB_err", sep = "_")

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

  ## vi plotting
  if (plot){
    loclEnv <- environment()

    # prepare plotting dataframe
    if (n != "all"){
      pltdfm <- tail(tmpdfm_vi, n)
    } else {
      pltdfm <- tmpdfm_vi
    }

    # plotting
    baseplt <- ggplot(pltdfm, aes(x = Target, y = Mean), environment = loclEnv) +
      geom_bar(position="dodge", stat="identity", color="black", fill = "gray66")+
      scale_x_discrete(expand = c(0.01, 0)) +
      scale_y_continuous(expand = c(0.01, 0)) +
      ggtitle(Title) +
      xlab(yLabel) + # the arguments for x and y labls are switched as the figure will be rotated
      ylab(xLabel) + # the arguments for x and y labls are switched as the figure will be rotated
      geom_hline(yintercept = 0) +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.text.x = element_text(size = xTxtSize, angle = 0, hjust = 0.5), # x and y not reversed as they are not associated with the roation of the axes.
            axis.text.y = element_text(size = yTxtSize, hjust = 0.5)) +
      coord_flip()

    if (errorbar == "SEM"){
      plt <- baseplt +
        geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = errorbarWidth,
                      position = position_dodge(0.9))
    } else if (errorbar == "SD") {
      plt <- baseplt +
        geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = errorbarWidth,
                      position = position_dodge(0.9))
    }

    ## add the right-side y axis
    grid.newpage()

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

    # export the file and draw a preview
    ggsave(filename = paste(objTitle,".vi.plot.pdf", sep = ""), plot = pltgtb,
           width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
    grid.draw(pltgtb) # preview
  }

  ## return the vi ranking and OOB err dataframes for the initial feature elimination
  outlst <- list(matrix_initial_FS = training_initFS,
                 feature_initial_FS = feature_initFS,
                 recur_vi_summary = outdfm_vi,
                 recur_OOB_err_summary = outdfm_OOB_err)

  sink(file = paste(objTitle,".initialFS.txt",sep = ""), append=FALSE) # dump the results to a file
  print(outlst)
  sink() # end dump

  return(assign(paste(objTitle, "_initial_FS", sep = ""), outlst, envir = .GlobalEnv)) # return a dataframe with the vi ranking dataframe
}
