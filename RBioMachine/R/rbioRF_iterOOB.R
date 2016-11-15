#' @title rbioRF_iterOOB
#'
#' @description Iterative nested random froest variable importance (vi) and OOB error rate computation. (to be completed)
#' @param x Input dataframe or matrix. Make sure to arrange the data with features as column names.
#' @param targetVar The target variable for random forest feature selection. This is a factor object.
#' @param nTimes Number of iteration of random forest vi computation. Default is \code{50} times.
#' @param nTree Number of trees generated for each random forest iteration. Default is \code{1001} trees.
#' @param multicore If to use parallel computing. Default is \code{TRUE}.
#' @param plot If to plot a bargraph to visualize vi and the ranking. Default is \code{TRUE}
#' @param Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param n Number of features to show. Takes integer numbers. Default is \code{"all"} (make sure to include quotation marks).
#' @param errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param symbolSize Size of the symbol. Default is \code{2}.
#' @param xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @return Outputs a \code{list} object with  OOB error rate summary, and a joint-point curve in \code{csv} format.
#' @details Make sure to arrange data (dfm) with feature (e.g., gene) as variables (i.e., columns), and rownames as sample names.
#' @import ggplot2
#' @importFrom randomForest randomForest importance
#' @importFrom parallel detectCores makeCluster stopCluster parApply parLapply
#' @examples
#' \dontrun{
#' rbioRF_iterOOB(training_HCvTC, tgtVar_HCvTC, multicore = TRUE)
#' }
#' @export
rbioRF_iterOOB <- function(x, targetVar, nTimes = 50, nTree = 1001,
                           multicore = TRUE,
                           plot = TRUE, n = "all",
                           Title = NULL, xLabel = NULL, yLabel = NULL,
                           errorbar = "SEM", errorbarWidth = 0.2,
                           symbolSize = 2, xTxtSize = 10, yTxtSize =10,
                           plotWidth = 170, plotHeight = 150){

  ## mark time
  start <- Sys.time()


  ## prepare the dataframe
  training <- data.frame(x)

  ### pepare the target variable
  tgt <- factor(as.character(targetVar), levels = unique(targetVar))

  ### prepare draw size. this uses down-sampling if the samples are unbalanced
  nlvl <- length(levels(tgt))
  size <- min(as.vector(table(tgt))) # down-sampling
  drawSize <- rep(size, nlvl)


  ## prepare blank tree OOB error matrics
  singleerrmtx <- matrix(nrow = 1, ncol = nTimes) # for the iterative OOB error rates from a single tree
  ooberrmtx <- matrix(nrow = ncol(training), ncol = nTimes) # for the iterative OOB error rates from all trees.


  if (!multicore){

    ## signle core computing: recursive structure
    tmpFunc <- function(n, m, tmptimes, tmperrmtx, tmpTraining, tmpTgt,
                        tmpTree, tmpSize){


      if (n == 0){
        return(tmperrmtx)

      } else {
        if (ncol(tmpTraining) < 4){
          rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, importance = TRUE,
                             proximity = TRUE, drawSize = tmpSize)
        } else {
          rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, mtry = max(ceiling(ncol(tmpTraining) / 3), 2),
                             importance = TRUE,
                             proximity = TRUE, drawSize = tmpSize)
        }

        tmperrmtx[, m] <- rf$err.rate[tmptimes, 1] # fill the OOB error rate
        tmpFunc(n - 1, m + 1, tmptimes, tmperrmtx, tmpTraining, tmpTgt,
                tmpTree, tmpSize)
      }
    }

    tmpFunc2 <- function(i, j, tmp2mtx){
      if (i == 0){
        rownames(tmp2mtx) <- seq(j - 1)
        colnames(tmp2mtx) <- c(paste("OOB_error_tree_rep", seq(nTimes), sep = "_"))

        return(tmp2mtx)
      } else {
        tmp2mtx[j, ] <- tmpFunc(n =nTimes, m = 1, tmptimes = nTree, tmperrmtx = singleerrmtx,
                                tmpTraining = training[1:j], tmpTgt = tgt, tmpTree = nTree,
                                tmpSize = drawSize)
        tmpFunc2(i - 1, j + 1, tmp2mtx)
      }
    }

    mtxforfunc2 <- ooberrmtx

    ooberrmtx <- tmpFunc2(i = ncol(training), j = 1, tmp2mtx = mtxforfunc2) # j is the tree index

  } else {

    ## parallel computing
    # set up cpu cluster
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    tmpfunc4 <- function(j, ...){

      n_cores2 <- parallel::detectCores() - 1
      cl2 <- parallel::makeCluster(n_cores2)
      on.exit(parallel::stopCluster(cl2)) # close connect when exiting the function


      errmtx <- matrix(nrow = 1, ncol = nTimes) # for the iterative OOB error rates from a single tree

      # iterative RF using par-apply functions
      tmpfunc3 <- function(i, ...){
        rf <- randomForest::randomForest(x = training[1:j], y = tgt, ntree = nTree, importance = TRUE,
                                         proximity = TRUE, drawSize = drawSize)

        tmperrmtx <- rf$err.rate[nTree, 1] # compute the OOB error rate
        lst <- list(tmperrmtx = tmperrmtx)
      }

      tmp <- parallel::parLapply(cl2, X = 1:nTimes, fun = tmpfunc3)

      for (i in 1:nTimes){
        errmtx[, i] <- tmp[[i]]$tmperrmtx # fill the OOB error rate
      }

      list = list(errmtx = errmtx)
    }

    l <- parLapply(cl, X = 1:ncol(training), tmpfunc4)

    for (p in 1:ncol(training)){
      ooberrmtx[p, ] <- l[[p]]$errmtx
    }

    rownames(ooberrmtx) <- seq(ncol(training))
    colnames(ooberrmtx) <- c(paste("OOB_error_tree_rep", seq(nTimes), sep = "_"))

  }

  ## perpare the summary dataframe for OOB error rates
  ooberrnames <- rownames(ooberrmtx)
  ooberrmean <- rowMeans(ooberrmtx)
  ooberrSD <- apply(ooberrmtx, 1, sd)
  ooberrSEM <- sapply(ooberrSD, function(x)x / sqrt(ncol(ooberrmtx)))
  ooberrsummary <- data.frame(Features = ooberrnames, Mean = ooberrmean, SD = ooberrSD,
                              SEM = ooberrSEM, stringsAsFactors = FALSE)
  ooberrsummary$Features <- factor(ooberrsummary$Features, levels = unique(ooberrsummary$Features))


  ## plot
  if (plot){

    loclEnv <- environment()

    # prepare plotting dataframe
    if (n != "all"){
      pltdfm <- head(ooberrsummary, n)
    } else {
      pltdfm <- ooberrsummary
    }

    # plotting
    baseplt <- ggplot(ooberrsummary, aes(x = Features, y = Mean, group = 1), environment = loclEnv) +
      geom_line() +
      geom_point(size = symbolSize) +
      scale_x_discrete(expand = c(0, 0)) +
      ggtitle(Title) +
      xlab(xLabel) + # the arguments for x and y labls are switched as the figure is rotated
      ylab(yLabel) + # the arguments for x and y labls are switched as the figure is rotated
      geom_hline(yintercept = 0) +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.text.x = element_text(size = xTxtSize),
            axis.text.y = element_text(size = yTxtSize, hjust = 0.5))

    if (errorbar == "SEM"){
      plt <- baseplt +
        geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = errorbarWidth) +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(0, with(ooberrsummary, max(Mean + SEM) * 1.2)))
    } else if (errorbar == "SD") {
      plt <- baseplt +
        geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = errorbarWidth) +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(0, with(ooberrsummary, max(Mean + SD) * 1.2)))
    }

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
    ggsave(filename = paste(deparse(substitute(x)),".OOBplot.pdf", sep = ""), plot = pltgtb,
           width = plotWidth, height = plotHeight, units = "mm",dpi = 600) # deparse(substitute(x)) converts object name into a character string
    grid.draw(pltgtb) # preview


  }


  ## mark time
  end <- Sys.time()

  ## output
  outlst <- list(OOB_error_rate_summary = ooberrsummary,
                 runtime = end - start)

  return(outlst)

}
