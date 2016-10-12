.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.")
  return(TRUE)
}


#' @title rbioRF_vi
#'
#' @description Iterative random froest variable importance (vi) computation.
#' @param dfm Input data frame.
#' @param targetVar The target variable for random forest feature selection.
#' @param nTimes Number of iteration of random forest vi computation. Default is \code{50} times.
#' @param transpo If the dataframe needs to be transposed before random forest. Default is \code{TRUE}.
#' @param nTree Number of trees generated for each random forest iteration. Default is \code{1001} trees.
#' @return Outputs a \code{matrix} object with vi values for all the random forest iterations for each feature.
#' @importFrom randomForest randomForest importance
#' @examples
#' \dontrun{
#' rbioRF(dataframe, tgtvar_CntlvsStress, transpo = FLASE, nTree = 501)
#' }
#' @export
rbioRF_vi <- function(dfm,  targetVar, nTimes = 50, transpo = TRUE, nTree = 1001){

  ### load the dataframe/matrix
  if (transpo == TRUE){
    training <- t(dfm) # load the dataframe/matrix and transpose
  } else {
    training <- dfm
  }

  ### pepare the target variable
  tgt <- factor(as.character(targetVar), levels = unique(targetVar))

  ### prepare draw size. this uses down-sampling if the samples are unbalanced
  nlvl <- length(levels(tgt))
  size <- min(as.vector(table(tgt))) # down-sampling
  drawSize <- rep(size, nlvl)

  ### repeating random forest - iterative approach
  # pre-set an empty matrix with the number of columns same as the number of RF iterations
  # note that nrow is the number of features, hence the ncol of the traning set
  tmpMtx <- matrix(nrow = ncol(training), ncol = nTimes)

  tmpFunc <- function(n, m, mtx, tmpTraining, tmpTgt,
                      tmpTree, tmpSize){
    if (n == 0){
      rownames(mtx) <- colnames(tmpTraining)
      colnames(mtx) <- c(paste("accuracy", seq(m - 1), sep = "_"))
      return(mtx)
    } else {
      rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, importance = TRUE,
                         proximity = TRUE, drawSize = tmpSize)
      impt <- importance(rf, type = 1)
      mtx[, m] <- impt[, 1]
      tmpFunc(n - 1, m + 1, mtx, tmpTraining, tmpTgt, tmpTree, tmpSize)
    }
  }

  tmpFunc(n = nTimes, m = 1, mtx = tmpMtx, tmpTraining = training, tmpTgt = tgt,
          tmpTree = nTree, tmpSize = drawSize)

}


#' @title rbioRF_viplot
#'
#' @description Generate bargraph for vi values obtained from itarative random forest
#' @param fsdfm Input data frame.
#' @param n Number of features to show. Takes integer numbers. Default is \code{"all"}.
#' @param errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @return Outputs a \code{pdf} figure file.
#' @import ggplot2
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @examples
#' \dontrun{
#' rbioRF(dataframe, tgtvar_CntlvsStress, transpo = FLASE, nTree = 501)
#' }
#' @export
rbioRF_viplot <- function(fsdfm, n = "all",
                          errorbar = "SEM", errorbarWidth = 0.2,
                          xTxtSize = 10, yTxtSize =10,
                          plotWidth = 170, plotHeight = 150){

  ## prepare the dataframe
  fName <- rownames(fsdfm)
  fMean <- rowMeans(fsdfm)
  fSD <- apply(fsdfm, 1, sd)
  fSEM <- sapply(fSD, function(x)x/sqrt(ncol(fsdfm)))
  pltdfm <- data.frame(Targets = fName, Mean = fMean, SD = fSD, SEM = fSEM, stringsAsFactors = FALSE)
  pltdfm <- pltdfm[order(pltdfm$Mean), ]
  pltdfm$Targets <- factor(pltdfm$Targets, levels = unique(pltdfm$Targets))

  if (n != "all"){
    pltdfm <- tail(pltdfm, n)
  }

  ## boxplot
  # prepare plotting dataframe (draft only)
  loclEnv <- environment()
  baseplt <- ggplot(pltdfm, aes(x = Targets, y = Mean), environment = loclEnv) +
    geom_bar(position="dodge", stat="identity", color="black")+
    scale_x_discrete(expand = c(0.005, 0)) +
    ggtitle(NULL) +
    xlab(NULL) + # we can hide it using NULL
    ylab("Mean Decrease in Accuracy") +
    geom_hline(yintercept = 0) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_text(size = xTxtSize, angle = 90, hjust = 1),
          axis.text.y = element_text(size = yTxtSize, hjust = 0.5)) +
    coord_flip()

  if (errorbar == "SEM"){
    plt <- baseplt +
      geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = errorbarWidth,
                    position = position_dodge(0.9)) +
      scale_y_continuous(limits = c(with(pltdfm, min(Mean - SEM) * 1.1), with(pltdfm, max(Mean + SEM) * 1.1)))
  } else if (errorbar == "SD") {
    plt <- baseplt +
      geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = errorbarWidth,
                    position = position_dodge(0.9)) +
      scale_y_continuous(limits = c(with(pltdfm, min(Mean - SD) * 1.1), with(pltdfm, max(Mean + SD) * 1.1)))
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
  ggsave(filename = paste(deparse(substitute(fsdfm)),".plot.pdf", sep = ""), plot = pltgtb,
         width = plotWidth, height = plotHeight, units = "mm",dpi = 600) # deparse(substitute(dfm)) converts object name into a character string
  grid.draw(pltgtb) # preview

}
