#' @title rbioRF_iterOOB
#'
#' @description Iterative nested random froest variable importance (vi) and OOB error rate computation. (to be completed)
#' @param x Input dataframe or matrix. Make sure to arrange the data with features as column names.
#' @param targetVar The target variable for random forest feature selection. This is a factor object.
#' @param nTimes Number of iteration of random forest vi computation. Default is \code{50} times.
#' @param nTree Number of trees generated for each random forest iteration. Default is \code{1001} trees.
#' @param multicore If to use parallel computing. Default is \code{TRUE}.
#' @return Outputs a \code{matrix} object with  OOB error rates.
#' @details Make sure to arrange data (dfm) with feature (e.g., gene) as variables (i.e., columns), and rownames as sample names.
#' @importFrom randomForest randomForest importance
#' @importFrom parallel detectCores makeCluster stopCluster parApply parLapply
#' @examples
#' \dontrun{
#' rbioRF_iterOOB(training_HCvTC, tgtVar_HCvTC, multicore = TRUE)
#' }
#' @export
rbioRF_iterOOB <- function(x, targetVar, nTimes = 50, nTree = 1001,
                   multicore = TRUE){


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
        rownames(tmp2mtx) <- c(paste("tree", seq(j - 1), sep = "_"))
        colnames(tmp2mtx) <- c(paste("OOB_error_tree_rep", seq(nTimes), sep = "_"))

        return(tmp2mtx)
      } else {
        tmp2mtx[j, ] <- tmpFunc(n =nTimes, m = 1, tmptimes = nTree, tmperrmtx = singleerrmtx,
                                tmpTraining = training[1:j], tmpTgt = tgt, tmpTree = nTree,
                                tmpSize = drawSize)
        tmpFunc2(i - 1, j + 1, tmp2mtx)
      }
    }

    tmpFunc2(i = ncol(training), j = 1, tmp2mtx = ooberrmtx) # j is the tree index

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

        tmperrmtx <- rf$err.rate[nTree, 1] # fill the OOB error rate
        lst <- list(tmperrmtx = tmperrmtx)
      }

      tmp <- parallel::parLapply(cl2, X = 1:nTimes, fun = tmpfunc3)

      for (i in 1:nTimes){
        errmtx[, i] <- tmp[[i]]$tmperrmtx
      }

      list = list(errmtx = errmtx)
    }

    l <- parLapply(cl, X = 1:ncol(training), tmpfunc4)

    for (p in 1:ncol(training)){
      ooberrmtx[p, ] <- l[[p]]$errmtx
    }

    rownames(ooberrmtx) <- c(paste("tree", seq(ncol(training)), sep = "_"))
    colnames(ooberrmtx) <- c(paste("OOB_error_tree_rep", seq(nTimes), sep = "_"))

    return(ooberrmtx)
  }
}
