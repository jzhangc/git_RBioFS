% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbioClass_plsda.R
\name{rbioClass_plsda_roc_auc}
\alias{rbioClass_plsda_roc_auc}
\title{rbioClass_plsda_roc_auc()}
\usage{
rbioClass_plsda_roc_auc(
  object,
  newdata,
  newdata.label,
  center.newdata = TRUE,
  rocplot = TRUE,
  plot.comps = 1:object$ncomp,
  plot.smooth = FALSE,
  multi_plot.ncol = length(plot.comps),
  multi_plot.nrow = 1,
  multi_plot.legend.pos = c("bottom", "top", "left", "right"),
  plot.rightsideY = TRUE,
  plot.SymbolSize = 2,
  plot.lineSize = 1,
  plot.display.Title = TRUE,
  plot.titleSize = 10,
  plot.fontType = "sans",
  plot.xLabel = "1 - specificity",
  plot.xLabelSize = 10,
  plot.xTickLblSize = 10,
  plot.yLabel = "sensitivity",
  plot.yLabelSize = 10,
  plot.yTickLblSize = 10,
  plot.legendSize = 9,
  plot.Width = 170,
  plot.Height = 150,
  verbose = TRUE
)
}
\arguments{
\item{object}{A \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.}

\item{newdata}{Newdata (test data) for ROC-AUC analysis, excluding labels, i.e. y. If missing, the function will use the transformed data from the model object.}

\item{newdata.label}{Newdata label vector (i.e. test data y). If missing, the function will use the training data and its corresponding labels.}

\item{center.newdata}{Only set when both \code{newdata} and \code{newdata.label} are set, if to center.scale newdata. Default is \code{TRUE}.}

\item{rocplot}{If to generate a ROC plot. Default is \code{TRUE}.}

\item{plot.comps}{Number of comps to plot. Default is \code{1:object$ncomp}}

\item{plot.smooth}{If to smooth the curves. Uses binomial method to smooth the curves. Default is \code{FALSE}.}

\item{multi_plot.ncol}{Number of columns on one figure page. Default is \code{length(plot.comps)}.}

\item{multi_plot.nrow}{Number of rows on one figure page. Default is \code{1}.}

\item{multi_plot.legend.pos}{The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.}

\item{plot.rightsideY}{If to show the right side y-axis. Default is \code{FALSE}. Note: doesn't seem to be necessary as PLS-DA always has at least two y classes.}

\item{plot.SymbolSize}{Symbol size. Default is \code{2}.}

\item{plot.lineSize}{Line size. Default is \code{1}.}

\item{plot.display.Title}{If to show the name of the y class. Default is \code{TRUE}.}

\item{plot.titleSize}{The font size of the plot title. Default is \code{10}.}

\item{plot.fontType}{The type of font in the figure. Default is "sans". For all options please refer to R font table, which is available on the website: \url{http://kenstoreylab.com/?page_id=2448}.}

\item{plot.xLabel}{X-axis label. Type with quotation marks. Could be NULL. Default is \code{"1 - specificity"}.}

\item{plot.xLabelSize}{X-axis label size. Default is \code{10}.}

\item{plot.xTickLblSize}{X-axis tick label size. Default is \code{10}.}

\item{plot.yLabel}{Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"sensitivity"}.}

\item{plot.yLabelSize}{Y-axis label size. Default is \code{10}.}

\item{plot.yTickLblSize}{Y-axis tick label size. Default is \code{10}.}

\item{plot.legendSize}{Legend size. Default is \code{9}.}

\item{plot.Width}{ROC width. Default is \code{170}.}

\item{plot.Height}{ROC height. Default is \code{150}.}

\item{verbose}{whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.}
}
\value{
Prints AUC values in the console. And a pdf file for ROC plot
}
\description{
ROC-AUC analysis and ploting for plsda model
}
\details{
Uses pROC module to calculate ROC.
}
\examples{
\dontrun{
rbioClass_plsda_roc_auc(object = model_binary, rocplot = TRUE, plot.comps = 1:2)
}
}
