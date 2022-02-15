#' Function for estimating consumption model
#'
#' @param data list of data inputs
#' @param parameters list of parameters
#' @param ... additional arguments for TMB::MakeADFun
#'
#' @return list of TMB object and optimization object
#' @export
#'
#' @examples \dontrun{
#' data("consumptionData")
#' par.list <- list(
#'   logCmax   = log(1.2),
#'   logChalf  = log(1e2),
#'   logalpha  = log(2),
#'   logbeta   = log(2),
#'   logSigma  = log(1e3)
#' )
#' cFit <- estimateConsumption(data = consumptionData, parameters = par.list, silent =TRUE,
#'                             map  = list(logalpha = factor(NA),
#'                                         logbeta = factor(NA)))
#' }
estimateConsumption <- function(data, parameters, ...){
  # .. initialize TMB object ..
  obj <- TMB::MakeADFun(data = c(model = "consumption", data),
                        parameters = parameters,
                        DLL = "bifrost_TMBExports", ...)
  # .. Optimize likelihood ..
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr)

  # ..Return results..

  return.list <- list()
  return.list$obj <- obj
  return.list$opt <- opt
  return.list$data <- data
  return.list$parameters <- parameters
  return.list$sdrep <- TMB::sdreport(obj)
  return.list$sumsdrep <- TMB::summary.sdreport(return.list$sdrep, p.value = TRUE)
  attr(return.list, "class") <- c("consumption","list")
  return(return.list)
}

#' Print consumption object
#'
#' @param x consumption object
#' @param ... additional arguments
#'
#' @return printout
#' @export
#'
print.consumption <- function(x,...){
  cat("-- CONSUMPTION MODEL: --",
    "\nConvergence? ", x$opt$convergence,": ", x$opt$message,
      "\n-------------",
      "\nEstimates? \n")
  par.names <- c("Cmax", "Chalf", "alpha", "beta", "sigma")
  print(x$sumsdrep[par.names, ])
  cat("-------------\n")

}

#' Summary of estimation of Consumption
#'
#' @param obj consumption, object from running estimation
#' @param ... additional arguments
#'
#' @return summary
#' @export
#'
#'
#'
summary.consumption <- function(obj, ...){
  summary.list <- list()
  class(summary.list)<-"summary.consumption"
  par.names <- c("Cmax", "Chalf", "alpha", "beta", "sigma")
  tab <-  obj$sumsdrep[par.names,]
  colnames(tab) <- c("Estimate", "Std. Error", "Zscore", "p-value*")
  rownames(tab) <- par.names
  summary.list$result.tab <- tab
  summary.list$convergence.code <- obj$opt$convergence
  summary.list$convergence.message <- obj$opt$message
  summary.list$loglikelihoodvalue <- obj$opt$objective
  summary.list$aic <- 2 * (obj$opt$objective +length(obj$opt$par))
 # summary.list$maturitytable <-tibble::tibble(
#    ml = x$data$meanlength,
#    r = maturing(meanlength = ml,
#                 p1 = tab["p1",1], p2 = tab["p2",1]))
 # summary.list$years <- x$data$start_year:x$data$end_year
 # summary.list$ages <-  x$data$min_age:x$data$max_age
  return(summary.list)
}

#' Print summary of consumption
#'
#' @param x object of type summary.consumption
#' @param ... addition arguments
#'
#' @return printout
#' @export
#'
#'
print.summary.consumption <- function(x,...) {
  print(x$result.tab)
  cat("\n* Using Gaussian approximation for p-values.\n")
  cat("\n-------------------------------------------",
      "\nConvergence code:             ", x$convergence.code,
      "\nCovergence message:           ", x$convergence.message,
      "\nNegative loglikelihood value: ", x$loglikelihoodvalue,
      "\nAkaike Information Criteria:  ", x$aic,
      "\n-------------------------------------------\n"
  )
  # ..plot consumption..
  #print(plotMaturity(x))
}

#' Plot consumption
#'
#' @param x object of type consumption
#' @param ... additional argument
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#'  data("consumptionData")
#' par.list <- list(
#'   logCmax   = log(1.2),
#'   logChalf  = log(1e2),
#'   logalpha  = log(2),
#'   logbeta   = log(2),
#'   logSigma  = log(1e3)
#' )
#' cFit <- estimateConsumption(data = consumptionData, parameters = par.list, silent =TRUE,
#'                             map  = list(logalpha = factor(NA),
#'                                         logbeta = factor(NA)))
#' plot(cFit)
#' }
plot.consumption <- function(x, ...){
  tmp <- tibble::tibble(
    "Years" = as.numeric(rownames(x$data$Nco)),
    "EmpiricalConsumption" = x$obj$report()$Econscolsum,
    "ModelConsumption" = x$sumsdrep[rownames(x$sumsdrep) == "conscolsum",1],
    "SE" = x$sumsdrep[rownames(x$sumsdrep) == "conscolsum",2],
    "low" =pmax(0, ModelConsumption - 2*SE),
    "high" = ModelConsumption + 2*SE
  )
  print(
  ggplot2::ggplot(tmp, ggplot2::aes_string(x = "Years",
                                    y = "ModelConsumption",
                                    ymin = "low",
                                    ymax = "high")) +
    ggplot2::geom_ribbon(alpha = .1, fill = "blue") +
    ggplot2::geom_line(lwd = 2, col = "blue") +
    ggplot2::geom_line(ggplot2::aes_string(y = "EmpiricalConsumption"))+
    ggplot2::geom_point(ggplot2::aes_string(y = "EmpiricalConsumption")) +
    ggplot2::scale_x_continuous(name = "Year", seq(1800,2200,2))+
    ggplot2::scale_y_continuous(name = "Total consumption of capelin by cod")
  )

}
