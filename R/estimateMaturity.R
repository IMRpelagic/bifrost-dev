#' Function for estimating maturity
#'
#' @param data list of data inputs
#' @param parameters list of parameters
#' @param ... additional arguments for TMB::MakeADFun
#'
#' @return list of TMB object and optimization object
#' @export
#'
#' @examples
#' \dontrun{
#' data(cap)
#' data(catch)
#' data(maturityInitialParameters)
#'
#' #..set up data list..
#' data.list <- createMaturityData(cap,
#'                                 catch,
#'                                 min_age = 2,
#'                                 max_age = 3,
#'                                 start_year =1972,
#'                                 end_year = 2010)
#' # ..set up parameter list..
#' par.list <- createMaturityParameters(parameter = maturityInitialParameters,
#'                                     year = data.list$start_year, agegr = "2-3")
#' mFit <- estimateMaturity(data = data.list, parameters = par.list, silent =TRUE)
#' }
estimateMaturity <- function(data, parameters, ...){
  NumbersAtLength <- data$Nl
  data$Nl <- data$Nl[, data$min_age:data$max_age]
  # .. initialize TMB object ..
  obj <- TMB::MakeADFun(data = c(model = "mature", data),
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
  return.list$NumbersAtLength = NumbersAtLength
  attr(return.list, "class") <- c("maturity","list")
  return(return.list)
}

#' Print maturity object
#'
#' @param x maturity
#' @param ... additional arguments
#'
#' @return printout
#' @export
#'
print.maturity <- function(x,...){
  cat("Convergence? ", x$opt$convergence,": ", x$opt$message,
      "\n-------------",
      "\nEstimates? \n")
  print(x$sumsdrep)
  cat("-------------\n")

}

#' Summary of estimation of Maturity
#'
#' @param obj maturity, object from running estimation
#' @param ... additional arguments
#'
#' @return summary
#' @export
#'
#'
#'
summary.maturity <- function(obj, ...){
  summary.list <- list()
  attr(summary.list, "class") <- c("summary.maturity","list")
  tab <- obj$sumsdrep[-(1:4),]
  #tab[,1:2] <- obj$sumsdrep[-(1:4),]
  #tab[,3] <- tab[,1]/tab[,2]
  #tab[,4] <- 2 * pnorm(abs(tab[,3]), lower.tail = FALSE)
  colnames(tab) <- c("Estimate", "Std. Error", "Test score", "p-value*")
  rownames(tab) <- rownames(obj$sumsdrep[-(1:4),])
  summary.list$result.tab <- tab
  summary.list$convergence.code <- obj$opt$convergence
  summary.list$convergence.message <- obj$opt$message
  summary.list$loglikelihoodvalue <- obj$opt$objective
  summary.list$aic <- 2 * (obj$opt$objective +length(obj$opt$par))
  summary.list$maturitytable <-tibble::tibble(
    ml = obj$data$meanlength,
    r = maturing(meanlength = ml,
                 p1 = tab["p1",1], p2 = tab["p2",1]))
  summary.list$years <- obj$data$start_year:obj$data$end_year
  summary.list$ages <-  obj$data$min_age:obj$data$max_age
  return(summary.list)
}

#' Print summary of maturity
#'
#' @param x object of type summary.maturity
#' @param ... addition arguments
#'
#' @return printout
#' @export
#'
#'
print.summary.maturity <- function(x,...) {
  print(x$result.tab)
  cat("\n* Using Gaussian approximation for p-values.\n")
  cat("\n-------------------------------------------",
      "\nConvergence code:             ", x$convergence.code,
      "\nCovergence message:           ", x$convergence.message,
      "\nNegative loglikelihood value: ", x$loglikelihoodvalue,
      "\nAkaike Information Criteria:  ", x$aic,
      "\n-------------------------------------------\n"
  )
  # ..plot maturity..
  print(plot(x))
}
