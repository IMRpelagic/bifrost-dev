#' Set up parameter list for maturity
#'
#' @param parameter data.frame, initial values for p1,p2,p3 and nu by year
#' @param agegr agegroup, should either "2-3" or "3-4"
#' @param method which method to be used for initiating parameters. Either "user-provided",
#' "start-at-vec" or "grid-search".
#' @param data.list data list, to be fed into maturity estimation
#' @param init numeric vector of initial values. Used for method = "user-provided".
#' @param year integer, year
#'
#' @return list, parameter object for running maturity estimation by TMB
#' @export
#'
#' @examples
#' data(maturityInitialParameters)
#' createMaturityParameters(parameter=maturityInitialParameters, year = 2017)
createMaturityParameters <- function(parameter=NULL,
                                     year=NULL,
                                     agegr = "3-4",
                                     method = "user-provided",
                                     data.list = NULL,
                                     init = numeric(4)
                                     ){
  if(method == "user-provided"){
    if(!is.data.frame(parameter)) stop("parameter is not a data.frame.")
    if(!(is.numeric(year)|is.integer(year))) stop("year is not an integer.")
    if(!all(c("p1","p2","p3","nu","year") %in% names(parameter)))
      stop("parameter data frame does not contain necessary columns:\np1, p2, p3, nu and year")
    if(!(agegr == "3-4"|agegr == "2-3"))
      stop("wrong age group.")
    i <- which(parameter$year == year & parameter$agegr == agegr)
    return(createMaturityParameters(
      method = "start-at-vec",
      init =as.numeric(log(parameter[i, c("p1","p2","p3","nu")]))))
  }else if(method =="start-at-vec"){
    return(list(
      lnp1 = init[1],
      lnp2 = init[2],
      lnp3 = init[3],
      lnnu = init[4]
    ))
  } else if(method == "grid-search"){
    par.list <- createMaturityParameters(method = "start-at-vec")
    obj <- TMB::MakeADFun(data = c(model = "mature", data.list),
                          parameters = par.list,
                          DLL = "bifrost_TMBExports")

    dummy1 <- seq(0.01,.99,length.out = 30)
    dummy1 <- log(dummy1)-log(1-dummy1)
    dummy2 <- log(seq(12,16,length.out = 30))
    dummy3 <- seq(-1, 30,length.out = 20)
    X<-expand.grid(dummy1, dummy2, dummy1, dummy3)
    return(createMaturityParameters(init =as.numeric(X[which.min(apply(X,1,obj$fn)),]),
                                                     method = "start-at-vec"))
  }
}

