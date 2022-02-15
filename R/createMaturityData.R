#' Title
#'
#' @param cap data.frame, cap table. Should contain the following column names:
#' meanlength(cm), 1,2,3,4,5 and year, where 1-4 are ages and 5 is a 5+ group.
#' @param catch data.frame, catches data. Should contain 14 columns including
#' age, year and the remaining 12 for each month.
#' @param min_age integer. Minimum age
#' @param max_age integer. Maximum age
#' @param start_year integer. Start year
#' @param end_year integer. End year
#'
#' @return list, data argument for TMB
#' @export
#'
#' @examples
#' data(cap)
#' data(catch)
#' data <- createMaturityData(
#' cap = cap,
#' catch = catch,
#' min_age = 3,
#' max_age = 4,
#' start_year = 2017,
#' end_year = 2018)
#' str(data)
#'
#'
createMaturityData <- function(cap,
                               catch,
                               min_age = 3,
                               max_age = 4,
                               start_year,
                               end_year)
{
  if(!all(c("meanlength(cm)",1:5, "year") %in% names(cap)))
    stop("cap does have the correct column names: \nmeanlength(cm), 1,2,3,4,5 and year")
  meanlength = cap$`meanlength(cm)`[cap$year == start_year]
  meanweight = cap$`meanweight(g)`[cap$year == end_year]
  start_index=start_year - min(cap$year) + 1
  end_index=end_year - min(cap$year) + 1

  # Nl is the reported number of capelin, in age classes min_age to max_age, start_year to end_year
  # It has the following strucutre: Nl(year,l,a)=Nl((year-1)*length_l+l,a)
  # dim(Nl): length(meanlength)*(end_year-start_year) x (max_age-min_age)
  Nl <- cap[cap$year %in% start_year:end_year, names(cap) %in% 1:5]
  # N is the total number of reported capelin for one age class and one year:
  # N(year,a)=sum_l Nl(year,l,a)
  N=colSums(Nl[1:length(meanlength),min_age:max_age])
  for(i in 2:(1+(end_index-start_index))){        #(start_index+1):end_index
    tmp=colSums(Nl[((i-1)*length(meanlength)+1):(i*length(meanlength)),min_age:max_age])
    N=rbind(N,tmp)
  }
  rownames(N) <- start_year:end_year

  # catchM are the monthly catches in years start_year to end_year; for all age classes and months october-september
  # It has the following strucutre: C(year,a,m)=C((year-1)*(min_age-max_age)+(a-min_age),m)
  C = subset(catch,
             subset = catch$year %in% start_year:end_year & catch$age %in% min_age:max_age,
             select = -c(age,year))

    # save dimensions, as well
  length_l=length(meanlength)
  index_Nl=start_year+dim(Nl)[1]-1
  index_C=start_year+dim(C)[1]-1
  data <- list(
    min_age = min_age,
    max_age = max_age,
    start_year = start_year,
    end_year = end_year,
    length_l = length_l,
    meanlength = meanlength,
    meanweight = meanweight,
    Nl = as.matrix(Nl),
    N = as.matrix(N),
    C = as.matrix(C)
  )
  return(data)
}

#' Update maturity data with new Nl object
#'
#' @param data list, maturity data to be updated
#' @param Nl matrix, numbers at length group
#'
#' @return updated maturity data list
#' @export
#'
updateMaturityData <- function(data, Nl){
  start_index = 1
  end_index = data$end_year - data$start_year + 1
  N=colSums(Nl[1:length(data$meanlength),data$min_age:data$max_age])
  for(i in 2:(1+(end_index-start_index))){        #(start_index+1):end_index
    tmp=colSums(Nl[((i-1)*length(data$meanlength)+1):(i*length(data$meanlength)),data$min_age:data$max_age])
    N=rbind(N,tmp)
  }
  data$Nl = Nl
  data$N = N
  return(data)
}
