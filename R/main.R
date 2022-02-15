# Script for å kjøre "ren" TMB-kode
# ---- 
#
#

library(TMB)
compile("src/mature.cpp")
dyn.load(dynlib("src/mature"))

# Data som inngår i modningsmodellen
load("data/cap.rda") # <-- Data fra tokt
head(cap)
# Fangst per alder, år og måned (kolonner oct-dec, jan-sep)
load("data/catch.rda")
head(catch)

# Initial verdier per år (satt sammen av Mahmood)
load("data/maturityInitialParameters.rda")
head(maturityInitialParameters)

source("R/createMaturityData.R")
source("R/createMaturityParameter.R")

year <- 2010
#.. Create data list: ..
data.list <- createMaturityData(cap,
                                catch,
                                min_age = 2, # estimerer modning fra alder 2-3
                                max_age = 3,
                                start_year =1972,
                                end_year = 2010)
par.list <- createMaturityParameters(parameter = maturityInitialParameters,
                                     year = data.list$end_year, agegr = "2-3")

NumbersAtLength <- data.list$Nl
data.list$Nl <- data.list$Nl[, data.list$min_age:data.list$max_age]
# .. initialize TMB object ..
obj <- TMB::MakeADFun(data = data.list,
                      parameters = par.list,
                      DLL = "mature", silent =TRUE)
# .. Optimize likelihood ..
opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
opt

# -----------------
# -- CONSUMPTION --
# -----------------
load("data/consumptionData.rda")

compile("src/consumption.cpp")
dyn.load(dynlib("src/consumption"))

par.list <- list(
  logCmax   = log(1.2),
  logChalf  = log(1e2),
  logalpha  = log(2),
  logbeta   = log(2),
  logSigma  = log(1e3)
)
obj <- TMB::MakeADFun(data = consumptionData,
                      parameters = par.list,
                      DLL = "consumption", silent =TRUE)
# .. Optimize likelihood ..
opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
opt
