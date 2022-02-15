Bifrost
================

## Bifrost uten R pakke-strukturen

Her er en enkel sammensetning av sentrale skript fra R-pakken bifrost. I
src/ folderen finner du TMB-koden (.cpp-filer). I data folderen ligger
dataene som .rda filer som kan lastes inn ved, f.eks.

``` r
load("data/cap.rda")
head(cap)
```

    ##   length.group     1 2 3 4 5 sum(10e9) biomass(10e3t) meanweight(g)
    ## 1            1  0.00 0 0 0 0      0.00            0.0           0.0
    ## 2            2  0.00 0 0 0 0      0.00            0.0           0.0
    ## 3            3  2.63 0 0 0 0      2.63            1.8           0.7
    ## 4            4  5.12 0 0 0 0      5.12            5.1           1.0
    ## 5            5 11.89 0 0 0 0     11.89           14.3           1.2
    ## 6            6 26.62 0 0 0 0     26.62           45.3           1.7
    ##   meanlength(cm) year
    ## 1           5.25 1972
    ## 2           5.75 1972
    ## 3           6.25 1972
    ## 4           6.75 1972
    ## 5           7.25 1972
    ## 6           7.75 1972

I R-mappen har jeg laget et skript *main.R* som kan kjøres som et
fungerende eksempel. Kjører det under her:

### R/Main.R

``` r
# Script for å kjøre "ren" TMB-kode
# - bruker noen funksjoner fra Rbifrost (fra R-mappen)
# - til å sette opp data og initial parametre

library(TMB)
compile("src/mature.cpp")
```

    ## [1] 0

``` r
dyn.load(dynlib("src/mature"))

# Data som inngår i modningsmodellen
load("data/cap.rda") # <-- Data fra tokt
head(cap)
```

    ##   length.group     1 2 3 4 5 sum(10e9) biomass(10e3t) meanweight(g)
    ## 1            1  0.00 0 0 0 0      0.00            0.0           0.0
    ## 2            2  0.00 0 0 0 0      0.00            0.0           0.0
    ## 3            3  2.63 0 0 0 0      2.63            1.8           0.7
    ## 4            4  5.12 0 0 0 0      5.12            5.1           1.0
    ## 5            5 11.89 0 0 0 0     11.89           14.3           1.2
    ## 6            6 26.62 0 0 0 0     26.62           45.3           1.7
    ##   meanlength(cm) year
    ## 1           5.25 1972
    ## 2           5.75 1972
    ## 3           6.25 1972
    ## 4           6.75 1972
    ## 5           7.25 1972
    ## 6           7.75 1972

``` r
# Fangst per alder, år og måned (kolonner oct-dec, jan-sep)
load("data/catch.rda")
head(catch)
```

    ##   age year        oct        nov        dec        jan      feb        mar apr
    ## 1   1 1972 0.00040000 0.00024000 0.00016000 0.00000000 0.000000 0.00000000   0
    ## 2   2 1972 0.22621739 0.13573043 0.09048696 0.00520000 0.026000 0.02080000   0
    ## 3   3 1972 0.16786677 0.10072006 0.06714671 0.25672112 1.283606 1.02688449   0
    ## 4   4 1972 0.03472727 0.02083636 0.01389091 0.02366479 0.118324 0.09465918   0
    ## 5   5 1972 0.00615000 0.00369000 0.00246000 0.00000000 0.000000 0.00000000   0
    ## 6   1 1973 0.00000000 0.00000000 0.00000000 0.00000000 0.000000 0.00000000   0
    ##   may jun jul     aug     sep
    ## 1   0   0   0 0.00000 0.00000
    ## 2   0   0   0 1.02296 0.25574
    ## 3   0   0   0 7.41944 1.85486
    ## 4   0   0   0 3.31168 0.82792
    ## 5   0   0   0 0.25320 0.06330
    ## 6   0   0   0 0.01000 0.00250

``` r
# Initial verdier per år (satt sammen av Mahmood)
load("data/maturityInitialParameters.rda")
head(maturityInitialParameters)
```

    ##   year agegr  p1   p2   p3 nu unknown
    ## 1 1972   2-3 0.3 14.0 0.02 10       0
    ## 2 1973   2-3 0.6 13.5 0.05 10       0
    ## 3 1974   2-3 0.4 12.5 0.02  2       0
    ## 4 1975   2-3 0.5 12.5 0.02  1       0
    ## 5 1976   2-3 0.3 15.0 0.02  6       0
    ## 6 1977   2-3 0.6 13.0 0.02  2       0

``` r
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
```

    ## Warning in stats::nlminb(obj$par, obj$fn, obj$gr): NA/NaN function evaluation

``` r
opt
```

    ## $par
    ##      lnp1      lnp2      lnp3      lnnu 
    ## -2.110836  2.741552 -2.713003  0.849294 
    ## 
    ## $objective
    ## [1] 148.5777
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $iterations
    ## [1] 37
    ## 
    ## $evaluations
    ## function gradient 
    ##       44       38 
    ## 
    ## $message
    ## [1] "relative convergence (4)"

``` r
# -----------------
# -- CONSUMPTION --
# -----------------
load("data/consumptionData.rda")

compile("src/consumption.cpp")
```

    ## [1] 0

``` r
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
```

    ## $par
    ##    logCmax   logChalf   logalpha    logbeta   logSigma 
    ##  20.777513   4.288593 -20.060854 -28.577762  -2.022785 
    ## 
    ## $objective
    ## [1] -16.30385
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $iterations
    ## [1] 52
    ## 
    ## $evaluations
    ## function gradient 
    ##       60       53 
    ## 
    ## $message
    ## [1] "relative convergence (4)"
