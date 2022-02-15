// --   TMB object for estimating maturity    --
// -- Written by Sondre Holleland --

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // -- COD data --
  DATA_MATRIX(Cco);  // Empirical consumption by year (row) and cod age A (col)
  DATA_MATRIX(Nco);  // Number of cod at start of year (row) and cod age A (col)
  DATA_MATRIX(Oco);  // Proportion of mature cod at start of year (row) and cod age A (col)
  DATA_MATRIX(SV);   // Proportion of immature cod in Svalbard area (year, age)
  DATA_MATRIX(Wco);  // Weight at age cod (year, age)
  DATA_MATRIX(Fco);  // Fishing mortality of cod (year, age)
  DATA_MATRIX(Mco);  // Natural mortality of cod (year, age)

  // -- Capelin data --
  DATA_MATRIX(Nmc); // Number at age of capelin (year, age)
  DATA_MATRIX(Wmc); // Weight at age capelin (year, age)
  DATA_VECTOR(p3);  // Natural mortality for sep-jan (by year)



  // -- Parameters --
  PARAMETER(logCmax);  // C_{max}
  PARAMETER(logChalf); // C_{1/2}
  PARAMETER(logalpha);
  PARAMETER(logbeta);
  PARAMETER(logSigma);


  // -- Parameter transformations --
  Type Cmax   = exp(logCmax);
  Type Chalf  = exp(logChalf);
  Type beta   = exp(logbeta);
  Type alpha  = exp(logalpha);
  Type sigma  = exp(logSigma);

  // Constants
  int nyears  = Mco.rows(); // number of years
  int nagecod = Mco.cols(); // number of cod ages
  int nt      = 3;          // number of time steps
  int nagecap = Nmc.cols(); // number of capelin ages

  // Initialize

  matrix<Type> MatBio(nyears,nt);     // Biomass of maturing capelin
  matrix<Type> PredAbility(nyears,nt);// Predation ability
  matrix<Type> Cons(nyears,nt);       // Consumption
  matrix<Type> ECons(nyears,nt);      // Empirical consumption
  matrix<Type> Zcod(nyears, nagecod); // Mortality cod
  matrix<Type> Mmc(nyears,nt);        // Mortality of maturing capelin

  // Cod suitability
  // - assumed to follow a beta distribution
  // - Should be zero for ages 1-2
  // - constant for age > 10
  vector<Type> CodSuit(nagecod);
  Type Atmp = 0.0;
  for(int A=0; A < nagecod; A++){
    if(A <= 1){
      CodSuit(A) = Type(0.0);
    }else if(A >= 10){
      CodSuit(A) = CodSuit(10);
    } else {
      Atmp = (Type(A+1) - Type(2)) / (Type(8.5));
      CodSuit(A) = dbeta(Atmp, alpha, beta, false);
    }
  }

  Type nll =  Type(0.0);         // Objective function
  vector<Type> conscolsum(nyears);  // consumption per year
  vector<Type> Econscolsum(nyears); // empirical consumption per year
  for(int y = 0; y<nyears; y++){
    Econscolsum(y) = Type(0.0);
    conscolsum(y)  = Type(0.0);
    for(int t = 0; t<nt; t++){

      if(t ==0){
        Mmc(y,t) = Type(0.0);
      }else{
        Mmc(y,t) = -log(Type(1.0) - Cons(y,t-1)/MatBio(y,t-1));
      }
      // MatBio:
      MatBio(y,t) =  Type(0.0);
      for(int a = 0; a < nagecap; a++){
        //Zcap(y,a) = Mmc(y,a) + Fmc(y,a); // mortality rate capelin
        MatBio(y,t) += Nmc(y,a) * exp(-Type(3) * p3(y)/Type(12)-Mmc(y,t)) * Wmc(y,a);
      }

      // PredAbility and empirical consumption:
      PredAbility(y,t) =  Type(0.0);
      ECons(y,t) = Type(0.0);
      for(int A = 0; A<nagecod; A++){
        Zcod(y,A) = Mco(y,A) + Fco(y,A); // mortality rate cod
        PredAbility(y,t) += Nco(y,A) * exp(-(t-Type(0.5)) * Zcod(y,A)/ Type(12.0)) * CodSuit(A) *(Type(1.0) - Oco(y,A))*(Type(1.0)-SV(y,A)) * pow(Wco(y,A), Type(0.801));
        ECons(y,t) += Nco(y,A) * exp(-(t-Type(0.5)) * Zcod(y,A) / Type(12.0)) *(Type(1.0) - Oco(y,A)) * (Type(1.0) - SV(y,A)) * Cco(y,A);
      }

      Cons(y,t) = Cmax * PredAbility(y,t) * MatBio(y,t)/(Chalf + MatBio(y,t));
      // aggregate over months (total per year)
      Econscolsum(y) += ECons(y,t);
      conscolsum(y)  += Cons(y,t);
      // Contribution to objective function per year
      //nll -= dnorm(ECons(y,t),Cons(y,t), sigma, true);
    }
    nll -= dnorm(Econscolsum(y),conscolsum(y), sigma, true);
  }

  // Reports:
  ADREPORT(Cmax);
  ADREPORT(Chalf);
  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(sigma);
  ADREPORT(conscolsum);
  REPORT(CodSuit);
  REPORT(Cons);
  REPORT(ECons);
  REPORT(MatBio);
  REPORT(Econscolsum);
  REPORT(conscolsum);
  REPORT(Mmc);
  return(nll);
}

