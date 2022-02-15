// --   TMB object for estimating maturity    --
// -- Written by Mahmood Jokar and Sam Subbey --
#include <TMB.hpp>


template <class Type>
Type posfun(Type x, double eps){
  if (x>=eps)
  {
    return x;
  }
  else
  {
    Type y=1.0-x/eps;
    return 0.5*eps*(1.+exp(-y));
  }
}



template<class Type>
Type objective_function<Type>::operator() ()
{
  // data----------------------------------------------------------
  DATA_INTEGER(min_age);
  DATA_INTEGER(max_age);
  DATA_INTEGER(start_year);
  DATA_INTEGER(end_year);
  DATA_INTEGER(length_l);
  DATA_VECTOR(meanlength);
  DATA_MATRIX(Nl);
  DATA_MATRIX(N);
  DATA_MATRIX(C);

  // parameters----------------------------------------------------
  PARAMETER(lnp1);
  PARAMETER(lnp2);
  PARAMETER(lnp3);
  PARAMETER(lnnu);

  // Convert parameter values from log scale-----------------------
  Type p1=exp(lnp1);
  Type p2=exp(lnp2);
  Type p3=exp(lnp3);
  Type nu=exp(lnnu);

  // define counters and other temporary variables-----------------
  int i,k,a,l;
  int years=end_year-start_year+1;
  int ages=max_age-min_age+1;

  //Define some needed vectors and matrices -----------------------
  vector<Type> m(length_l);
  matrix<Type> Nimmature(years,ages);
  matrix<Type> simN(years,ages);
  matrix<Type> simNmonth(years,12);

  //----------------------------------------------------------------
  for (i=0;i<length_l;i++){
    m(i)=1.0/(1.0+exp(4*p1*(p2-meanlength(i))));
  }

  for(i=0;i<years;i++){
    a=0;                     //min_age
    Type tmp=0;
    for(l=0;l<length_l;l++){
      int b=i*length_l + l;
      tmp += (1-m(l))*Nl(b,a);
    }
    Nimmature(i,0)=tmp;
  }
  for(a=1;a<ages;a++){
    Type tmp=0;
    for(l=0;l<length_l;l++){
      tmp += (1-m(l))*Nl(l,a);
    }
    Nimmature(0,a)=tmp;
  }

  // Initialize nll-----------------------------------------------
  Type nll = 0.0;

  // Initialize simN for the first year and age class-------------
  for (i=0;i<years;i++){
    simN(i,0)=Nimmature(i,0);
  }
  for (a=1;a<ages;a++){
    simN(0,a)=Nimmature(0,a);
  }
  for(int a = 0; a < ages-1; a++){//int a=0; // a should not be fixed
  for (i=0;i<(years-1);i++) {

    int d=i*ages+a;
    simNmonth(i,0)=(simN(i,a)*exp(-0.5*p3)-C(d,0))*exp(-0.5*p3);
    for (k=2;k<=3;k++) {
      simNmonth(i,k-1)=(simNmonth(i,k-2)*exp(-0.5*p3)-C(d,k-1))*exp(-0.5*p3);
    }
    for (k=4;k<=12;k++) {
      int d=(i+1)*ages+(a+1);
      simNmonth(i,k-1)=(simNmonth(i,k-2)*exp(-0.5*p3)-C(d,k-1))*exp(-0.5*p3);
    }
    simN(i+1,a+1)=posfun(simNmonth(i,11),0.00005);
  }
  }

  // Compute negative loglikelihood----------------------------------
  for (i=1;i<years;i++) {
    for (a=1;a<ages;a++) {
      nll += - nu*log(nu*N(i,a)/simN(i,a)) + lgamma(nu) + nu*N(i,a)/simN(i,a) + log(N(i,a));
    }
  }
  ADREPORT(p1);
  ADREPORT(p2);
  ADREPORT(p3);
  ADREPORT(nu);
  REPORT(N);
  REPORT(simN);
  return(nll);
}

