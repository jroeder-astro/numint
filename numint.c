#include<math.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>


// LIBRARY FOR NUMERICAL INTEGRATION


double gaussian(double x, void *params){

  double *par = (double *)params;
  double mu = par[0];
  double sigma = par[1];

  return 1./(sigma*sqrt(2.*M_PI))*exp(-pow(((x-mu)/sigma), 2.)/2.);
}


double somecos(double x, void *params){

  return x*pow(cos(2.*M_PI*x*x), 2.);
}

//double int_1(double a, double b, void *p, double (*f)(double, void *), double e){

//}




double int_left_riemann(double a, double b, void *p, double (*f)(double, void *)) {  // e is error
  // HOW DOES ONE PUT THE PARAMETERS IN HERE???   -> this apparently works lol  
  // double *p = (double*)params;  // Line not needed if void *p instead of void *params

  // double N = 10000.;
  double N =10000.;
  double h = (b - a) / N;

  // left Riemann sum
  double L = (*f)(a, p);                    // this bitch is the reason you gotta start with a+h in the for loop
  for (double i = a + h; i <= b - h; i += h)
  {
    L += (*f)(i, p);
  }

  L *= h;
  printf("Left sum of function is %+6.10lf\n", L);
}

double int_right_riemann(double a, double b, void *p, double (*f)(double, void *)) {
  // right Riemann sum
  double N =10000.;
  double h = (b - a) / N;
  double R = 0;
  
  for (double i = a + h; i <= b; i += h)
  {
    R += (*f)(i, p);
  }

  R *= h;
  printf("Right sum of function is %+6.10lf\n", R);
}

double int_trapezoidal(double a, double b, void *p, double (*f)(double, void *)) {
  // trapezoidal rule
  double N =10000.;
  double h = (b - a) / N;
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated

  
  for (double i = a + h; i <= b - h; i += h)
  {
    T += 2. * (*f)(i, p);
  }

  T *= h/2.;

  printf("Trapezoidal int. of function is %+6.10lf\n", T);

}

double int_simpson_one_loop(double a, double b, void *p, double (*f)(double, void *)){
  // Simpson's rule
  double N =100000.;
  double h = (b - a) / N;
  double S = (*f)(a,p)+(*f)(b,p);
  for (double i = a+h; i <= b-h; i += h){
    S += 2. * (*f)(i,p);
    S += 4. * (*f)(i-h/2.,p);   // we encouter at this postion that there is a difference between this an the older Veriosn
    // with two for loops.

  }

  S += 4. * (*f)(b-h/2.,p);

  S *= h/6.;

  printf("Simpson's rule gives us %+6.10lf\n", S);

}

double int_simpson_two_loop(double a, double b, void *p, double (*f)(double, void *)){
  // Simpson's rule
  double N =100000.;
  double h = (b - a) / N;
  double S = (*f)(a,p)+(*f)(b,p);
  for (double i = a+h; i <= b-h; i += h) {
    S += 2. * (*f)(i, p);
  }
  for (double i = a+h/2.; i<=b-h/2.; i+=h){


    S += 4. * (*f)(i,p);   // we encouter at this postion that there is a difference between this an the older Veriosn
    // with two for loops.

  }

  S *= h/6.;

  printf("Simpson's rule gives us %+6.10lf\n", S);

}



int main(){

  double x;
  double p[2] = {0., 1.}; // array with mu and sigma
  gaussian(x, p);
//  int_1(-1., 1., p, gaussian);
  int_left_riemann(-1.,1.,p, gaussian);
  int_right_riemann(-1.,1.,p, gaussian);
  int_trapezoidal(-1.,1,p, gaussian);
  int_simpson_one_loop(-1.,1,p,gaussian);
  int_simpson_two_loop(-1.,1,p,gaussian);

//  int_1(0., 2., NULL, somecos);

}

