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


double adapt_step_mid(double a, double b, void *p, double (*f)(double, void *), double e){

  double rel = 1.;
  double N = 1000.;
  double K = 3.;
  double h = (b-a)/N;
  
  double M = 0.;
  double M1 = 0.;

  printf("start first for");

  for (double i = a+h/2.; i <= b-h/2.; i += h)
  {
  	M += h*(*f)(i,p);
  }

  printf("start while");

  while (rel >= e)
  {
	M1 = 1./3. * M;
 	for (double i=a+h/(2.* K); i <= b-h/(2.* K); i += 2.*h/K)
	{
		M1 = h/K * (*f)(i,p);
	}

  	rel = fabs(M-M1)/fabs(M1);  // calculate relative error

  	K *= 3.; // halving the stepsize
	M = M1; 
  }
  printf("Midoint rule. WITH STEP TRIPLING OMG!!! This gives us: M = %+6.10lf \n", M);

}



/*
double adapt_step_trap(double a, double b, void *p, double (*f)(double, void *), double e){

  // relative error e>0 
  
  double rel = 1.; 		// initialize relative error 
  double N = 1.;   		// initialize number of steps to start with for initial stepwidth
  double K = 2.;   		// initialize halving parameter
  double h = (b - a) / N;	// initialize stepwidth
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated start value 

  double T1 = 0.;  // initialize halved stepsize value

  for (double i = a + h; i <= b - h; i += h)  // value of integral before stepsize halving, trapezoidal method
  {
    	T += 2. * (*f)(i, p);
  }

  T *= h/2.;

  while (rel >= e) // while loop for error control; runs while relative error is greater than given error 
  {  
 	T1 = 1./2. * T;    // calculate value with halved stepsize value 
	for (double i = a + h/K; i <= b-h/K; i += 2. * h/K)
	{
		T1 += h/K * (*f)(i,p);
	} 

	rel = fabs(T-T1)/fabs(T1);  // calculate relative error

  	K *= 2.; // halving the stepsize
	T = T1; 
  }
  printf("Trapezoidal method, enhanced with stepsize halving (super amazing!), gives us T = %+6.10lf \n", T);
}
*/


/*

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

*/

int main(){
  printf("start of main");
  double x;
  double p[2] = {0., 1.}; // array with mu and sigma
  gaussian(x, p);
  adapt_step_mid(-1., 1., p, gaussian, 0.01);

//  int_1(-1., 1., p, gaussian);
//  int_left_riemann(-1.,1.,p, gaussian);
//  int_right_riemann(-1.,1.,p, gaussian);
//  int_trapezoidal(-1.,1,p, gaussian);
//  int_simpson_one_loop(-1.,1,p,gaussian);
//  int_simpson_two_loop(-1.,1,p,gaussian);

//  int_1(0., 2., NULL, somecos);

//  adapt_step_trap(-1., 1., p, gaussian, 0.0001);
  return 0 ; 
}

