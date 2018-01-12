#include<math.h>
#include<time.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>


// =================================
// LIBRARY FOR NUMERICAL INTEGRATION
// =================================


/*

Commit message:
added return statements, turned double loops into int loops, removed old calls in main

Fixed loops in:    (sing_int & infty_bound have no for loops)

trapezoidal
left sum i
right sum
simpson rule (1 loop version)
monte_carlo

adapt_step_mid (hopefully)
adapt_step_trap, first for loop
adapt_step_simp, first for loop

*/


// FUNCTIONS

double gaussian(double x, void *params){

  double *par = (double *)params;
  double mu = par[0];
  double sigma = par[1];

  return 1./(sigma*sqrt(2.*M_PI))*exp(-pow(((x-mu)/sigma), 2.)/2.);
}


double somecos(double x, void *params){

  return x*pow(cos(2.*M_PI*x*x), 2.);
}


double quadexp(double x, void *params){

  return exp(-pow(x, 2.));
}


double inverse_sqrt(double x, void *params){

  return 1./sqrt(x);
}


// TRANSFORMATIONS


double identity_trafo(double(*f)(double, void *), double x, void *p){
  // Gives back the function itself
  // Needed for non-infinite boundaries
  // For normal integration, GIVE THIS TO ADAPT_STEP_MID!
  return (*f)(x,p);
}


double inverse_square_trafo(double(*f)(double, void *), double x, void *p){
  // Gives back function called with inverse variable multiplied by inverse variable squared
  // used for infinity boundary algorithm
  return 1./(x*x)* (*f)(1./x, p);
}


double singularity_trafo(double(*f)(double, void *), double x, void *p){
  // Implementing given transformation to not integrate over "1/0" 
  double *par = (double *)p;
  double g = par[0]; // parameter use...
  double a = par[1];
  return  pow(x, g/(1.-g)) * (*f)(pow(x, 1/(1.-g))+a,p);
}


// INTEGRATION METHODS


double adapt_step_mid(double a, double b, void *p, double (*f)(double, void *), double e, 
                      double(*trafo)(double(*)(double, void *), double, void *)){

  // This function now takes an extra argument "trafo"
  // Which usually is an indentity_trafo which only gives back the function itself.
  if (a==b)
  {
      return 0. ;
  }

  if(b < a)
  {
      return -1. * adapt_step_mid(b, a, p, f, e, trafo);
  }

  double rel = 1.;
  int N = 1000;
  double K = 3.;
  double h = (b-a)/(double)N;
  
  double M = h*(*trafo)(*f,a+h/2.,p);
  double M1 = 0.;
  
  if (a == b){
        // "empty" integralls shall not be evaluated numerically

	return 0;
  }  

  for (int i = 1; i <= N-1; i ++) // initializes M with midpoint rule
  {
        double m = a+i*h/2.+i*h;
  	M += h*(*trafo)(*f,m,p);
  }
  printf("this is M = %6.10lf \n", M);

  int cnt = 1; // counter used for stepsize in iterations > h/3. "eval, not, eval, eval, not, e, e, n,...."  

  while (rel >= e)
  {
	cnt = 1;
	M1 = 1./3. * M;
 	for (double i=a+h/(2.* K); i <= b-h/(2.* K); i+=0 ) // analytically derived formula for step tripling
	{
		// We do not actually loop here, the "loop" is executed by an if statement to
		// implement the e,n,e,e,n,e,e,n,...,n,e law. It might not matter here if
		// there is a double in the for loop.

		M1 += h/K * (*trafo)(*f,i,p);
		if (cnt%2==0)  // doing the e, n, e, e, n, e, e, n,... stuff
		{
			i += h/K;
		}
		else
		{
			i += 2.*h/K; 
		}
		
		cnt += 1;
	}

  	rel = fabs(M-M1)/fabs(M1);  // calculate relative error

  	K *= 3.; // tripling the stepsize
	M = M1; 
  }
  printf("Midoint rule. WITH STEP TRIPLING OMG!!! This gives us: M = %+6.10lf \n", M);
  printf("cnt = %d\n" ,cnt);
  return M;
}


double infty_bound(double a, int isinf, void *p, double (*f)(double, void *), double e){
  // Used if upper bound is infinity
  
  if (isinf == 0) 
  {
  	// Check if upper bound really is infinity
	// One meaning it is infinity, when 0 it is not.
	printf("a must be greater than zero and the upper must equal infinity.");
  	return 0.;
  }

  double result = 0;

  if (a <= 0)
  {
	// Splitting Integral for non-zero lower bound
	// Also avoids upper < lower bound after trafo
	result += adapt_step_mid(a, 1, p, f, e, identity_trafo);
	a = 1.;
  }
  
  double b = 1./a;
  a = 0.;
  result += adapt_step_mid(a, 1, p, f, e, inverse_square_trafo);  
  
  printf("Midoint rule. WITH STEP TRIPLING OMG!!! This gives us: M = %+6.10lf \n", result);
  return result;
}


double adapt_step_trap(double a, double b, void *p, double (*f)(double, void *), double e){

  if(a==b)
  {
      return 0. ;
  }

  if(b < a)
  {
      return -1. *adapt_step_trap(b, a, p, f,e);
  }
  // relative error e>0

  double rel = 1.; 		// initialize relative error 
  int N = 1000;   		// initialize number of steps to start with for initial stepwidth
  double K = 2.;   		// initialize halving parameter
  double h = (b - a) / (double)N;	// initialize stepwidth
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated start value 

  double T1 = 0.;  // initialize halved stepsize value

  for (int i = 1; i <= N-1; i ++)  // value of integral before stepsize halving, trapezoidal method
  {
	double m = a+i*h;
    	T += 2. * (*f)(m, p);
  }

  T *= h/2.;

  while (rel >= e) // while loop for error control; runs while relative error is greater than given error 
  {

      T1 = 1./2. * T;    // calculate value with halved stepsize value
	    for (int i = 0.; i <= /* N-1   */ ; i ++)  // That .../K thing is difficult to get into an int loop
	    { 		      // This N-1 is wrong, but no idea how to put a+h/k......b-h/K into this loop,
			      // with stepsize 2h/K
		      double m = a + h/K + i * 2.*h/K  ;
		      T1 += h/K * (*f)(i,p);
	    }

	    rel = fabs(T-T1)/fabs(T1);  // calculate relative error

  	  K *= 2.; // halving the stepsize
	    T = T1;
  }
  printf("Trapezoidal method, enhanced with stepsize halving (super amazing!), gives us T = %+6.10lf \n", T);
  return T;
}


double int_left_riemann(double a, double b, void *p, double (*f)(double, void *)) {  // e is error
  // HOW DOES ONE PUT THE PARAMETERS IN HERE???   -> this apparently works lol  
  // double *p = (double*)params;  // Line not needed if void *p instead of void *params

  if (a==b)
  {
      return 0. ;
  }

  if ( b < a )
  {
      return -1. * int_left_riemann(b, a,p, f);
  }

  // double N = 10000.;
  int N =10000;
  double h = (b - a) / (double)N;

  // left Riemann sum
  double L = (*f)(a, p);           // this bitch is the reason you gotta start with a+h in the for loop
  for (int i = 1; i <=N-1; i ++)
  {
    double m = a+i*h;
    L += (*f)(m, p);
  }

  L *= h;
  printf("Left sum of function is %+6.30lf\n", L);
  return L;
}


double int_right_riemann(double a, double b, void *p, double (*f)(double, void *)) {
  // right Riemann sum
  if (a==b)
  {
    return 0. ;
  }

  if (b < a)
  {
      return -1. * int_right_riemann(b, a, p, f);
  }

  int N =10000;
  double h = (b - a) / (double)N;
  double R = 0;
  
  for (int i = 1; i <= N; i ++)
  {
    double m = a+i*h;
    R += (*f)(m, p);
  }

  R *= h;
  printf("Right sum of function is %+6.30lf\n", R);
  return R;
}


double int_trapezoidal_double(double a, double b, void *p, double (*f)(double, void *)) {
  // trapezoidal rule

  if (a==b)
  {
      return 0. ;
  }

  if ( b < a )
  {
      return -1. * int_trapezoidal_double(b, a, p, f);
  }

  int N =10000;
  double h = (b - a) / (double)N;
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated

  
  for (double i = a + h; i <= b - h; i += h)
  {
    T += 2. * (*f)(i, p);
  }

  T *= h/2.;

  printf("Trapezoidal int. of funct., doublez loop, is %+6.30lf\n", T);
  return T;
}


double int_trapezoidal_int(double a, double b, void *p, double (*f)(double, void *)) {

  // trapezoidal rule

  if (a == b)
  {
      return 0. ;
  }

  if (b < a)
  {
    return -1. * int_trapezoidal_int(b, a, p, f);
  }
  int N =10000;
  double h = (b - a) / (double)N;
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated

  
  for (int i = 1; i <= N-1; i ++)
  {
    double m = a+i*h;
    T += 2. * (*f)(m, p);
  }

  T *= h/2.;

  printf("Trapezoidal int. of funct., integer loop, is %+6.30lf\n", T);
  return T;
}


double int_simpson_one_loop(double a, double b, void *p, double (*f)(double, void *)){
  // Simpson's rule
  if (a==b)
  {
      return 0.;
  }

  if (b < a)
  {
      return -1. * int_simpson_one_loop(b, a, p, f);
  }

  int N =100000;
  double h = (b - a) / (double)N;
  double S = (*f)(a,p)+(*f)(b,p);
  for (int i = 1; i <= N-1; i ++){

    double m = a+i*h;
    S += 2. * (*f)(m,p);
    S += 4. * (*f)(m-h/2.,p);   // We encouter at this postion that there is a difference 
			 	// between this and the older version with two for loops.

  }

  S += 4. * (*f)(b-h/2.,p);

  S *= h/6.;

  printf("Simpson's rule gives us %+6.10lf\n", S);
  return S;
}


double int_simpson_two_loop(double a, double b, void *p, double (*f)(double, void *)){
  // Simpson's rule // OUTDATED!!!
  // VERY crappy to do with two loops using m instead of i, DO NOT TOUCH THIS FUNCTION
  if (a==b)
  {
       return 0.;
  }

  if (b < a)
  {
      return -1.* int_simpson_two_loop(b, a, p, f);
  }

  int N =100000;
  double h = (b - a) / (double)N;
  double S = (*f)(a,p)+(*f)(b,p);
  for (int i = 1; i <= N-1; i ++) {
    
    double m = a+i*h;
    S += 2. * (*f)(m, p);
  }
  for (double i = a+h/2.; i<=b-h/2.; i+=h){

    double m = a+i*h;

    S += 4. * (*f)(m,p);   // We encouter at this position that there is a difference 
			   // between this and the older version with two for loops.

  }

  S *= h/6.;

  printf("Simpson's rule gives us %+6.10lf\n", S);
  return S;
}


double montecarlo(double a, double b, void *p,  double (*f)(double, void *), double abe){

  if (a == b)
  {
      return 0.;
  }

  if (b < a)
  {
      return -1. * montecarlo(b, a, p, f, abe);
  }
  int N = 1000;
  double h = (b-a)/ (double)N;

  double result = 0.;
  double error = abe;
  double rndm = 0.;	// random number initialization
  double dummy_eval = 0.; // dummy to avoid too many evaluations in innermost for loop

  time_t t ;    // this is to give time 
  struct tm tm;

  srand(time(NULL));  // RNG seed

  while (error >= abe)  // error control
  {
	 t = time(NULL);  
	 tm = *localtime(&t);
	 error = 0;
	 result = 0;
 	 printf("We are at N = %d\n at : ", N);    // keeping track of calculations is awesome
         printf("now %d:%d:%d\n",tm.tm_hour, tm.tm_min, tm.tm_sec);	 

 	 for (int i = 0; i <= N; i ++) 
  	 {
	    	double m = a+i*h;

		double partial_sum = 0.;    // sum part of <f> 
		double partial_sum_square = 0.;   // sum part of <f²>
		double partial_error = 0.;   // sqrt part of error
				// above three lines are for one interval only 
		for (int l = 1; l <= N; l++)
		{
			rndm = (double)rand()/RAND_MAX*h;    // actual RNG
			dummy_eval = (*f)(m + rndm,p);         // dummy calculation

			partial_sum += dummy_eval;  // building up sum part of <f>
 			partial_sum_square +=  pow(dummy_eval,2.);  // building <f²>
		
		}
		partial_sum /= N;   // calculating actual <f>
		partial_sum_square /= N;   // calculating <f²>
		partial_error = sqrt((partial_sum_square - pow(partial_sum,2.))/N); // calculating error
		
		result += partial_sum * h;   // calculation end result
		error += partial_error *h;
 	 }
  N *= 2.;   // increasing number of samples to get a more precise result
  }
  
  printf("Monte Carlo integral gives us %6.10lf +- %6.10lf \n", result, error);
  return result;
}


// Optional task: singularity integration
// Should at some point be revamped as parameters are technically being misused

double sing_int(double a, double b, void *p, double(*f)(double, void*), double e){

  if (a==b)
  {
    return 0.;
  }

  if (b < a)
  {
      return -1. * sing_int(b, a, p, f, e);
  }

  double *par = (double *)p; // parameters should not be used that way
  double g = par[0];

  if (g == 1.) // if g=0, one would divide by zero
  {
  	printf("g is not allowed to be one.\n");
	return 0.;
  }

  double result = 0.;
  b = pow((b-a), (1.-g)); // new upper bound after trafo
  a = 0.; 		  // new lower bound after trafo
  
  result += 1./(1.-g) * adapt_step_mid(a, b, p, f, e, singularity_trafo);
  printf("Singularity integration gives us: S = %+6.10lf\n", result);
  return result;

}


// Optional task: Simpson's rule with stepsize halving

double adapt_step_simp(double a, double b, void *p, double(*f)(double, void *), double e){

  if (a == b)
  {
    return 0. ;
  }

  if (b < a)
  {
      return -1. * adapt_step_simp(b, a, p, f, e);
  }

  double S = (*f)(a,p)+(*f)(b,p); // intital value, derived analytically
  int N = 100000;
  double K = 2.;
  double h = (b - a) / (double)N;
  double S1 = 0.;

  double dummy_subs = 0.; // we need a substitute variable since we derived analytically that the following
			  // while loop can produce stepsize halfing with the former value minus this subs
  double dummy_eval = 0.; // need this for evaluation efficiency

  double rel = 1.;	  // relative error initiation
  // nor
  for (int i = 1; i <= N-1; i ++)
  {
	double m = a+i*h;

	S += 2. * (*f)(m,p);
	dummy_eval = (*f)(m - h/2.,p);
  	S += 4. * dummy_eval;
	dummy_subs += dummy_eval;

  }
  // Trick needed to avoid second loop
  dummy_eval = (*f)(b-h/2.,p);
  
  S += 4. * dummy_eval;
  dummy_subs += dummy_eval;

  S *= h/6.;
  while (rel >= e)
  {
	// analytically derived formula
	S1 = 1./2. * (S - h/3. * dummy_subs);
	for (double i = a+h/(2.*K); i <= b-h/(2.*K); i += h/K)
	{			// Same problem as in the previous adapt_step_... function.
				// How does one convert this for loop into an int loop?
		dummy_eval += (*f)(i,p);
	}

 	S1 += 2. * h/ (3. * K) * dummy_eval; // analytically derived
  	dummy_subs += h/3. *  dummy_eval; // analytically derived
	// normal stepsize halfing process
	rel = fabs(S1-S)/fabs(S1);
	S = S1;
	K *= 2.;
  }

  printf("Simpson's method, now with stepsize halving! Only 99 Cents!! Only Here!! Gives: %+6.10lf\n", S);
  return S;

}



int main(){
  
  double x;

  double p[2] = {0., 1.}; // array with mu and sigma
  
  double q[2] = {0.5,0.}; // order of singularity


//  adapt_step_simp(-1., 1., p, gaussian, 0.01);

  
//  sing_int(0., 1., q, inverse_sqrt, 0.000001);


//  infty_bound(0, 1, NULL, quadexp, 0.01);


//  adapt_step_mid(0., 2., NULL, somecos, 0.99, identity_trafo);
//  adapt_step_mid(0., 2., NULL, somecos, 0.0001, identity_trafo);

  
//  montecarlo(-1., 1., p, gaussian, 0.000001);


  int_left_riemann(-1.,1.,p, gaussian);
  int_right_riemann(-1.,1.,p, gaussian);

  int_trapezoidal_int(-1.,1,p, gaussian);
  int_trapezoidal_double(-1.,1,p, gaussian);


//  int_simpson_one_loop(-1.,1,p,gaussian);
//  int_simpson_two_loop(-1.,1,p,gaussian);


//  adapt_step_trap(-1., 1., p, gaussian, 0.0001);
 
 return 0 ; 
}

