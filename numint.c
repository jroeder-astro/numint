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


double quadexp(double x, void *params){

  return exp(-pow(x, 2.));
}


double inverse_sqrt(double x, void *params){

  return 1./sqrt(x);
}


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


double adapt_step_mid(double a, double b, void *p, double (*f)(double, void *), double e, 
                      double(*trafo)(double(*)(double, void *), double, void *)){
  // This function now takes an extra argument "trafo"
  // Which usually is an indentity_trafo which only gives back the function itself.
  double rel = 1.;
  double N = 1000.;
  double K = 3.;
  double h = (b-a)/N;
  
  double M = 0;
  double M1 = 0.;
  
  if (a == b){
        // "empty" integralls shall not be evaluated numerically

	return 0;
  }  

  for (double i = a+h/2.; i <= b-h/2.; i += h) // initializes M with midpoint rule
  {
  	M += h*(*trafo)(*f,i,p);
  }
  printf("this is M = %6.10lf \n", M);

  int cnt = 1; // counter used for stepsize in iterations > h/3. "eval, not, eval, eval, not, e, e, n,...."  

  while (rel >= e)
  {
	cnt = 1;
	M1 = 1./3. * M;
 	for (double i=a+h/(2.* K); i <= b-h/(2.* K); i+=0 ) // analytically derived formula for step tripling
	{
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


// Optional task: singularity integration
// Should at some point be revamped as parameters are technically being misused

double sing_int(double a, double b, void *p, double(*f)(double, void*), double e){
    
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
  printf("start of main\n");
  double x;
  double p[2] = {0., 1.}; // array with mu and sigma
  gaussian(x, p);

  double q[2] = {0.5,0.}; // order of singularity

  sing_int(0., 1., q, inverse_sqrt, 0.000001);

//  infty_bound(0, 1, NULL, quadexp, 0.01);

//  adapt_step_mid(0., 2., NULL, somecos, 0.99, identity_trafo);
//  adapt_step_mid(0., 2., NULL, somecos, 0.0001, identity_trafo);

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

