#include<math.h>
#include<time.h>
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
  double N = 300000.;
  double K = 3.;
  double h = (b-a)/N;
  
  double M = 0;
  double M1 = 0.;

  for (double i = a+h/2.; i <= b-h/2.; i += h) // initializes M with midpoint rule
  {
  	M += h*(*f)(i,p);
  }
  printf("this is M = %6.10lf", M);

  int cnt = 1; // counter used for stepsize in iterations > h/3. "eval, not, eval, eval, not, e, e, n,...."  

  while (rel >= e)
  {
	cnt = 1;
	M1 = 1./3. * M;
 	for (double i=a+h/(2.* K); i <= b-h/(2.* K); i+=0 ) // analytically derived formula for step tripling
	{
		M1 += h/K * (*f)(i,p);
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


double montecarlo(double a, double b, void *p,  double (*f)(double, void *), double abe){

  double N = 1000.;
  double h = (b-a)/N;

  double result = 0.;
  double error = abe;
  double rndm = 0.;
  double dummy_eval = 0.;

  time_t t ;
  struct tm tm;

  srand(1.);

  while (error >= abe)
  {
	 t = time(NULL);
	 tm = *localtime(&t);
	 error = 0;
	 result = 0;
 	 printf("We are at N = %lf\n at : ", N);
         printf("now %d:%d:%d\n",tm.tm_hour, tm.tm_min, tm.tm_sec);	 

 	 for (double i = a; i <= b; i += h) 
  	 {
		double partial_sum = 0.;
		double partial_sum_square = 0;
		double partial_error = 0.;
		for (double l = 1; l <= N; l++)
		{
			rndm = (double)rand()/RAND_MAX*h;
			dummy_eval = (*f)(i+rndm,p);

			partial_sum += dummy_eval;
 			partial_sum_square +=  pow(dummy_eval,2.);
		
		}
		partial_sum /= N;
		partial_sum_square /= N;
		partial_error = sqrt((partial_sum_square - pow(partial_sum,2.))/N);
		
		result += partial_sum * h;
		error += partial_error *h;
 	 }
  N *= 2.;
  }
  
  printf("Monte Carlo integral gives us %6.10lf +- %6.10lf \n", result, error);
  return result;
}




int main(){
  printf("start of main\n");
  double x;
  double p[2] = {0., 1.}; // array with mu and sigma
  gaussian(x, p);

  montecarlo(-1., 1., p, gaussian, 0.000001);


//  adapt_step_mid(0., 2., NULL, somecos, 0.99);
//  adapt_step_mid(0., 2., NULL, somecos, 0.0001);

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

