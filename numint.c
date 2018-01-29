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

adapt_step_mid  done
adapt_step_trap, done
adapt_step_simp, done

*/


// FUNCTIONS

double gaussian(double x, void *params) {

  double *par = (double *) params;
  double mu = par[0];
  double sigma = par[1];

  return 1. / (sigma * sqrt(2. * M_PI)) * exp(-pow(((x - mu) / sigma), 2.) / 2.);
}


double somecos(double x, void *params) {

  return x * pow(cos(2. * M_PI * x * x), 2.);
}


double quadexp(double x, void *params) {

  return exp(-pow(x, 2.));
}


double inverse_sqrt(double x, void *params) {

  return 1. / sqrt(x);
}


// TRANSFORMATIONS


double identity_trafo(double(*f)(double, void *), double x, void *p) {
  // Gives back the function itself
  // Needed for non-infinite boundaries
  // For normal integration, GIVE THIS TO ADAPT_STEP_MID!
  return (*f)(x, p);
}


double inverse_square_trafo(double(*f)(double, void *), double x, void *p) {
  // Gives back function called with inverse variable multiplied by inverse variable squared
  // used for infinity boundary algorithm
  return 1. / (x * x) * (*f)(1. / x, p);
}


double singularity_trafo(double(*f)(double, void *), double x, void *p) {
  // Implementing given transformation to not integrate over "1/0" 
  double *par = (double *) p;
  double g = par[0]; // parameter use...
  double a = par[1];
  return pow(x, g / (1. - g)) * (*f)(pow(x, 1 / (1. - g)) + a, p);
}


// INTEGRATION METHODS


double adapt_step_mid(double a, double b, void *p, double (*f)(double, void *), double e,
                      double(*trafo)(double(*)(double, void *), double, void *)) {

  // This function now takes an extra argument "trafo"
  // Which usually is an indentity_trafo which only gives back the function itself.
  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * adapt_step_mid(b, a, p, f, e, trafo);
  }

  double rel = 1.;
  int N = 1000;
  double K = 3.;
  double h = (b - a) / (double) N;

  double M = h * (*trafo)(*f, a + h / 2., p);
  double M1 = 0.;

  if (a == b) {
    // "empty" integralls shall not be evaluated numerically

    return 0;
  }

  double m = a + h/2.;

  for (int i = 1; i <= N - 1; i++) // initializes M with midpoint rule
  {

    M += h * (*trafo)(*f, m, p);
    m += h;
  }


  // counter used for stepsize in iterations > h/3. "eval, not, eval, eval, not, e, e, n,...."

  while (rel >= e) {

    M1 = 1. / 3. * M;
    m = a+ h/(2.*K);
    // for (double i = a + h / (2. * K); i <= b - h / (2. * K); i += 0) // analytically derived formula for step tripling

    for (int i = 1; i <= 2*N*K/3. ;i++)
    {
      printf("i=%d", i);
      // We do not actually loop here, the "loop" is executed by an if statement to
      // implement the e,n,e,e,n,e,e,n,...,n,e law. It might not matter here if
      // there is a double in the for loop.

      M1 += h / K * (*trafo)(*f, m, p);
      if (i % 2 == 0)  // doing the e, n, e, e, n, e, e, n,... stuff
      {
        m += h / K;
      }
       else {
        m += 2. * h / K;
      }

      cnt += 1;
    }

    rel = fabs(M - M1) / fabs(M1);  // calculate relative error

    K *= 3.; // tripling the stepsize
    M = M1;
  }

  return M;
}


double infty_bound(double a, int isinf, void *p, double (*f)(double, void *), double e) {
  // Used if upper bound is infinity

  if (isinf == 0) {
    // Check if upper bound really is infinity
    // One meaning it is infinity, when 0 it is not.
    printf("a must be greater than zero and the upper must equal infinity.");
    return 0.;
  }

  double result = 0;

  if (a <= 0) {
    // Splitting Integral for non-zero lower bound
    // Also avoids upper < lower bound after trafo
    result += adapt_step_mid(a, 1, p, f, e, identity_trafo);
    a = 1.;
  }

  double b = 1. / a;
  a = 0.;
  result += adapt_step_mid(a, 1, p, f, e, inverse_square_trafo);

  return result;
}


double adapt_step_trap(double a, double b, void *p, double (*f)(double, void *), double e) {

  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * adapt_step_trap(b, a, p, f, e);
  }
  // relative error e>0

  double rel = 1.;    // initialize relative error
  int N = 1000;      // initialize number of steps to start with for initial stepwidth
  double K = 2.;      // initialize halving parameter
  double h = (b - a) / (double) N;  // initialize stepwidth
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated start value 

  double T1 = 0.;  // initialize halved stepsize value

  for (int i = 1; i <= N - 1; i++)  // value of integral before stepsize halving, trapezoidal method
  {
    double m = a + (double) i * h;
    T += 2. * (*f)(m, p);
  }

  T *= h / 2.;

  double m;

  while (rel >= e) // while loop for error control; runs while relative error is greater than given error 
  {

    T1 = 1. / 2. * T;    // calculate value with halved stepsize value

    h = (b - a) / (double) N; // recalculate the stepsize
    m = a + h;

    for (int i = 0; i <= N / 2 - 1; i++)   // That .../K thing is difficult to get into an int loop
    {
      // this formula is used to determine which Maximum i we need
      //  static double m = a + ( 1. + 2.* (double) i )*h ;
      // for i being N-1 ==> a + (1+ 2N-2)h = a + N-1 *h + N * h= b-h + N* h = b-h + b - a
      // for i being N/2 ==> a + (1 + N)h = b +h
      // for i being N/2-1 ==> a +(1 + N -2)h = b-h !!!!
      // but are there all needed steps in it ???
      T1 += h * (*f)(m, p);
      m += 2 * (double) h;
    }

    rel = fabs(T - T1) / fabs(T1);  // calculate relative error

    N *= 2; // halving the stepsize
    T = T1;
  }

  return T;
}


double int_left_riemann(double a, double b, void *p, double (*f)(double, void *)) {  // e is error
  // HOW DOES ONE PUT THE PARAMETERS IN HERE???   -> this apparently works lol  
  // double *p = (double*)params;  // Line not needed if void *p instead of void *params

  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * int_left_riemann(b, a, p, f);
  }

  // double N = 10000.;
  int N = 1000;
  double h = (b - a) / (double) N;

  // left Riemann sum
  double L = (*f)(a, p);           // this bitch is the reason you gotta start with a+h in the for loop
  for (int i = 1; i <= N - 1; i++) {
    double m = a + i * h;
    L += (*f)(m, p);
  }

  L *= h;

  return L;
}


double int_right_riemann(double a, double b, void *p, double (*f)(double, void *)) {
  // right Riemann sum
  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * int_right_riemann(b, a, p, f);
  }

  int N = 1000;
  double h = (b - a) / (double) N;
  double R = 0;

  for (int i = 1; i <= N; i++) {
    double m = a + i * h;
    R += (*f)(m, p);
  }

  R *= h;

  return R;
}


double int_trapezoidal_double(double a, double b, void *p, double (*f)(double, void *)) {
  // trapezoidal rule

  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * int_trapezoidal_double(b, a, p, f);
  }

  int N = 1000;
  double h = (b - a) / (double) N;
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated


  for (double i = a + h; i <= b - h; i += h) {
    T += 2. * (*f)(i, p);
  }

  T *= h / 2.;

  return T;
}


double int_trapezoidal_int(double a, double b, void *p, double (*f)(double, void *)) {

  // trapezoidal rule

  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * int_trapezoidal_int(b, a, p, f);
  }
  int N = 1000;
  double h = (b - a) / (double) N;
  double T = (*f)(a, p) + (*f)(b, p); //analytically evaluated


  for (int i = 1; i <= N - 1; i++) {
    double m = a + i * h;
    T += 2. * (*f)(m, p);
  }

  T *= h / 2.;

  return T;
}


double int_simpson_one_loop(double a, double b, void *p, double (*f)(double, void *)) {
  // Simpson's rule
  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * int_simpson_one_loop(b, a, p, f);
  }

  int N = 1000;
  double h = (b - a) / (double) N;
  double S = (*f)(a, p) + (*f)(b, p);
  for (int i = 1; i <= N - 1; i++) {

    double m = a + i * h;
    S += 2. * (*f)(m, p);
    S += 4. * (*f)(m - h / 2., p);                           // We encouter at this postion that there is a difference
    // between this and the older version with two for loops.

  }

  S += 4. * (*f)(b - h / 2., p);

  S *= h / 6.;

  return S;
}


double int_simpson_two_loop(double a, double b, void *p, double (*f)(double, void *)) {
  // Simpson's rule // OUTDATED!!!
  // VERY crappy to do with two loops using m instead of i, DO NOT TOUCH THIS FUNCTION
  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * int_simpson_two_loop(b, a, p, f);
  }

  int N = 1000;
  double h = (b - a) / (double) N;
  double S = (*f)(a, p) + (*f)(b, p);
  for (int i = 1; i <= N - 1; i++) {

    double m = a + i * h;
    S += 2. * (*f)(m, p);
  }
  for (double i = a + h / 2.; i <= b - h / 2.; i += h) {

    double m = a + i * h;

    S += 4. * (*f)(m, p);                                // We encouter at this position that there is a difference
    // between this and the older version with two for loops.

  }

  S *= h / 6.;

  return S;
}


double montecarlo(double a, double b, void *p, double (*f)(double, void *), double abe) {

  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * montecarlo(b, a, p, f, abe);
  }
  int N = 1000;
  double h = (b - a) / (double) N;

  double result = 0.;
  double error = abe;
  double rndm = 0.;  // random number initialization
  double dummy_eval = 0.; // dummy to avoid too many evaluations in innermost for loop

  time_t t;    // this is to give time
  struct tm tm;

  srand(time(NULL));  // RNG seed

  while (error >= abe)  // error control
  {
    t = time(NULL);
    tm = *localtime(&t);
    error = 0;
    result = 0;

    for (int i = 0; i <= N; i++) {
      double m = a + i * h;

      double partial_sum = 0.;    // sum part of <f>
      double partial_sum_square = 0.;   // sum part of <f²>
      double partial_error = 0.;   // sqrt part of error
      // above three lines are for one interval only
      for (int l = 1; l <= N; l++) {
        rndm = (double) rand() / RAND_MAX * h;    // actual RNG
        dummy_eval = (*f)(m + rndm, p);         // dummy calculation

        partial_sum += dummy_eval;  // building up sum part of <f>
        partial_sum_square += pow(dummy_eval, 2.);  // building <f²>

      }
      partial_sum /= N;   // calculating actual <f>
      partial_sum_square /= N;   // calculating <f²>
      partial_error = sqrt((partial_sum_square - pow(partial_sum, 2.)) / N); // calculating error

      result += partial_sum * h;   // calculation end result
      error += partial_error * h;
    }
    N *= 2.;   // increasing number of samples to get a more precise result
  }

  return result;
}


// Optional task: singularity integration
// Should at some point be revamped as parameters are technically being misused

double sing_int(double a, double b, void *p, double(*f)(double, void *), double e) {

  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * sing_int(b, a, p, f, e);
  }

  double *par = (double *) p; // parameters should not be used that way
  double g = par[0];

  if (g == 1.) // if g=0, one would divide by zero
  {
    printf("g is not allowed to be one.\n");
    return 0.;
  }

  double result = 0.;
  b = pow((b - a), (1. - g)); // new upper bound after trafo
  a = 0.;      // new lower bound after trafo

  result += 1. / (1. - g) * adapt_step_mid(a, b, p, f, e, singularity_trafo);
  return result;

}


// Optional task: Simpson's rule with stepsize halving

double adapt_step_simp(double a, double b, void *p, double(*f)(double, void *), double e) {

  if (a == b) {
    return 0.;
  }

  if (b < a) {
    return -1. * adapt_step_simp(b, a, p, f, e);
  }

  double S = (*f)(a, p) + (*f)(b, p); // intital value, derived analytically
  int N = 1000;
  double K = 2.;
  double h = (b - a) / (double) N;
  double S1 = 0.;

  double dummy_subs = 0.; // we need a substitute variable since we derived analytically that the following
  // while loop can produce stepsize halfing with the former value minus this subs
  double dummy_eval = 0.; // need this for evaluation efficiency

  double rel = 1.;    // relative error initiation
  // nor
  for (int i = 1; i <= N - 1; i++) {
    double m = a + i * h;

    S += 2. * (*f)(m, p);
    dummy_eval = 4 * (*f)(m - h / 2., p);
    S += dummy_eval;
    dummy_subs += dummy_eval;

  }
  // Trick needed to avoid second loop
  dummy_eval = 4 * (*f)(b - h / 2., p);

  S += dummy_eval;
  dummy_subs += dummy_eval;

  S *= h / 6.;
  dummy_subs *= h / 12.;

  double m = 0;
  while (rel >= e) {
    // analytically derived formula
    S1 = 1. / 2. * (S - dummy_subs);
    m = a + h / (2. * K);
    dummy_eval = 0.;
    //  for (double i = a+h/(2.*K); i <= b-h/(2.*K); i += h/K)
    for (int i = 0; i <= N * K - 1; i++) {                                  // m = (a + h/K(1/2+i)); N*K -1
      dummy_eval += (*f)(m, p);
      m += h / K;
    }

    dummy_eval = (2. * h / (3. * K)) * dummy_eval;
    S1 += dummy_eval;                     // analytically derived
    dummy_subs = 1. / 2. * dummy_eval;

    // normal stepsize halfing process

    rel = fabs(S1 - S) / fabs(S1);
    S = S1;
    K *= 2.;
  }

  return S;

}



int main() {

  double x;

  double p[2] = {0., 1.}; // array with mu and sigma

  double q[2] = {0.5, 0.}; // order of singularity

  double result;
  printf("Numerical Integration Programm started ... \n");
  printf("All results will be shown till 15^⁻10.\n");
  printf("1) Integrate the Function x cos(2 pi x²)² \n");
  printf("------------------------------------------\n");

  printf("Riemann Leftsum :\n");
  result = int_left_riemann(0.,2.,NULL,somecos);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Riemann Rightsum :\n");
  result = int_trapezoidal_int(0.,2.,NULL,somecos);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Trapezoidal Rule :\n");
  result = int_trapezoidal_int(0.,2.,NULL,somecos);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Trapezoidal Rule with semiadaptive stepsizes :\n");
  printf("Relative error e = 0.00001\n");
  result = adapt_step_trap(0.,2.,NULL,somecos,0.00001);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Simpsons Rule :\n");
  result = int_simpson_one_loop(0.,2.,NULL,somecos);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Midpoint Rule with semiadaptive stepsizes :\n");
  printf("Relative error e = 0.00001\n");
  result = adapt_step_mid(0.,2.,NULL,somecos, 0.00001, identity_trafo);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Open boundary Integration :\n");
  printf("We integrate the function exp(x²) from 0 to infinity.\n");
  printf("Relative error e = 0.00001\n");
  result = infty_bound(0.,1,NULL,quadexp,0.00001);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Monte Carlo Integration :\n");
  printf("The absoloute error abse = 0.001\n");
  printf("This Integration is O(N²)\n");
  result = montecarlo(0.,2.,NULL,somecos,0.001);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("2) Integrate the Function the gaussian with mean = 0 and standart deviation = 1 \n");
  printf("------------------------------------------\n");

  printf("Riemann Leftsum :\n");
  result = int_left_riemann(-1.,1.,p,gaussian);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Riemann Rightsum :\n");
  result = int_trapezoidal_int(-1.,1.,p,gaussian);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Trapezoidal Rule :\n");
  result = int_trapezoidal_int(-1.,1.,p,gaussian);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Trapezoidal Rule with semiadaptive stepsizes :\n");
  printf("Relative error e = 0.00001\n");
  result = adapt_step_trap(-1.,1.,p,gaussian,0.00001);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Simpsons Rule :\n");
  result = int_simpson_one_loop(-1.,1.,p,gaussian);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Midpoint Rule with semiadaptive stepsizes :\n");
  printf("Relative error e = 0.00001\n");
  result = adapt_step_mid(-1.,1.,p,gaussian, 0.00001, identity_trafo);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Open boundary Integration : skipped, allready done.\n");


  printf("------------------------------------------\n");

  printf("Monte Carlo Integration :\n");
  printf("The absoloute error abse = 0.001\n");
  printf("This Integration is O(N²)\n");
  result = montecarlo(-1.,1.,p,gaussian,0.001);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");
  printf("------------------------------------------\n");

  printf("Optional Tasks :\n");

  printf("Integrate functions with integrateable Sin:\n");
  printf("Relative error e = 0.00001\n");
  result = sing_int(0., 1.,q,inverse_sqrt,0.00001);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  printf("Simpsons Rule with semiadaptive stepsizes:\n");
  printf("We integrate the cos function from Part 1.");
  printf("Relative error e = 0.00001\n");
  result = adapt_step_simp(0., 2.,NULL,somecos,0.00001);
  printf("The result is : %+2.15lf\n", result);

  printf("------------------------------------------\n");

  return 0;
}

