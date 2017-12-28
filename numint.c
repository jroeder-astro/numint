#include<math.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>


// LIBRARY FOR NUMERICAL INTEGRATION



double gaussian(double x, void *params){

  double *par = (double *)params;
  double mu=par[0];
  double sigma=par[1];

  return 1/(sigma*sqrt(2*M_PI))*exp(-pow(((x-mu)/sigma), 2)/2);
}

double integration_1(double a, double b, void *params, double (*f)(double, void *)){
                                  // HOW DOES ONE PUT THE PARAMETERS IN HERE???   -> this apparently works lol  
  double *p=(double*)params;  

  double N=4;
  double h=(fabs(a)+fabs(b))/N;

  // left Riemann sum
  double L=(*f)(a, p);                    // this bitch is the reason you gotta start with a+h in the for loop
  for (double i=a+h; i<=b-h; i+=h){
  	L+=(*f)(i, p);
  }
  printf("Left sum of function is %f\n", h*L);   

  // right Riemann sum
  double R=(*f)(a+h,p);
  for (double i=a+2*h; i<=b; i+=h){
	R+=(*f)(i, p);
  }
  printf("Right sum of function is %f\n", h*R); 

}


/*   // works just fine!
double f(double x){ return x*x; }         //   1/(sqrt(2*M_PI))*exp(-pow(x,2)/2); }

double integration_1(double a, double b, double (*f)(double)){

  // Riemann left sum
  double N=4;
  double h=(fabs(a)+fabs(b))/N;
  printf("h=%f\n", h);
  double T=(*f)(a);
  for (double i=a+h; i<=b-h; i+=h){
  	T+=(*f)(i);} 
  printf("Left sum of function is %f\n", h*T);
}
*/

int main(){
  
  double x;
  double p[2]={0., 1.}; // array with mu and sigma
  gaussian(x, p);
  // printf("it is %f\n", gaussian(x,p));
  integration_1(-1., 1., p, gaussian);

  // integration_1(0, 2., f);
}

