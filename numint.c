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
  double h=(b-a)/N;

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

  // trapezoidal rule
  double T=(*f)(a,p)+(*f)(b,p);
  for (double i=a+h; i<=b-h; i+=h){
	T+=2*(*f)(i, p);
  }
  printf("Trapezoidal int. of function is %f\n", h/2*T); 

  // Simpson's rule
  double S=(*f)(a,p)+(*f)(b,p);
  for (double i=a+h; i<=b-h; i+=h){
  	S+=2*(*f)(i, p);
  }
  for (double i=a+h/2; i<=b-h/2; i+=h){
	S+=4*(*f)(i, p);
  }
  printf("Simpson's rule gives us %f\n", h/6*S); 

}

 
int main(){
  
  double x;
  double p[2]={0., 1.}; // array with mu and sigma
  gaussian(x, p);
  integration_1(-1., 1., p, gaussian);
  
}

