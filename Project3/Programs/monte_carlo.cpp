#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
using namespace std;

//     Here we define various functions called by the main program

double integrand_a(double x1, double y1, double z1, double x2, double y2, double z2);
double integrand_b(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double gammln(double);

//   Main function begins here
void GQ()
{
     int N;
     double a, b;

     cout << "Read in the number of integration points" << endl;
     cin >> N;
     cout << "Read in integration limits" << endl;
     cin >> a >> b;

     double *r = new double[N+1];               // radial mesh points [0,inf] from Laguerre polynomials
     double *w_r = new double[N+1];             // corresponding weights
     double *theta = new double[N];             // mesh points for theta-angle [0,pi] from Legendre polynomials
     double *w_theta = new double [N];           // corresponding weights
     double *phi = new double[N];               // mesh points for phi-angle [0,2*pi] from Legendre polynomials
     double *w_phi = new double[N];             // corresponding weights
     double *x = new double[N];                 // mesh points for "brute force" Legendre approach
     double *w_x = new double[N];               // corresponding wights

     double exact = (5*M_PI*M_PI)/(16.0*16.0);

     gauss_laguerre(r, w_r, N, 2);      // filling in mesh and weight vectors
     gauleg(0,M_PI,theta ,w_theta,N);
     gauleg(0,2*M_PI,phi ,w_phi,N);
     gauleg(a,b, x, w_x, N);

//   loops for computations

     double int_gauss_leg = 0.;

     for ( int i = 0;  i < N; i++){     // Legendre loop
         for (int j = 0; j < N; j++){
         for (int k = 0; k < N; k++){
         for (int l = 0; l < N; l++){
         for (int m = 0; m < N; m++){
         for (int n = 0; n < N; n++){

            double weights = w_x[i]*w_x[j]*w_x[k]*w_x[l]*w_x[m]*w_x[n];
            int_gauss_leg += weights*integrand_a(x[i],x[j],x[k],x[l],x[m],x[n]);

          }}}}}} // end Legendre loop

     double int_gauss_lag = 0.;
     double norm_fact = 1.0/1024;   // (2*alpha)^-5

     for ( int i = 1;  i <= N; i++){    // Laguerre loop
         for (int j = 1; j <= N; j++){
         for (int k = 0; k < N; k++){
         for (int l = 0; l < N; l++){
         for (int m = 0; m < N; m++){
         for (int n = 0; n < N; n++){

             double weights = w_r[i]*w_r[j]*w_theta[k]*w_theta[l]*w_phi[m]*w_phi[n];
             int_gauss_lag += weights*integrand_b(r[i],r[j],theta[k],theta[l],phi[m],phi[n]);

          }}}}}} // end Laguerre loop

     int_gauss_lag *= norm_fact;

//   end of computation loops

//    final output
      cout  << setiosflags(ios::showpoint | ios::uppercase);
      cout << "Exact solution = "<< setw(20) << setprecision(15)  << exact << endl;
      cout << "Gaussian-Legendre quad = "<< setw(20) << setprecision(15)  << int_gauss_leg << endl;
      cout << "Gaussian-Laguerre quad = " << setw(20) << setprecision(15) << int_gauss_lag  << endl;
}  // end of main program

//  this function defines the integrand used for task a)
double integrand_a(double x1, double y1, double z1, double x2, double y2, double z2)
{
   double alpha = 2.;
   double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
   double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
   double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));

  if(deno < 1e-6) { return 0;}
  else return exp(exp1+exp2)/deno;
} // end of integrand a function

//  this function defines the integrand used for task b)
double integrand_b(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    double x1 = r1;
    double x2 = r2;
    //double B = cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2);
    //double B = theta1*theta2 + sqrt(1-theta1*theta1)*sqrt(1-theta2*theta2)*cos(phi1-phi2);
    double B = cos(theta1)*cos(theta2)+sqrt(1-cos(theta1)*cos(theta1))*sqrt(1-cos(theta2)*cos(theta2))*cos(phi1-phi2);
    double deno = x1*x1+x2*x2-2*x1*x2*B;

    if(deno < 1e-6) {return 0;}
    else return (sin(theta1)*sin(theta2))/sqrt(deno);
} // end of integrand b function

       /*
       ** The function gauleg() takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauleg(double x1, double x2, double x[], double w[], int N)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (N + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + N - 1;
   w_low  = w;
   w_high = w + N - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(N + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= N; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

     pp = N * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()


void gauss_laguerre(double *x, double *w, int N, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=N;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*N+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*N);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=N;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(N*p1-(N+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+N)-gammln((double)N))/(pp*N*p2);
    }
}
// end function gaulag

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
//#undef EPS
//#undef MAXIT
