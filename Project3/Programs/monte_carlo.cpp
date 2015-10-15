#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace arma;

double integrand_c(double x_1, double y_1,double z_1,double x_2,double y_2,double z_2);
double integrand_d(double r1, double r2,double costheta1,double costheta2,double phi1,double phi2);

void MC()
{
    int N;
    double a, b, fx,r1,r2,costheta1,costheta2,phi1,phi2, exact,lambda, int_bMC, int_MC, sum_sigma, var_bMC, var_MC;
    const double PI = 3.141592653589793238463;

    rowvec y = zeros<rowvec>(6);    // making space for 6 random numbers in [a,b]
    exact = (5*PI*PI)/(16*16);

    std::cout << "Read in the number Monte Carlo cycles" << std::endl;
    std::cin >> N;
    std::cout << "Read in integration limits" << std::endl;
    std::cin >> a >> b;

//  computation loop: Brute force monte carlo
    fx = 0; sum_sigma= 0;
    for (int n=1; n <= N; n++)
    {
            rowvec x = randu<rowvec>(6);      // six random numbers in [0,1]
            y = a + x*(b-a);                  // transformation to [a,b] (general normal dist. )

            fx += integrand_c(y(0),y(1),y(2),y(3),y(4),y(5));
            sum_sigma += fx*fx;
    }

    int_bMC = (fx/N)*pow((b-a),6);
    sum_sigma = sum_sigma/N;
    var_bMC = sum_sigma-(fx/N)*(fx/N);

    arma_rng::set_seed_random(); // change RNG seed

//  computation loop: Monte carlo using importance sampling
    lambda = 1.0/4;     // exp. weight
    fx = 0;  sum_sigma= 0;
    for (int n=1; n <= N; n++)
    {
            rowvec x = randu<rowvec>(6);    // draw six random numbers in [0,1]

            r1 = -lambda*log(1-x(0));       // map to exp.dist [0,inf]
            r2 = -lambda*log(1-x(1));       // map to exp.dist [0,inf]
            costheta1 = 2*x(2)-1;           // map to uni. dist [-1,1]
            costheta2 = 2*x(3)-1;           // map to uni. dist [-1,1]
            phi1 = x(4)*2*PI;               // map to uni. dist [0,2pi]
            phi2 = x(5)*2*PI;               // map to uni. dist [0,2pi]

            fx += integrand_d(r1,r2,costheta1,costheta2,phi1,phi2);
            sum_sigma += fx*fx;
    }

    int_MC = (fx/N)*2*2*(2*PI)*(2*PI)*lambda*lambda;
    sum_sigma = sum_sigma/N;
    var_MC = sum_sigma-(fx/N)*(fx/N)*lambda*lambda;

    std::cout << "Method:                              " << "    I            " << " Var." << std::endl;
    std::cout << "Brute force Monte carlo:             " << int_bMC<< "       " << var_bMC << std::endl;
    std::cout << "Monte carlo(imprtance sampling):     " << int_MC << "       " << var_MC << std::endl;
    std::cout << "Exact solution:                      " << exact <<  std::endl;
}

double integrand_c(double x_1, double y_1,double z_1,double x_2,double y_2,double z_2)
{
  double alpha = 2;
  double num = exp(-2*alpha*(sqrt(x_1*x_1+y_1*y_1+z_1*z_1)+sqrt(x_2*x_2+y_2*y_2+z_2*z_2)));
  double deno = sqrt(pow((x_1-x_2),2)+pow((y_1-y_2),2)+pow((z_1-z_2),2));
  if(deno < pow(10.,-6.)) // to avoid zero division
    {return 0;}
  else return num/deno;
}
double integrand_d(double r1, double r2,double costheta1,double costheta2,double phi1,double phi2)
{
    double num = r1*r1*r2*r2;
    double B = costheta1*costheta2 + sqrt(1-costheta1*costheta1)*sqrt(1-costheta2*costheta2)*cos(phi1-phi2);
    double deno = r1*r1+r2*r2-2*r1*r2*B;
    if(deno < pow(10.,-6.)) // to avoid zero division
      {return 0;}
    else return num/sqrt(deno);
}

