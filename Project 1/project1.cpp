#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "tridiag_solver.h"

double f(double x);

double u_theory(double x);

bool valid_args(int argc, char *argv[]);


int main( int argc, char *argv[] )
{
  if (!valid_args(argc, argv))
  {
    std::cout << "Usage: " << argv[0] << " N filename <type>" << std::endl;
    std::cout << "With N > 1, and <type> as M for matrix, S for sparse matrix\n" <<
     "and T for Thomas algorithm, L for LU-algorithm" << std::endl;
    exit(1);
  }
  
  const int N = atoi(argv[1]);
  char alg_type = argv[3][0];
  
  Writer fileprinter(argv[2]);
  
  fileprinter.print("N", N);
  
  // Analytical solution
  arma::Col<double> u = f_column(0, 1, N, u_theory);
  // Prints the analytical to the file
  fileprinter.print(u);
  
  
  
  std::clock_t t0, t1;
  double err;
  
  if (alg_type=='M')
  {
    t0 = std::clock();
    arma::Col<double> v = matrix_alg(0, 1, N, f);
    t1 = std::clock();
    
    err = std::log10(max_relative_error(v, u));
    fileprinter.print(v);
  } else if (alg_type == 'S')
  {
    t0 = std::clock();
    arma::Col<double> v = matrix_alg(0, 1, N, f, true);
    t1 = std::clock();
    
    err = std::log10(max_relative_error(v, u));
    fileprinter.print(v);
  } else if (alg_type == 'T')
  {
    t0 = std::clock();
    double *v = new double[N];
    thomas_alg(v, 0, 1, N, f);
    t1 = std::clock();
    
    fileprinter.print(N, v);
    err = std::log10(max_relative_error(v, u));
    delete[] v;
  } else //if (alg_type == L)
  {
    t0 = std::clock();
    arma::Col<double> v = LU_alg(0, 1, N, f);
    t1 = std::clock();
    err = std::log10(max_relative_error(v, u));
    fileprinter.print(v);
  }
  

  double time_used = (t1-t0)/(double)CLOCKS_PER_SEC;

  fileprinter.newline();
  fileprinter.print("Time used", time_used);
  fileprinter.print("max(log10(rel_error))", err);
  
  return 0;
}

double f(double x)
{
  return 100.0*std::exp(-10.0*x);
}

double u_theory(double x)
{
  return 1.0 - (1.0 - std::exp(-10.0))*x - std::exp(-10.0*x);
}

bool valid_args(int argc, char *argv[])
{
  if (argc < 4) {
    return false;
  }
  if (atoi(argv[1]) < 2) {
    return false;
  } 
  char alg_type = argv[3][0];
  return (alg_type == 'M' || alg_type == 'S' || alg_type == 'T'
            || alg_type == 'L');
}