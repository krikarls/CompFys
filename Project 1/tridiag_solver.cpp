#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "tridiag_solver.h"




arma::Col<double> f_column(double x0, double x1, int N, double (*f)(double))
{
  arma::Col<double> b(N);
  double h = (x1-x0)/(double)N;
  
  for (int iii=0; iii<N; iii++)
  {
    double x = x0 + h*iii;
    b[iii] = f(x);
  }
  return b;
}

arma::Col<double> f_column_(double x0, double x1, int N, double (*f)(double))
{
  arma::Col<double> b(N);
  double h = (x1-x0)/(double)N;
  double neg_h_square = -h*h;
  
  for (int iii=0; iii<N; iii++)
  {
    double x = x0 + h*iii;
    b[iii] = f(x)*neg_h_square;
  }
  return b;
}

arma::Mat<double> second_deriv_matr(int N)
{
  arma::Mat<double> A(N, N);
  A.diag(0) += -2.0;
  A.diag(1) += 1;
  A.diag(-1) += 1;
  return A;
}

arma::SpMat<double> spsecond_deriv_matr(int N)
{
  arma::SpMat<double> A(N, N);
  for (int diag=0; diag < N; diag++)
  {
    A(diag, diag) = -2; //diagonal
    if (diag < N - 1) //upper band
    {
      A(diag, diag+1) = 1;
    }
    if (diag > 0) //lower band
    {
      A(diag, diag-1) = 1;
    }
  }
  return A;
}


arma::Col<double> matrix_alg(double x0, double x1, int N, double (*f)(double), bool SPARSE)
{
  arma::Col<double> f_col = f_column_(x0, x1, N, f);
  
  if (SPARSE)
  {
    arma::SpMat<double> A = spsecond_deriv_matr(N);
    arma::Col<double> v = spsolve(A, f_col);
    return v;
  } else {
    arma::Mat<double> A = second_deriv_matr(N); 
    arma::Col<double> v = solve(A, f_col);
    return v;
  }
}

arma::Col<double> LU_alg(double x0, double x1, int N, double (*f)(double))
{
  arma::Col<double> v;
  arma::Mat<double> X = second_deriv_matr(N);
  arma::Mat<double> L, U, P;
  arma::Col<double> f_col = f_column_(x0, x1, N, f);
  
  lu( L, U, P, X );
  
  arma::Col<double> y = inv(P.t()*L)*f_col;
  v = solve(U, y);
  
  return v;  
}

void thomas_alg(double *v, double x0, double x1, int N, double (*f)(double))
{
  
  double *gamma = new double[N];
  double *beta= new double[N];
  
  double *b_tilde = new double[N];
  
  double h = (x1-x0)/(double)N;
  double h_square = h*h;
  
  
  for (int iii = 0; iii<N; iii++)
  {
    double x = x0 + h*iii;
    b_tilde[iii] = f(x)*h_square;
  }
  
  // the diagonal band does not change
  const double a = -1;
  const double b = 2;
  const double c = -1;
  
  
  beta[0] = b_tilde[0]/b;
  gamma[0] = -c/b;
  
  for (int iii=1; iii<N; iii++)
  {
    beta[iii] = (b_tilde[iii] - a*beta[iii-1])/(b+a*gamma[iii-1]);
    gamma[iii] = -c/(b+a*gamma[iii-1]);
  }
  delete[] b_tilde;
  
  v[N-1] = beta[N-1];
  for (int iii=N-2; iii>0; iii--)
  {
    v[iii] = beta[iii] + gamma[iii]*v[iii+1];
  }
  v[0] = 0;
  delete[] gamma;
  delete[] beta;
}

double max_relative_error(const arma::Col<double> &v, const arma::Col<double> &u)
{
  const int N = v.n_elem;
  double max_error = 0;
  for (int iii=0; iii<N; iii++)
  {
    double err = std::abs((v[iii] - u[iii])/u[iii]);
    if (err > max_error && u[iii]) // Checking not null element
    {
      max_error = err;
    }
  }
  return max_error;
}

double max_relative_error(const double* v, const arma::Col<double>& u)
{
  const int N = u.n_elem;
  double max_error = 0;
  for (int iii=0; iii<N; iii++)
  {
    double err = std::abs((v[iii] - u[iii])/u[iii]);
    if (err > max_error && u[iii]) // Checking not null element
    {
      max_error = err;
    }
  }
  return max_error;
}


Writer::Writer(const char* name) : outf(name) 
{
  outf.precision(15);
}

void Writer::print(const char *text, double value)
{
  outf << value << "\t" << text << std::endl;
}

void Writer::print(const char *text, int value)
{
  outf << value << "\t" << text << std::endl;
}

void Writer::print(const arma::Col<double> &vec)
{
  outf << "\n";
  vec.raw_print(outf);
}

void Writer::print(int N, const double *vec)
{
  outf << "\n";
  for (int iii=0; iii<N; iii++)
  {
    outf << vec[iii] << "\n";
  }
}

void Writer::newline()
{
  outf << "\n";
}
