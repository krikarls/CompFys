#ifndef tridiag_H
#define tridiag_H

#include <armadillo>
#include <fstream>


arma::SpMat<double> spsecond_deriv_matr(int N);

arma::Mat<double> second_deriv_matr(int N);

arma::Col<double> f_column(double x0, double x1, int N, double (*f)(double));

arma::Col<double> f_column_(double x0, double x1, int N, double (*f)(double));

void thomas_alg(double *v, double x0, double x1, int N, double (*f)(double));

arma::Col<double> matrix_alg(double x0, double x1, int N, double (*f)(double), bool SPARSE=false);

arma::Col<double> LU_alg(double x0, double x1, int N, double (*f)(double));

double max_relative_error(const arma::Col<double> &v, const arma::Col<double> &u);

double max_relative_error(const double* v, const arma::Col<double> &u);

class Writer
{
private:
  std::ofstream outf;
public:
  Writer(const char* name);
  void print(const char *text, const double value);
  void print(const char* text, int value);
  void print(const arma::Col<double> &vec);
  void print(int N, const double *vec);
  void newline();
};

#endif
