#include <iostream>
#include <cstdlib>
#include <cmath>

#include <armadillo>

#include "helper_files.hpp"
 
void check_args(int argc, char *argv[], int &N, double &rho_0, double &rho_inf,
		bool &two_electron, double &omega_r, bool &eigv)
{
  
  if (argc > 3)
  {
    N = std::atoi(argv[1]);
    rho_inf = std::atof(argv[2]);
    rho_0 = 0;
    eigv = std::atoi(argv[3]);
    two_electron = false;
  } 
  if (argc > 4)
  {
    two_electron = true;
    omega_r = std::atof(argv[4]);
  }
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <N> <rho_inf> <eigv> <optional: omega_r> \n\n";
    std::cerr << "<N> an integer which gives the size and approximation\n";
    std::cerr << "<rho_inf> max limit of rho to simulate infinity\n";
    std::cerr << "<eigv> as 1 for eigenvectors, and 0 to not compute these\n";
    std::cerr << "<omega_r> given if it is a two_electron problem" << std::endl;
    std::exit(1);
  }
}

arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf,
			     double (*V)(double))
{
  arma::Mat<double> A(N, N, arma::fill::zeros);
  
  double h = (rho_inf-rho_0)/(double)N;
  double inv_h_square = 1/(h*h);
  A.diag(1) -= inv_h_square;
  A.diag(-1) -= inv_h_square;
  
  
  for (int iii=0; iii<N; iii++)
  {
    double rho = rho_0 + iii*h;
    A(iii,iii) = 2*inv_h_square + V(rho);
  }
  
  return A;
}

arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf,
			     bool two_electrons, double omega_r)
{
  arma::Mat<double> A(N, N, arma::fill::zeros);
  
  double h = (rho_inf-rho_0)/(double)N;
  double inv_h_square = 1/(h*h);
  A.diag(1) -= inv_h_square;
  A.diag(-1) -= inv_h_square;
  
  
  if (two_electrons)
  {
    double omega_square = omega_r*omega_r;
    for (int iii=0; iii<N; iii++)
    {
      double rho = rho_0 + iii*h;
      A(iii,iii) = 2*inv_h_square + omega_square*rho*rho + 1.0/rho;
    }
  } else
  {
    for (int iii=0; iii<N; iii++)
    {
      double rho = rho_0 + iii*h;
      A(iii,iii) = 2*inv_h_square + rho*rho;
    }
  }
  
  return A;
}



Energies min_three_diag(const arma::Mat<double> &A)
{
  Energies E;
  
  double *min = E.Energy;
  int *ind = E.indexes;
  
  for (int iii=0; iii<3; iii++)
  {
    min[iii] = arma::datum::inf;
    ind[iii] = -1;
  }
  
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    double num = A(iii,iii);
    if (num<min[2])
    {
      min[2] = num;
      ind[2] = iii;
      if (num < min[1])
      {
	min[2] = min[1];
	ind[2] = ind[1];
	min[1] = num;
	ind[1] = iii;
	if (num < min[0])
	{  
	  min[1] = min[0];
	  ind[1] = ind[0];
	  min[0] = num;
	  ind[0] = iii;
	}
      }
    }
  }
  
  return E;
}


arma::vec get_eigv(const arma::Mat<double> &S, int index)
{
  int N = S.n_cols;
  arma::vec x(N);
  
  double square_sum = 0;
  
  for (int jjj=0; jjj<N; jjj++)
  {
    double S_ij = S(index, jjj);
    x(jjj) = S_ij;
    square_sum += S_ij*S_ij;
  }
  
  double sum = std::sqrt(square_sum);
  x /= sum;
  
  return x;
}

arma::Mat<double> identity(int N)
{
  arma::Mat<double> I(N, N, arma::fill::eye);
  return I;
}
