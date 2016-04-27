#ifndef helper_files_h
#define helper_files_h

#include <armadillo>

struct Energies
{
  double Energy[3];
  int indexes[3];
};

void check_args(int argc, char *argv[], int &N, double &rho_0, double &rho_inf,
		bool &two_electron, double &omega_r, bool &eigv);

Energies min_three_diag(const arma::Mat<double> &A);

arma::vec get_eigv(const arma::Mat<double> &S, int index);

arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf,
			     bool two_electrons=false, double omega_r=1);

arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf,
			     double (*V)(double));

arma::Mat<double> identity(int N);

#endif