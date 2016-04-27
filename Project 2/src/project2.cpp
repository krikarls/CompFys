#include <iostream>
#include <cstdlib>

#include <armadillo>
#include "Jacobi_rotation.hpp"
#include "filewriter.hpp"
#include "helper_files.hpp"

#include "unittest++/UnitTest++.h"

inline double hydrogen(double r)
{
  return -10/r;  
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    UnitTest::RunAllTests();
    std::cout << std::endl;
  }
  
  int N;
  double rho_0, rho_inf, omega_r;
  bool two_electrons, eigv;
  
  rho_0 = 0;
  check_args(argc, argv, N, rho_0, rho_inf, two_electrons, omega_r, eigv);
  
  
  arma::Mat<double> A = ham_matrix(N, rho_0, rho_inf, two_electrons, omega_r);
  //arma::mat A = ham_matrix(N, rho_0, rho_inf, hydrogen);
  int n_rotations;
  
  if (eigv)
  {
    arma::Mat<double> S = identity(N);
    n_rotations = rotate_to_diag_with_eigvec(A, S, 1e-8);
    std::cout << "Number of rotations = " << n_rotations << std::endl;

    auto E = min_three_diag(A);
    
    
    FileWriter file("test.txt");
    file.print(N, rho_0, rho_inf, two_electrons, omega_r);
    file.print(E);
    
    for (int iii=0; iii<3; iii++)
    {
      auto Evec = get_eigv(S, E.indexes[iii]);
      file.print(Evec);
    }
  } else {
    n_rotations = rotate_to_diag(A, 1e-8);
    std::cout << "Number of rotations = " << n_rotations << std::endl;

    auto E = min_three_diag(A);
    
    
    FileWriter file("test.txt");
    file.print(N, rho_0, rho_inf, two_electrons, omega_r);
    file.print(E);
    
  }
  return 0;
}


