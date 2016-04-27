
#include <armadillo>
#include <cmath>

#include "Jacobi_rotation.hpp"
#include "helper_files.hpp"
#include "unittest++/UnitTest++.h"

const double tolerance = 1e-10;

SUITE(Jacobi)
{
  TEST(Rotate)
  {
    int N = 5;
    
    arma::arma_rng::set_seed_random();
    arma::Mat<double> A(N, N, arma::fill::randu);
    
    
    for (int iii=0; iii<N; iii++) // Symmetric matrix
    {
      for (int jjj=0; jjj<N; jjj++)
      {
        if (jjj>iii)
	{
	  A(jjj, iii) = A(iii, jjj);
	}
      }
    }
    
    double cos = 0.5;
    double sin = std::sqrt(1-cos*cos);
    int k = 2;
    int l = 4;
    
    arma::Mat<double> S(N, N, arma::fill::eye);
    S(k, k) = cos;
    S(l, l) = cos;
    
    S(k, l) = sin;
    S(l, k) = -sin;
    
    arma::Mat<double> B = S.t()*A*S;
    rotate(A, cos, sin, k, l);
    
    
    for (int iii=0; iii<N; iii++)
    {
      for (int jjj=0; jjj<N; jjj++)
      {
	if ((iii != k) && (iii != l)) // ignore the hard-coded zero-values
	{
	  CHECK_CLOSE(A(iii,jjj), B(iii,jjj), tolerance);
	}
      }
    }
  }
  TEST(max_index)
  {
    arma::Mat<double> A;
    A << 1 << 0 << 6 << arma::endr
      << 0 << 3 << 9;
    int k = 0;
    int l = 0;
    double err = 0;
    
    max_err_offdiag(A, k, l, err);
    
    CHECK_EQUAL(k, 1);
    CHECK_EQUAL(l, 2);
    CHECK_EQUAL(err, 9);
  }
  TEST(sin_cos)
  {
    
   /* A = [1        off_number]
    *     [off_number        3]
    */
   
    double off_number = 1/std::sqrt(2);
    double tau = (3-1)/((int)2*off_number);
    double t = -tau - std::sqrt(1+ tau*tau);
    double c = 1/std::sqrt(1+t*t);
    double s = t*c;

    
    double x,y;
    find_cos_sin(1, 3, off_number, x, y);
    
    CHECK_CLOSE(x, c, tolerance);
    CHECK_CLOSE(y, s, tolerance);
  }
  TEST(abs_sum_offdiag)
  {
    arma::Mat<double> A;
    
    A << 1 << 6 << -1 << arma::endr
      << 4 << 999 << 9;
    
    CHECK_CLOSE( abs_sum_offdiag(A), 6+1+4+9, tolerance);
    
  }
}


SUITE(helper_files)
{
  TEST(Sort)
  {
    int N = 6;
    arma::Mat<double> A(N, N, arma::fill::zeros);
    A(0,0) = 9;
    A(1,1) = 4;
    A(2,2) = 5;
    A(3,3) = -6;
    A(5,5) = 400;
    
    A(5,3) = -78;
    
    auto E = min_three_diag(A);
    CHECK_CLOSE(E.Energy[0], -6, tolerance);
    CHECK_CLOSE(E.Energy[1], 0, tolerance);
    CHECK_CLOSE(E.Energy[2], 4, tolerance);
    
    CHECK_EQUAL(E.indexes[0], 3);
    CHECK_EQUAL(E.indexes[1], 4);
    CHECK_EQUAL(E.indexes[2], 1);
  }
  
  TEST(Identity)
  {
    int N = 3;
    auto A = identity(N);
    for (int iii=0; iii<N; iii++)
    {
      for (int jjj=0; jjj<N; jjj++)
      {
	if (iii!=jjj)
	{
	  CHECK_CLOSE(A(iii,jjj), 0, tolerance);
	} else
	{
	  CHECK_CLOSE(A(iii,jjj), 1, tolerance);
	}
      }
    }
    
  }

  TEST(Hamiltonian)
  {
    int N = 5;
    const auto A = ham_matrix(N, 0, 10);
    double h = 10/(double)N;
    
    for (int iii=0; iii<N; iii++) // diagonal elements
    {
      double compare = 2/(h*h) + (h*iii)*(h*iii);
      CHECK_CLOSE(compare, A(iii,iii), tolerance);
    }
    
    for (int iii=2; iii<N; iii++) // two of the diagonal
    {
      CHECK_CLOSE(0, A(0,iii), tolerance);
    }
    
    for (int iii=0; iii<N-1; iii++) // one of the diagonal
    {
      CHECK_CLOSE(-1/(h*h), A(iii, iii+1), tolerance);
    }
  }
  
  TEST(Eigenvector)
  {
    int N = 5;
    arma::Mat<double> A(N, N, arma::fill::randu);
    

    
    int index = 3;
    
    auto eigv = get_eigv(A, index);
    
    
    double sum = 0;
    for (int iii=0; iii<N; iii++)
    {
      sum += (A(index, iii))*(A(index, iii));
    }
    A /= std::sqrt(sum);
    
    
    for (int iii=0; iii<N; iii++)
    {
      CHECK_CLOSE(A(index, iii), eigv(iii), tolerance);
    }
    
  }
}