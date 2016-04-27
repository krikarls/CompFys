#include <cmath>

#include <armadillo>

#include "Jacobi_rotation.hpp"
#include "unittest++/UnitTest++.h"



double abs_sum_offdiag(const arma::Mat<double> &A) // Can be used for tolerance check
{
  double sum = 0;
  for (int iii=0; iii < (int)A.n_rows; iii++)
  {
    for (int jjj=0; jjj < (int)A.n_cols; jjj++)
     {
       if(iii != jjj){
         sum += std::abs(A(iii,jjj));
       }
     }
  }
  return sum;
}

void find_cos_sin(double a_kk, double a_ll, double a_kl, double &c, double &s)
{
  double tau = (a_ll - a_kk)/(2*a_kl);
  double t = tau <= 0 ? -tau + std::sqrt(1+tau*tau) : -tau -std::sqrt(1+tau*tau);
  
  c = 1/std::sqrt(1+t*t);
  s = t*c;  
}

void rotate(arma::Mat<double> &A, double c, double s, int k, int l)
{
  double a_kk = A(k, k);
  double a_ll = A(l, l);
  double a_kl = A(k, l);
  
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    double a_ik = A(iii, k);
    A(iii, k) = A(iii, k)*c - A(iii, l)*s;
    A(k, iii) = A(iii, k);
    A(iii, l) = A(iii, l)*c + a_ik*s;
    A(l, iii) = A(iii, l);
  }
  
  A(k, k) = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s;
  A(l, l) = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s;
  A(k,l) = 0; // (a_kk-a_ll)*c*s + a_kl*(c*c-s*s); // DEBUG
  A(l,k) = 0; // (a_kk-a_ll)*c*s + a_kl*(c*c-s*s); // DEBUG
}

void rotate_with_eigvec(arma::Mat<double> &A, arma::Mat<double> &S,
			double c, double s, int k, int l)
{
  
  int N = A.n_rows;
  double a_kk = A(k, k);
  double a_ll = A(l, l);
  double a_kl = A(k, l);
  
  for (int iii=0; iii<N; iii++)
  {
    double a_ik = A(iii, k);
    A(iii, k) = A(iii, k)*c - A(iii, l)*s;
    A(k, iii) = A(iii, k);
    A(iii, l) = A(iii, l)*c + a_ik*s;
    A(l, iii) = A(iii, l);
  }
  
  A(k, k) = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s;
  A(l, l) = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s;
  A(k,l) = 0; // (a_kk-a_ll)*c*s + a_kl*(c*c-s*s); // DEBUG
  A(l,k) = 0; // (a_kk-a_ll)*c*s + a_kl*(c*c-s*s); // DEBUG
  
  
  //eigenvector part
  for (int jjj=0; jjj<N; jjj++)
  {
    double s_kj = S(k, jjj);
    double s_lj = S(l, jjj);
    S(k, jjj) = c*s_kj - s*s_lj;
    S(l, jjj) = s*s_kj + c*s_lj;
  }
  
}

void max_err_offdiag(const arma::Mat<double> &A, int &k, int &l, double &err)
{
  err = 0;
  l = -1;
  k = -1;
  
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    for (int jjj=0; jjj<(int)A.n_cols; jjj++)
    {
      if ((jjj != iii) && (std::abs(A(iii,jjj)) > err))
      {
        err = std::abs(A(iii,jjj));
        k = iii;
        l = jjj;
      }
    }
  }
}

int rotate_to_diag(arma::Mat<double> &A, double tolerance)
{
  
  int n_rotations = 0;
  double max_err = tolerance + 1.0;
  
  while (max_err > tolerance)
  {
    int k, l;
    max_err_offdiag(A, k, l, max_err);
    
    double cos, sin;
    find_cos_sin(A(k,k), A(l,l), A(k,l), cos, sin);
    
    rotate(A, cos, sin, k, l);
    n_rotations++;
  }
  return n_rotations;
}

int rotate_to_diag_with_eigvec(arma::Mat<double> &A, arma::Mat<double> &S, double tolerance)
{
  int n_rotations = 0;
  double max_err = tolerance + 1.0;
  
  while (max_err > tolerance)
  {
    int k, l;
    max_err_offdiag(A, k, l, max_err);
    
    double cos, sin;
    find_cos_sin(A(k,k), A(l,l), A(k,l), cos, sin);
    
    rotate_with_eigvec(A, S, cos, sin, k, l);
    n_rotations++;
  }
  return n_rotations;
}