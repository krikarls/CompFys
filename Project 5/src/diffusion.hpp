#ifndef integrate_hpp
#define integrate_hpp

#include "vector.hpp"

/*! \file diffusion.hpp
 * 
 * @brief Functions to solve the diffusion problem
 * 
 * These functions solves the diffusion problem with schemes based
 * on tridiagonal matrix problems and MC-schemes
 * 
 * 
 */

/*!
 * @brief Contains functions to solve the diffusion equation
 * 
 * The Euler and Crank-Nocolson methods use matrices to compute the next
 * derivative, and there is also a MC method
 * 
 * The boundary conditions need to be 0, so init_vec needs to be the 
 * negative steady solution, as shown in the report
 * 
 * The MC method will advance the system according to the number of MC-steps and 
 * the length of the interval. These methods require a seed
 * 
 */
namespace diffusion{

  /*!
   * @brief Does forward euler steps on the initial vector
   * 
   * @param init_vec The initial vector with ends
   * @param alpha \f$ \frac{\Delta t}{\Delta x^2} \f$
   * @param steps Number of forward Euler moves
   * 
   * @return Returns a vector where 'steps' euler steps has been completed, with ends equal to 0
   * 
   * 
   * The internal method is based on a multiplication
   * with a sparse matrix
   */
  Vector<double> forward_euler(const Vector<double> &init_vec,
                              const double alpha, const int steps);


  /*!
   * @brief Does backward euler steps on the initial vector
   * 
   * @param init_vec The initial vector with ends
   * @param alpha \f$ \frac{\Delta t}{\Delta x^2} \f$
   * @param steps Number of backward Euler moves
   * 
   * @return Returns a vector where 'steps' euler steps has been completed.
   * 
   * The internal method is based on an inverse multiplication
   * with a sparse matrix
   */
  Vector<double> backward_euler(const Vector<double> &init_vec,
                                const double alpha, const int steps);

  /*!
  * @brief Evolves the system using the Crank Nicolson scheme
  * 
  * @param init_vec The initial input vector with ends
  * @param alpha \f$ \frac{\Delta t}{\Delta x^2} \f$
  * @param steps Number of Crank Nicolson steps
  * 
  * @return A vector with the solution after the number of
  * steps specified
  * 
  * The internal method is based on an inverse multiplication
  * with a sparse matrix and a multiplication with sparse matrix
  */
  Vector<double> Crank_Nicolson(const Vector<double> &init_vec,
                                const double alpha, const int steps);
  
  /*!
   * @brief Monte Carlo simulation of diffusion equation
   * 
   * @param steps number of MC-steps
   * @param fill_rigth The amount of particles which should be ''refilled''
   * on the right side after one MC-step.
   * @param seed The seed used to initialize the random number generator
   * 
   * @return A vector with the elements as the number of particles currently 
   * on this lattice point
   */
  Vector<int> Monte_Carlo(const int N, const int steps, const int fill_rigth,
                          const long seed);
  
  /*!
   * @brief Monte Carlo simulation of diffusion equation
   * 
   * @param N Amount of points along axis
   * @param steps number of MC-steps
   * @param fill_rigth The amount of particles which should be ''refilled''
   * on the right side after one MC-step.
   * @param seed The seed used to initialize the random number generator
   * 
   * @return A vector with the elements as the number of particles currently 
   * on this lattice point
   */
  Vector<int> Monte_Carlo(const Vector<int> &init_vec,
                          const int steps, const int fill_rigth,
                          const long seed);
  
  /*! 
   * @brief Monte Carlo simulation of diffusion equation using uneven steps
   * 
   * Sets up particles on the left side, and lets them move both ways with the 
   * same probability. The step length is variable, and chosen from a normal 
   * distribution.
   * 
   * 
   * @param steps Number of timesteps
   * @param bins Number of bins the particles should be placed in after simulation
   * @param Nparticles Number of particles simulated, in the system at the same time
   * @param L0 Typical length, dependent upon time as \f$ \sqrt{2D\Delta t} \f$
   * @param seed The seed used to initialize the random number generator
   * 
   * @return A binned distribution of the particles after the random walk
   * 
   */
  Vector<int> Monte_Carlo_gaussian(const int steps, const int bins,
                                   const int Nparticles, const double L0,
                                   const long seed);
  
  /*!
   * @brief Analytical solution of diffusion equation for our special problem
   * 
   * The function solved is x-1, the negative of the steady state solution,
   * with boundary conditions 0
   * 
   * @param t Time after initial start-time
   * @param N Number of points the interval is split into
   * @param order Order which the analytical solution is scaled to. The 
   * default value is 2000
   * 
   */
  Vector<double> Analytical(const double t, const int N, const int order=2000);
  
  
  /*!
   * @brief Computes the discrete L2-norm 
   * @param vector Vector which is compared to the analytical solution
   * @param time Time of simulations (dimensionless)
   * @return The L2-norm error
   */
  double Error(const Vector<double> &vector, const double time);
  
  /*!
   * @brief Computes the discrete L2-norm 
   * @param vector Vector which is compared to the analytical solution
   * @param ANAvector Analytical solution
   * @return The L2-norm error
   */
  double Error(const Vector<double> &vector, const Vector<double> &ANAvector);
  
}
#endif