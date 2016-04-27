#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <map>

#include <omp.h>

#include "lattice.hpp"
#include "ising.hpp"

//! Makes uniformly distributed points from T0 to T1 inclusive
void linspace(const double T0, const double T1, const int N, double* T);

void print_to_file(const double *T, const double *E, const double* Cv, const double* M, const double* chi, const double* ar,
                   const int N, const char* filename = "test.txt");

void print_to_file(const int N, const int *Eo, const int *Mo, const int* ao,
                   const int *Er, const int *Mr, const int* ar,
                   const char* filename="test.txt");


int main(const int argc, const char **argv)
{
  if (argc<2 || !std::strcmp(argv[1],"h")){
    std::cout << "Usage of " << argv[0] << ":\n"
              << "All the problems in the exercise set is solved by the alphabetical "
              << "numbering given in the project problem set.\n"
              << "input : output\n"
              << "if ouput is 'test.txt' the file test.txt contains the output" << std::endl;
    std::cout << argv[0] << " b <T> : <E> <sigmaE> <M> <sigmaM>\n"
              << argv[0] << " c <T> : file 'test.txt'\n"
              << argv[0] << " d <T> <Nmeasurements> : pairs of E and the"
              << " number of occurences\n"
              << argv[0] << " e <dim> : file 'test,txt'" << std::endl;
              
  }
  
  const auto global_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

  
  
  //Problem b)
  if (argc>2 && !std::strcmp(argv[1], "b")){
    
    const double T = std::atof(argv[2]);
    const double beta = 1.0/T;    
    const int dim = 2;
    Ising model(dim, global_seed, beta, 'u');
    

    
    double E, M, sigmaE, sigmaM, ar;
    
    model.thermalise(50*dim*dim);
    
    model.find_statistics(10*dim*dim, 1000, E, sigmaE, M, sigmaM, ar);
    
    
    
    std::cout << E << " " << M << " " << cv(beta, 4, sigmaE) << " "
              << chi(beta, 4, sigmaM) << std::endl;
    
  }
  
  //Problem c) energy and magnetisation reaching steady state
  if (argc>2 && !std::strcmp(argv[1],"c")){
    
    const double T = std::atof(argv[2]);
    const int L = 20;
    
    const int Ncycles = 100000;
    int *Er = new int[Ncycles];
    int *Mr = new int[Ncycles];
    int *ar = new int[Ncycles+1];
    ar[0] = 0;
    
    
    int *Eo = new int[Ncycles];
    int *Mo = new int[Ncycles];
    int *ao = new int[Ncycles+1];
    ao[0] = 0;
    
    
    Ising random(L, global_seed, 1/T, 'r');
    Ising ordered(L, global_seed+1, 1/T, 'u');
    
    
    for (int ii=0; ii<Ncycles; ii++){
      Er[ii] = random.get_energy();
      Mr[ii] = random.get_magnetisation();
      const int accept = random.try_flip();
      ar[ii+1] = (accept==FlipCodes::NOT_FLIPPED) ? ar[ii] : ar[ii] + 1;
    }
    
    for (int ii=0; ii<Ncycles; ii++){
      Eo[ii] = ordered.get_energy();
      Mo[ii] = ordered.get_magnetisation();
      const int accept = ordered.try_flip();
      ao[ii+1] = (accept==FlipCodes::NOT_FLIPPED) ? ao[ii] : ao[ii] + 1;
    }
    
    print_to_file(Ncycles, Eo, Mo, ao, Er, Mr, ar);
    
    delete[] Er;
    delete[] Mr;
    delete[] Eo;
    delete[] Mo;
    
  }
  
  
  //Problem d) probability for energy state
  if (argc>3 && !std::strcmp(argv[1],"d")){
    const double T = std::atof(argv[2]);
    const int L = 20;
    
    const int M = std::atoi(argv[3]);
    const int Ntherm = 5000*L*L;
    
    Ising random(L, global_seed, 1.0/T, 'u');
    
    random.thermalise(Ntherm);
    
    std::map<int,int> E;
    
    for (int ii=0; ii<M; ii++){
      E[random.get_energy()] += 1;
      random.thermalise(L*L);
    }
    for (auto elem : E){
      std::cout << elem.first << "," << elem.second << " ";
    }
    std::cout << std::endl;
  }
  
  //Problem e) study of critical temperature for differently sized matrices
  if (argc>2 && !std::strcmp(argv[1], "e")){
    const int dim = std::atoi(argv[2]);
    
    const int Ntemps = 15;
    
    
    const int tau = dim*dim; // a MC-step is one flip per site
    const int Nmeasurements = 100000; // 10000 measurments for each average
  

  
    double *T, *E, *cv, *M, *chi, *ar;
    T = new double[Ntemps];
    E = new double[Ntemps];
    cv = new double[Ntemps];
    M = new double[Ntemps];
    chi = new double[Ntemps];
    ar = new double[Ntemps];
    
    linspace(2.0, 2.4, Ntemps, T);

  
    #pragma omp parallel for default(shared)
    for (int ii = 0; ii<Ntemps; ii++){
      const auto seed = global_seed + ii;
      Ising model(dim, seed, 1.0/T[ii], 'u');
      model.thermalise(100000*tau);
      
      model.find_statistics(tau, Nmeasurements, E[ii], cv[ii], M[ii], chi[ii], ar[ii]);
    }
    
    #pragma omp parallel for default(shared)
    for (int ii=0; ii<Ntemps; ii++){ // modifying cv and chi to be dependent on beta and N
      cv[ii] *= dim*dim/(T[ii]*T[ii]);
      chi[ii] *= (dim*dim)/(T[ii]*T[ii]);
    }
      
    
    print_to_file(T, E, cv, M, chi, ar, Ntemps);
    
    delete[] T;
    delete[] E;
    delete[] cv;
    delete[] M;
    delete[] chi;
    delete[] ar;
  }
  
}



void linspace(const double T0, const double T1, const int N, double* T)
{
  const double deltaT = (T1 - T0)/(N-1);
  for (int ii=0; ii<N; ii++){
    T[ii] = T0 + deltaT*ii;
  }
}

std::string int_array_to_string(const int* int_array, const int size_of_array) {
  std::ostringstream stringy;
  for (int temp = 0; temp < size_of_array; temp++)
    stringy << int_array[temp] << " ";
  return stringy.str();
}

void print_to_file(const int N, const int* Eo, const int* Mo, const int* ao,
                   const int* Er, const int* Mr, const int* ar, const char* filename)
{
  std::ofstream outf(filename);
  outf << int_array_to_string(Eo, N) << "\n";
  outf << int_array_to_string(Mo, N) << "\n";
  outf << int_array_to_string(ao, N) << "\n";
  
  outf << int_array_to_string(Er, N) << "\n";
  outf << int_array_to_string(Mr, N) << "\n";
  outf << int_array_to_string(ar, N) << "\n";
}


void print_to_file(const double *T, const double* E, const double* Cv, const double* M, const double* chi, const double* ar,
                   const int N, const char* filename)
{
  std::ofstream outf(filename);
  outf << N << " N" << "\n\n\n\n\n\n\n\n\n";
  outf << "T E M Cv chi acceptance_rate\n";
  for (int ii=0; ii<N; ii++){
    outf << T[ii] << " " << E[ii] << " " << M[ii] << " " << Cv[ii] << " " << chi[ii] << " " << ar[ii] << "\n";
  }
}


