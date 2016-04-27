#include "CImg.h"
using namespace cimg_library;

#include <ctime>
#include <cstdlib>
#include <iostream>

#include "ising.hpp"

/*! \file animate.cpp
 * \brief Program to make an animation of the Ising system
 * 
 * This program sets up a window and prints the Ising system.
 * The program used is CImg, which is available from 
 * https://github.com/dtschump/CImg and http://cimg.eu/.
 * This program is installed as "sudo apt install cimg-dev" on
 * debian systems.
 * 
 * The system is set up randomly, and shows domains when enough
 * cycles has been completed.
 */ 




//! Main program of animate
int main(int argc, char **argv) {
  int L = 100;
  double beta = 1/2.4;
  
  if ((argc > 1) && (argv[1][0] == 'h')){
    std::cout << "This program plots the spins in the Ising model in two dimensions\n"
              << "The default settings are L = " << L << " and beta = " << beta << "\n"
              << "The parameters beta and L can be selected by passing these by the command line\n"
              << "Supported inputs are '" << argv[0] << " beta' and '" << argv[0] << " beta L'" << std::endl;
    std::exit(0);
  }

  if (argc > 1){
    beta = std::atof(argv[1]);
  }
  if (argc > 2){
    L = std::atoi(argv[2]);
  }
  
  
  
  const long int seed = std::clock();
  const int wait_time = 1;// time in milliseconds between redraws
  
  
  
  Ising model(L, seed, beta);
  model.init_rand();
  
  const int calculate_N = L*L/10;
  
  // Getting the future buffer
  const lat_t* image_buffer = model.buffer();

  //making a window
  CImgDisplay window(900,600,"Ising model",0);
  
  //binding a buffer to image
  CImg<lat_t> image(image_buffer,L,L,1,1,true);
  
  
  while (!window.is_closed()) { // event binding
    if (window.is_resized()) {
      window.resize();
    }
    window.display(image);
    
    
    // updating buffer live here
    model.try_flip(calculate_N);
    
    window.wait(wait_time); // limits the framerate
  }
  
    
  
  
}
