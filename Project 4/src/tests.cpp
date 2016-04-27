#include "catch.hpp"

#include <cmath>


#include "lattice.hpp"
#include "ising.hpp"

/*! \file tests.cpp
 * \brief Contains unit tests of lattice.hpp and ising.hpp
 */

//! Testing if catch makes a TEST_CASE
TEST_CASE( "SANITIY" )
{
  CHECK( 1 == 1 );
  CHECK( 2 == 3-1 );
}

//! Testing basic operations on a lattice
TEST_CASE( "Setting lattice, indexing and changing elements" , "[lattice]")
{
  int max_x = 2;
  int max_y = 4;
  Lattice lattice(max_x, max_y, init::up);
  
  SECTION( "Check successful init") {
    for (int ii=0; ii<max_y; ii++){
      for (int jj=0; jj<max_x; jj++){
        REQUIRE( lattice(ii,jj) == 1 );
      }
    }
  }
  
  SECTION( "Change element by index") {
    lattice(0, 1) = -5;
    CHECK( lattice(0,1) == -5 );
  }
  
  SECTION( "Indexing periodic boundary" ) {
    lattice(1, 0) = -5;
    lattice(1, max_y - 1) = 3; // last y-element
    
    SECTION( "Wrapping around x-axis" ){
      CHECK(lattice(1+max_x, 0) == -5);
    }
    SECTION( "Wrapping around y-axis" ){
      CHECK(lattice(1, max_y) == -5);
    }
    SECTION( "Wrapping around negative values" ){
      CHECK(lattice(1, -1) == 3); // wrap around negative
    }
  }
}

//! Testing energy and magnetism of a lattice
TEST_CASE( "Energies and spins in lattice", "[lattice]" )
{
  Lattice lattice(2, 2, init::up);
  Lattice lat(5, 5, init::random);
  
  SECTION( "Energy of 2x2 case" ){
    CHECK( lattice.energy() == -8 );
    lattice(1,0) *= -1;
    CHECK( lattice.energy() == 0 );
  }
  
  SECTION( "Check total spins" ) {
    CHECK( lattice.sum_spins() == 4 );
    lattice(1,0) = -1;
    CHECK( lattice.sum_spins() == 2 );
  }
  
  SECTION( "Change in energy by flipping single spin given by dE" ) {
    for (int ii=0; ii<5; ii++){
      const int deltaE = lat.dE(1,ii);
      
      const int Ebefore = lat.energy();
      lat(1,ii) *= -1; //flipping
      const int Eafter = lat.energy();
      
      CHECK( deltaE == Eafter - Ebefore );
      
      CHECK( lat.dE(1,ii) == Ebefore - Eafter ); // Reverse flip
    }
  }
  
  SECTION( "dE should give numbers in the range {-8,-4,0,4,8}"){
    for (int ii=0; ii<5; ii++){
      for (int jj=0; jj<5; jj++){
        const int dE = lat.dE(ii,jj);
        const bool cheq = ((dE == -8) || (dE == -4) || (dE == 0) || (dE == 4) || (dE == 8));
        CHECK( cheq );
      }
    }
    
  }
}
    
//! Testing the Ising class
TEST_CASE( "Initialising and testing energies and mag during flips", "[ising]" )
{
  const int seed = 4; // Random seed, decided by throw of a dice
  const int L = 5;
  const double beta = 0;
  Ising model(L, seed, beta);
  model.init_rand();
  
  SECTION( "Energy correctly init" ) {
    CHECK( model.get_energy() == model.recompute_energy() );
  }
  
  SECTION( "Is the energy correctly accumulated" ) {
    model.try_flip(100);
    const int Estored = model.get_energy();
    const int Ecomp = model.recompute_energy();
    
    CHECK( Estored == Ecomp );
    
  }
  
  SECTION( "Magnetisation correctly init" ) {
    CHECK( model.get_magnetisation() == model.recompute_magnetisation() );
  }
  
  SECTION( "Is the magnetisation correctly accumulated" ) {
    model.try_flip(100);
    const int Mstored = model.get_magnetisation();
    const int Mcomp = model.recompute_magnetisation();
    
    CHECK( Mstored == Mcomp );
    
  }    
    
  
}