#ifndef lattice_hpp
#define lattice_hpp

#include <iostream>
#include <cstdint>

typedef std::int_fast8_t lat_t; // lattice type

//! A lattice which contains spins for the Ising model.
/*!
 * This class is used in the Ising model, and so has functions
 * which computes the energy for a single lattice point, and 
 * total energy in the entire system along with the total magnetisation
 */
class Lattice
{
private:
  //! Array to hold all the spin values
  lat_t *lattice;
  
  //! Dimension along x-axis (horisontal size)
  const int Lx;
  
  //! Dimension along y-axis (vertical size)
  const int Ly;
  
  //! Function pointer to a print formatter
  std::string (*print_format)(lat_t);
  
public:
  //! Does not initialise any values of the lattice
  /*
   * Default print format is print_t::numbers
   */
  Lattice(const int lx, const int ly);
  
  //! Initialises with the function specified
  /*
   * The init function sets each element with the input
   * given by the indexes of that element
   * 
   * Default print format is print_t::numbers
   */
  Lattice(const int lx, const int ly, lat_t (*init)(int, int));
  
  //! Destructor
  ~Lattice();
  
  
  //! Returns the size of the lattice
  void get_size(int &Lx, int&Ly) const;
  
  //! Returns a reference to the element x, y of the lattice
  /*!
   * This function is periodic, so element x + Lx = x, and
   * x - Lx = x for both dimensions
   */
  lat_t &operator() (int x, int y);
  
  //! Returns the element x,y of the lattice
  lat_t operator() (int x, int y) const;
  
  //! Changes the print format used
  void set_print_format(std::string (*print_func)(lat_t));
  
  //! Prints the lattice
  /*! 
   * The lattice is printed with x along the horizontal axis
   * and y along the vertical. The format it is printed in is 
   * decided by print_format, which takes the number in element (x,y)
   * and transforms it into the suitable format. 
   * 
   * The row-elements (x) are printed horizontally, and a newline is 
   * appended to the y-elements to the next row.
   */
  friend std::ostream& operator<< (std::ostream &out, const Lattice &lat);
  
  //! Sums all the elements of the lattice
  int sum_spins() const;
  

  //! Calculates the energy for a single spin surrounded by the four neighbours
  int energy(const int x, const int y) const;
  
  //! Sums all the individual energies for all the lattice points
  int energy() const;
  
  //! Gives the change in energy if the spin in (x,y) were to be flipped
  int dE(const int x, const int y) const;
  
  //! Returns the pointer to the lattice
  lat_t *buffer();
};

//! This namespace contains print formatters for the Lattice class
namespace print_t {
  //! Arrows with up for 1, and down for -1
  std::string arrows(lat_t num);
  //! Space separated numbers
  std::string numbers(lat_t num);
  //! Try for yourself!
  std::string crazy(lat_t num);
}

//! Contains initialisers for the Lattice class
namespace init {
  //! Sets all spins to spin up (+1). 
  /*!
   * This is equal to a "cold" start, as this gives the 
   * minimal energy
   */
  lat_t up(int, int);
  
  //! Sets all spins randomly to +1 or -1
  /*!
   * This is equal to a "warm" start, as this almost gives the
   * maximal energy
   */
  lat_t random(int, int);
}

#endif