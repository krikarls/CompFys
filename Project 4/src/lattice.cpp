#include <iostream>
#include <sstream>
#include <random>
#include <chrono>
#include <functional>

#include "lattice.hpp"

// The following functions makes the compiler compile so the
// most (un)likely gets a preference, used in periodic boundary.
// This may not be supported in the architecture, and might be superfluous
// anyway as the CPU has autobranching itself

#define likely(x)      __builtin_expect((x), 1)
#define unlikely(x)    __builtin_expect((x), 0)

Lattice::Lattice(const int lx, const int ly) : Lx(lx), Ly(ly)
{
  lattice = new lat_t[Ly*Lx];

  print_format = print_t::numbers;
}

Lattice::Lattice(const int lx, const int ly, lat_t (*init)(int, int)) : Lx(lx), Ly(ly)
{
  lattice = new lat_t[Ly*Lx];

  for (int yy=0; yy<Ly; yy++){
    for (int xx=0; xx<Lx; xx++){
      lattice[Lx*xx + yy] = init(xx, yy);
    }
  }
  print_format = print_t::numbers;
}

Lattice::~Lattice()
{
  delete[] lattice;
}

void Lattice::get_size(int& Lx, int& Ly) const
{
  Lx = this->Lx;
  Ly = this->Ly;
}


lat_t Lattice::operator()(int x, int y) const
{
  if (unlikely(x >= Lx)) x = x%Lx;
  if (unlikely(x < 0)) x = x%Lx + Lx;
  
  if (unlikely(y >= Ly)) y = y%Ly;
  if (unlikely(y < 0)) y = y%Ly + Ly;
  
  return lattice[Lx*x + y];
}

lat_t& Lattice::operator()(int x, int y)
{
  if (unlikely(x >= Lx)) x = x%Lx;
  if (unlikely(x < 0)) x = x%Lx + Lx;
  
  if (unlikely(y >= Ly)) y = y%Ly;
  if (unlikely(y < 0)) y = y%Ly + Ly;
  
  return lattice[Lx*x + y];
}

lat_t* Lattice::buffer()
{
  return lattice;
}


std::ostream& operator<<(std::ostream& out, const Lattice &lat)
{
  for (int ii=0; ii<lat.Ly; ii++){
    for (int jj=0; jj<lat.Lx; jj++){
      out << lat.print_format(lat(jj,ii));
    }
    out << "\n";
  }
  return out;
}


void Lattice::set_print_format(std::string (*print_func)(lat_t))
{
  print_format = print_func;
}

std::string print_t::numbers(lat_t num)
{
  std::ostringstream stringy;
  stringy << (int)num << " ";
  return stringy.str();
}

std::string print_t::arrows(lat_t num)
{
  if (num > 0) {
    return std::string("▲");
  } else if (num < 0) {
    return std::string("▼");
  } else {
    return std::string("?");
  }
}

std::string print_t::crazy(lat_t num)
{
  std::ostringstream stringy;
  if (num>0) {
    return std::string("  ໒( • ͜ʖ • )७  ");
  } else {
    return std::string(" ლ( ◕ 益 ◕ ) ლ ");
  }
}



lat_t init::up(int, int)
{
  return 1;
}

namespace init { namespace rand {
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  const std::uniform_int_distribution<lat_t> distribution(0,1);
  auto rand = std::bind(distribution, generator);
}}

lat_t init::random(int, int)
{  
  if (init::rand::rand() == 0) {
    return -1;
  } else {
    return 1;
  }
}


int Lattice::sum_spins() const
{
  int sum = 0;
  for (int y=0; y<Ly; y++) {
    for (int x=0; x<Lx; x++) {
      sum += this->operator()(x, y);
    }
  }
  return sum;
}

int Lattice::energy() const
{
  int sum_E = 0;
  for (int y = 0; y<Ly; y++) {
    for (int x = 0; x<Lx; x++) {
      sum_E += this->operator()(x,y)*(this->operator()(x+1,y) + this->operator()(x,y+1));
    }
  }
  return -sum_E;
}

int Lattice::energy(const int x, const int y) const
{
  const int S0 = this->operator()(x, y);
  const int S1 = this->operator()(x+1, y);
  const int S2 = this->operator()(x-1, y);
  const int S3 = this->operator()(x, y+1);
  const int S4 = this->operator()(x, y-1);
  
  const int dE = S0*(S1 + S2 + S3 + S4);
  return -dE;
}

int Lattice::dE(const int x, const int y) const
{
  return -2*this->energy(x, y);
}

